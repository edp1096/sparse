package sparse // import "sparse"

import (
	"fmt"
)

func checkConfig(config *Configuration) *Configuration {
	defaultConfig := Configuration{
		Real:                    true,
		Complex:                 true,
		SeparatedComplexVectors: false,
		Expandable:              true,
		Translate:               true,
		ModifiedNodal:           true,
		DefaultThreshold:        1.0e-3,
		TiesMultiplier:          5,
		DefaultPartition:        AUTO_PARTITION,
		PrinterWidth:            80,
		Annotate:                0,
	}

	if config == nil {
		config = &defaultConfig
	}

	if config.DefaultThreshold < 0.0 {
		config.DefaultThreshold = 1.0e-3
	}
	if config.TiesMultiplier <= 0 {
		config.TiesMultiplier = 5
	}
	if config.DefaultPartition == DEFAULT_PARTITION {
		config.DefaultPartition = AUTO_PARTITION
	}
	if config.PrinterWidth <= 0 {
		config.PrinterWidth = 80
	}

	return config
}

func Create(size int64, config *Configuration) (*Matrix, error) {
	if size < 0 || (size == 0 && !config.Expandable) {
		return nil, fmt.Errorf("invalid size: %d", size)
	}

	config = checkConfig(config)

	matrixSize := size + 1 // 1-based indexing

	m := &Matrix{
		Config:          *config,
		Size:            size,
		ExtSize:         size,
		Complex:         config.Complex,
		DoRealDirect:    make([]bool, matrixSize),
		DoComplexDirect: make([]bool, matrixSize),
		CurrentSize:     0,
		Elements:        0,
		Diags:           make([]*Element, matrixSize),
		FirstInRow:      make([]*Element, matrixSize),
		FirstInCol:      make([]*Element, matrixSize),
		Intermediate:    make([]float64, matrixSize),
		MarkowitzRow:    make([]int64, matrixSize),
		MarkowitzCol:    make([]int64, matrixSize),
		MarkowitzProd:   make([]int64, matrixSize),
		IntToExtRowMap:  make([]int64, matrixSize),
		IntToExtColMap:  make([]int64, matrixSize),
		ExtToIntRowMap:  make([]int64, matrixSize),
		ExtToIntColMap:  make([]int64, matrixSize),
		NeedsOrdering:   true,
		RelThreshold:    config.DefaultThreshold,
		AbsThreshold:    0.0,
		// TrashCan:        &Element{NextInRow: nil, NextInCol: nil},
	}

	if err := m.CreateInternalVectors(); err != nil {
		return nil, fmt.Errorf("failed to create internal vectors: %v", err)
	}

	for i := int64(1); i <= size; i++ {
		m.IntToExtRowMap[i] = i
		m.IntToExtColMap[i] = i
		m.ExtToIntRowMap[i] = -1
		m.ExtToIntColMap[i] = -1
	}
	m.ExtToIntRowMap[0] = 0
	m.ExtToIntColMap[0] = 0

	return m, nil
}

func (m *Matrix) GetInitInfo(element *Element) *ComplexNumber {
	return element.InitInfo
}

func (m *Matrix) SetInitInfo(element *Element, value *ComplexNumber) {
	element.InitInfo = &ComplexNumber{Real: value.Real, Imag: value.Imag}
}

func (m *Matrix) Initialize() error {
	for j := int64(1); j <= m.Size; j++ {
		element := m.FirstInCol[j]
		for element != nil {
			if element.InitInfo == nil {
				element.Real = 0.0
				element.Imag = 0.0
			} else {
				element.Real = element.InitInfo.Real
				element.Imag = element.InitInfo.Imag
			}
			element = element.NextInCol
		}
	}

	m.Factored = false
	m.SingularCol = 0
	m.SingularRow = 0

	// m.TrashCan.Real = 0.0
	// m.TrashCan.Imag = 0.0

	return nil
}

func (m *Matrix) Clear() {
	for i := m.Size; i > 0; i-- {
		element := m.FirstInCol[i]
		for element != nil {
			element.Real = 0.0
			if m.Complex {
				element.Imag = 0.0
			}
			element = element.NextInCol
		}
	}

	m.Factored = false
	m.SingularCol = 0
	m.SingularRow = 0

	// m.TrashCan.Real = 0.0
	// m.TrashCan.Imag = 0.0
}

func (m *Matrix) Destroy() {
	m.DoRealDirect = nil
	m.DoComplexDirect = nil

	m.IntToExtColMap = nil
	m.IntToExtRowMap = nil
	m.ExtToIntColMap = nil
	m.ExtToIntRowMap = nil

	m.Diags = nil
	m.FirstInRow = nil
	m.FirstInCol = nil
	m.Intermediate = nil
	m.MarkowitzRow = nil
	m.MarkowitzCol = nil
	m.MarkowitzProd = nil

	m.Elements = 0

	m.Size = 0
	m.Complex = false
	m.NeedsOrdering = false
	m.Partitioned = false
	m.Factored = false
	m.Reordered = false
	m.InternalVectorsAllocated = false

	m.SingularRow = 0
	m.SingularCol = 0

	m.Fillins = 0

	m.PivotsOriginalRow = 0
	m.PivotsOriginalCol = 0
	m.PivotSelectionMethod = 0

	m.Singletons = 0

	// m.TrashCan.Row = 0
	// m.TrashCan.Col = 0
	// m.TrashCan.Real = 0.0
	// m.TrashCan.Imag = 0.0
	// m.TrashCan.NextInRow = nil
	// m.TrashCan.NextInCol = nil
	// m.TrashCan.InitInfo = nil
}

func (m *Matrix) createElement(row, col int64, firstInRow, firstInCol **Element, fillin bool) *Element {
	const (
		LargestShortInteger = 32767
		LargestLongInteger  = 2147483647
	)

	current := *firstInCol
	var prev **Element = firstInCol
	for current != nil && current.Row < row {
		prev = &current.NextInCol
		current = current.NextInCol
	}

	if current != nil && current.Row == row {
		return current
	}

	var element *Element
	if fillin {
		element = &Element{Row: row, Col: col, Real: 0.0, Imag: 0.0}
		m.Fillins++

		// Update Markowitz Counts
		m.MarkowitzRow[row]++
		m.MarkowitzCol[col]++

		// Calculate Markowitz Product for Row
		if (m.MarkowitzRow[row] > LargestShortInteger && m.MarkowitzCol[row] != 0) ||
			(m.MarkowitzCol[row] > LargestShortInteger && m.MarkowitzRow[row] != 0) {
			product := float64(m.MarkowitzCol[row]) * float64(m.MarkowitzRow[row])
			if product >= float64(LargestLongInteger) {
				m.MarkowitzProd[row] = LargestLongInteger
			} else {
				m.MarkowitzProd[row] = int64(product)
			}
		} else {
			m.MarkowitzProd[row] = m.MarkowitzRow[row] * m.MarkowitzCol[row]
		}

		// Calculate Markowitz Product for Column
		if (m.MarkowitzRow[col] > LargestShortInteger && m.MarkowitzCol[col] != 0) ||
			(m.MarkowitzCol[col] > LargestShortInteger && m.MarkowitzRow[col] != 0) {
			product := float64(m.MarkowitzCol[col]) * float64(m.MarkowitzRow[col])
			if product >= float64(LargestLongInteger) {
				m.MarkowitzProd[col] = LargestLongInteger
			} else {
				m.MarkowitzProd[col] = int64(product)
			}
		} else {
			m.MarkowitzProd[col] = m.MarkowitzCol[col] * m.MarkowitzRow[col]
		}

		if m.MarkowitzRow[row] == 1 && m.MarkowitzCol[row] != 0 {
			m.Singletons--
		}
		if m.MarkowitzRow[col] != 0 && m.MarkowitzCol[col] == 1 {
			m.Singletons--
		}
	} else {
		element = &Element{Row: row, Col: col, Real: 0.0, Imag: 0.0}
		m.NeedsOrdering = true
	}

	if m.Complex {
		element.Imag = 0.0
	}
	if m.Config.Initialize {
		// element.InitInfo = &ComplexNumber{Real: element.Real, Imag: element.Imag}
		element.InitInfo = nil
	}

	m.Elements++

	element.NextInCol = current
	*prev = element

	if m.RowsLinked {
		current = *firstInRow
		prev = firstInRow
		for current != nil && current.Col < col {
			prev = &current.NextInRow
			current = current.NextInRow
		}
		element.NextInRow = current
		*prev = element
	}

	if row == col {
		m.Diags[row] = element
	}

	return element
}

func (m *Matrix) GetElement(row, col int64) *Element {
	if row < 0 || col < 0 {
		return nil
	}
	if row == 0 || col == 0 {
		return &Element{}
		// return m.TrashCan
	}

	internalRow, internalCol := row, col
	switch {
	case m.Config.Translate:
		err := m.Translate(&internalRow, &internalCol)
		if err != nil {
			return nil
		}
	default:
		if row > m.Size || col > m.Size {
			if m.Config.Expandable {
				newSize := max(row, col)
				err := m.EnlargeMatrix(newSize)
				if err != nil {
					return nil
				}
			} else {
				if m.Reordered {
					panic("Set Translate to add elements to a reordered matrix")
				}
				return nil
			}
		}
	}

	if internalRow == internalCol {
		if element := m.Diags[internalRow]; element != nil {
			return element
		}
	}

	element := m.FirstInCol[internalCol]
	for element != nil {
		if element.Row == internalRow {
			return element
		}
		element = element.NextInCol
	}

	return m.createElement(internalRow, internalCol, &m.FirstInRow[internalRow], &m.FirstInCol[internalCol], false)
}

func (m *Matrix) GetAdmittance(node1, node2 int64, template *Template) error {
	template.Element1 = m.GetElement(node1, node1)
	template.Element2 = m.GetElement(node2, node2)
	template.Element3Negated = m.GetElement(node2, node1)
	template.Element4Negated = m.GetElement(node1, node2)

	if template.Element1 == nil || template.Element2 == nil || template.Element3Negated == nil || template.Element4Negated == nil {
		return fmt.Errorf("memory allocation failed")
	}

	if node1 == 0 {
		template.Element1, template.Element2 = template.Element2, template.Element1
	}

	return nil
}

func (m *Matrix) LinkRows() {
	for col := m.Size; col >= 1; col-- {
		m.FirstInRow[col] = nil
	}

	for col := m.Size; col >= 1; col-- {
		element := m.FirstInCol[col]
		for element != nil {
			element.Col = col
			element.NextInRow = m.FirstInRow[element.Row]
			m.FirstInRow[element.Row] = element
			element = element.NextInCol
		}
	}

	m.RowsLinked = true
}

func (m *Matrix) Translate(row, col *int64) error {
	var err error
	var intRow, intCol, extRow, extCol int64

	extRow, extCol = *row, *col

	if extRow > m.ExtSize || extCol > m.ExtSize {
		err = m.ExpandTranslationArrays(max(extRow, extCol))
		if err != nil {
			return fmt.Errorf("failed to expand translation arrays: %v", err)
		}
		m.ExtSize = max(extRow, extCol)
	}

	// Translate row number: external to internal
	intRow = m.ExtToIntRowMap[extRow]
	if intRow == -1 {
		m.CurrentSize++
		m.ExtToIntRowMap[extRow] = m.CurrentSize
		m.ExtToIntColMap[extRow] = m.CurrentSize
		intRow = m.CurrentSize

		if intRow > m.Size {
			if !m.Config.Expandable {
				return fmt.Errorf("matrix size fixed")
			}
			err = m.EnlargeMatrix(intRow)
			if err != nil {
				return fmt.Errorf("failed to enlarge matrix: %v", err)
			}
		}

		m.IntToExtRowMap[intRow] = extRow
		m.IntToExtColMap[intRow] = extRow
	}

	// Translate column number: external to internal
	intCol = m.ExtToIntColMap[extCol]
	if intCol == -1 {
		m.CurrentSize++
		m.ExtToIntRowMap[extCol] = m.CurrentSize
		m.ExtToIntColMap[extCol] = m.CurrentSize
		intCol = m.CurrentSize

		if intCol > m.Size {
			if !m.Config.Expandable {
				return fmt.Errorf("matrix size fixed")
			}
			err = m.EnlargeMatrix(intCol)
			if err != nil {
				return fmt.Errorf("failed to enlarge matrix: %v", err)
			}
		}

		m.IntToExtRowMap[intCol] = extCol
		m.IntToExtColMap[intCol] = extCol
	}

	*row = intRow
	*col = intCol

	return nil
}

func (m *Matrix) EnlargeMatrix(newSize int64) error {
	if newSize <= m.Size {
		return nil
	}

	m.Size = newSize
	newMatrixSize := m.Size + 1 // 1-based indexing

	newDiags := make([]*Element, newMatrixSize)
	newFirstInRow := make([]*Element, newMatrixSize)
	newFirstInCol := make([]*Element, newMatrixSize)
	newIntToExtColMap := make([]int64, newMatrixSize)
	newIntToExtRowMap := make([]int64, newMatrixSize)

	copy(newDiags, m.Diags)
	copy(newFirstInRow, m.FirstInRow)
	copy(newFirstInCol, m.FirstInCol)
	copy(newIntToExtColMap, m.IntToExtColMap)
	copy(newIntToExtRowMap, m.IntToExtRowMap)

	for i := int64(len(m.IntToExtColMap)); i < newMatrixSize; i++ {
		newIntToExtColMap[i] = i
		newIntToExtRowMap[i] = i
	}

	m.Diags = newDiags
	m.FirstInRow = newFirstInRow
	m.FirstInCol = newFirstInCol
	m.IntToExtColMap = newIntToExtColMap
	m.IntToExtRowMap = newIntToExtRowMap

	m.MarkowitzRow = nil
	m.MarkowitzCol = nil
	m.MarkowitzProd = nil
	m.DoRealDirect = nil
	m.DoComplexDirect = nil
	m.Intermediate = nil
	m.InternalVectorsAllocated = false

	// Use at OrderAndFactor insetead of this
	// if err := m.CreateInternalVectors(); err != nil {
	// 	return fmt.Errorf("failed to create internal vectors: %v", err)
	// }

	return nil
}

func (m *Matrix) ExpandTranslationArrays(newSize int64) error {
	if newSize <= m.ExtSize {
		return nil
	}

	m.ExtSize = newSize
	newMatrixSize := newSize + 1 // 1-based indexing

	newExtToIntRowMap := make([]int64, newMatrixSize)
	newExtToIntColMap := make([]int64, newMatrixSize)

	copy(newExtToIntRowMap, m.ExtToIntRowMap)
	copy(newExtToIntColMap, m.ExtToIntColMap)

	for i := int64(len(m.ExtToIntRowMap)); i < (newMatrixSize); i++ {
		newExtToIntRowMap[i] = -1
		newExtToIntColMap[i] = -1
	}

	m.ExtToIntRowMap = newExtToIntRowMap
	m.ExtToIntColMap = newExtToIntColMap

	return nil
}

func (m *Matrix) CreateInternalVectorsNotUse() error {
	matrixSize := m.Size + 1 // 1-based indexing

	if m.Complex {
		m.Intermediate = make([]float64, 2*(matrixSize))
	} else {
		m.Intermediate = make([]float64, matrixSize)
	}

	m.InternalVectorsAllocated = true
	return nil
}

func (m *Matrix) CreateInternalVectors() error {
	matrixSize := m.Size + 1 // 1-based indexing

	m.MarkowitzRow = make([]int64, matrixSize)
	m.MarkowitzCol = make([]int64, matrixSize)
	m.MarkowitzProd = make([]int64, matrixSize+1)

	if m.Complex {
		m.DoComplexDirect = make([]bool, matrixSize)
	} else {
		m.DoRealDirect = make([]bool, matrixSize)
	}

	if m.Complex {
		m.Intermediate = make([]float64, 2*(matrixSize))
	} else {
		m.Intermediate = make([]float64, matrixSize)
	}

	m.InternalVectorsAllocated = true
	return nil
}
