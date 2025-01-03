package sparse // import "sparse"

import (
	"fmt"
)

func Create(size int64, config *Configuration) (*Matrix, error) {
	if size <= 0 {
		return nil, fmt.Errorf("invalid size: %d", size)
	}

	defaultConfig := Configuration{
		Real:                    true,
		Complex:                 true,
		SeparatedComplexVectors: false,
		Expandable:              true,
		Translate:               true,
		ModifiedNodal:           true,
		TiesMultiplier:          5,
		PrinterWidth:            80,
		Annotate:                0,
	}

	if config == nil {
		config = &defaultConfig
	}

	matrixSize := size + 1 // 1-based indexing

	m := &Matrix{
		Config:          *config,
		Size:            size,
		Complex:         config.Complex,
		DoRealDirect:    make([]bool, matrixSize+1),
		DoComplexDirect: make([]bool, matrixSize+1),
		CurrentSize:     0,
		Elements:        0,
		Diags:           make([]*Element, matrixSize+1),
		FirstInRow:      make([]*Element, matrixSize+1),
		FirstInCol:      make([]*Element, matrixSize+1),
		Intermediate:    make([]float64, matrixSize+1),
		MarkowitzRow:    make([]int64, matrixSize+1),
		MarkowitzCol:    make([]int64, matrixSize+1),
		MarkowitzProd:   make([]int64, matrixSize+1),
		IntToExtRowMap:  make([]int64, matrixSize+1),
		IntToExtColMap:  make([]int64, matrixSize+1),
		ExtToIntRowMap:  make([]int64, matrixSize+1),
		ExtToIntColMap:  make([]int64, matrixSize+1),
		NeedsOrdering:   true,
		RelThreshold:    0.001,
		AbsThreshold:    0.0,
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
	// element.InitInfo = value
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
				// if m.Complex {
				// 	element.Imag = element.InitInfo.Imag
				// }
			}
			element = element.NextInCol
		}
	}

	m.Factored = false
	m.SingularCol = 0
	m.SingularRow = 0

	return nil
}

// func (m *Matrix) ClearNotUse() {
// 	m.Elements = 0
// 	m.Fillins = 0

// 	matrixSize := m.Size + 1

// 	m.Diags = make([]*Element, matrixSize+1)
// 	m.FirstInRow = make([]*Element, matrixSize+1)
// 	m.FirstInCol = make([]*Element, matrixSize+1)
// 	m.Intermediate = make([]float64, matrixSize+1)
// 	if m.Complex {
// 		m.Intermediate = make([]float64, matrixSize*2)
// 	}
// 	m.MarkowitzRow = make([]int64, matrixSize+1)
// 	m.MarkowitzCol = make([]int64, matrixSize+1)
// 	m.MarkowitzProd = make([]int64, matrixSize+1)
// }

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
	// if row < 1 || col < 1 || row > m.Size || col > m.Size {
	// 	return nil
	// }

	if row < 0 || col < 0 {
		return nil
	}
	if row == 0 || col == 0 {
		return &Element{}
	}

	internalRow, internalCol := row, col
	if m.Config.Translate {
		if err := m.Translate(&internalRow, &internalCol); err != nil {
			return nil
		}
	} else {
		if row > m.Size || col > m.Size {
			return nil
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
	var intRow, intCol, extRow, extCol int64

	extRow, extCol = *row, *col

	// Translate external row number to internal row number
	intRow = m.ExtToIntRowMap[extRow]
	if intRow == -1 {
		m.CurrentSize++
		m.ExtToIntRowMap[extRow] = m.CurrentSize
		m.ExtToIntColMap[extRow] = m.CurrentSize
		intRow = m.CurrentSize

		if !m.Config.Expandable && intRow > m.Size {
			return fmt.Errorf("matrix size fixed")
		}

		m.IntToExtRowMap[intRow] = extRow
		m.IntToExtColMap[intRow] = extRow
	}

	// Translate external column number to internal column number
	intCol = m.ExtToIntColMap[extCol]
	if intCol == -1 {
		m.CurrentSize++
		m.ExtToIntRowMap[extCol] = m.CurrentSize
		m.ExtToIntColMap[extCol] = m.CurrentSize
		intCol = m.CurrentSize

		if !m.Config.Expandable && intCol > m.Size {
			return fmt.Errorf("matrix size fixed")
		}

		m.IntToExtRowMap[intCol] = extCol
		m.IntToExtColMap[intCol] = extCol
	}

	*row = intRow
	*col = intCol

	return nil
}

func (m *Matrix) CreateInternalVectors() error {
	size := m.Size
	if m.Complex {
		m.Intermediate = make([]float64, 2*(size+1))
	} else {
		m.Intermediate = make([]float64, size+1)
	}

	m.InternalVectorsAllocated = true
	return nil
}
