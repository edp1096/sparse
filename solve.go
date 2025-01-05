package sparse

import (
	"fmt"
)

func (m *Matrix) Solve(rhs []float64) (solution []float64, err error) {
	solution = make([]float64, len(rhs))

	if !m.Factored {
		return nil, fmt.Errorf("matrix is not factored")
	}
	if len(rhs) < int(m.Size) || len(solution) < int(m.Size) {
		return nil, fmt.Errorf("rhs or solution array size(%d,%d) is smaller than matrix size(%d)",
			len(rhs), len(solution), m.Size)
	}
	if m.Complex {
		return nil, fmt.Errorf("complex matrix not implemented yet")
	}
	if m.Intermediate == nil {
		return nil, fmt.Errorf("intermediate vector not allocated")
	}

	size := m.Size
	intermediate := m.Intermediate
	intToExtRowMap := m.IntToExtRowMap
	intToExtColMap := m.IntToExtColMap
	diags := m.Diags

	for i := size; i > 0; i-- {
		intermediate[i] = rhs[intToExtRowMap[i]]
	}

	// Forward elimination - Solves Lc = b
	for i := int64(1); i <= size; i++ {
		temp := intermediate[i]
		if temp != 0.0 {
			pivot := diags[i]
			if pivot == nil {
				return nil, fmt.Errorf("nil diagonal element at %d", i)
			}
			temp *= pivot.Real
			intermediate[i] = temp

			for element := pivot.NextInCol; element != nil; element = element.NextInCol {
				intermediate[element.Row] -= temp * element.Real
			}
		}
	}

	// Backward Substitution - Solves Ux = c
	for i := size; i > 0; i-- {
		temp := intermediate[i]

		for element := diags[i].NextInRow; element != nil; element = element.NextInRow {
			temp -= element.Real * intermediate[element.Col]
		}
		intermediate[i] = temp
	}

	// Unscramble Intermediate vector - reorder from internal to external ordering
	for i := size; i > 0; i-- {
		solution[intToExtColMap[i]] = intermediate[i]
	}

	return solution, nil
}

func (m *Matrix) SolveTransposed(rhs []float64) (solution []float64, err error) {
	solution = make([]float64, len(rhs))

	if !m.Factored {
		return nil, fmt.Errorf("matrix is not factored")
	}
	if len(rhs) < int(m.Size) || len(solution) < int(m.Size) {
		return nil, fmt.Errorf("rhs or solution array size(%d,%d) is smaller than matrix size(%d)",
			len(rhs), len(solution), m.Size)
	}
	if m.Complex {
		return nil, fmt.Errorf("complex matrix not implemented yet")
	}
	if m.Intermediate == nil {
		return nil, fmt.Errorf("intermediate vector not allocated")
	}

	size := m.Size
	intermediate := m.Intermediate
	intToExtRowMap := m.IntToExtRowMap
	intToExtColMap := m.IntToExtColMap
	diags := m.Diags

	// Initialize Intermediate vector - Convert from external to internal ordering
	for i := size; i > 0; i-- {
		intermediate[i] = rhs[intToExtColMap[i]]
	}

	// Forward elimination
	for i := int64(1); i <= size; i++ {
		temp := intermediate[i]
		if temp != 0.0 {
			pivot := diags[i]
			if pivot == nil {
				return nil, fmt.Errorf("nil diagonal element at %d", i)
			}

			for element := pivot.NextInRow; element != nil; element = element.NextInRow {
				intermediate[element.Col] -= temp * element.Real
			}
		}
	}

	// Backward Substitution
	for i := size; i > 0; i-- {
		pivot := diags[i]
		temp := intermediate[i]

		for element := pivot.NextInCol; element != nil; element = element.NextInCol {
			temp -= element.Real * intermediate[element.Row]
		}

		intermediate[i] = temp * pivot.Real
	}

	for i := size; i > 0; i-- {
		solution[intToExtRowMap[i]] = intermediate[i]
	}

	return solution, nil
}

func (m *Matrix) SolveComplex(rhs, irhs []float64) ([]float64, []float64, error) {
	size := m.Size
	matrixSize := size + 1 // 1-based indexing

	solution := make([]float64, 2*(matrixSize))

	if m.Config.SeparatedComplexVectors {
		for i := int64(1); i <= size; i++ {
			extIdx := m.IntToExtRowMap[i]
			m.Intermediate[i*2] = rhs[extIdx]
			m.Intermediate[i*2+1] = irhs[extIdx]
		}
	} else {
		for i := int64(1); i <= size; i++ {
			extIdx := m.IntToExtRowMap[i]
			temp := &Element{
				Real: rhs[extIdx*2],
				Imag: rhs[extIdx*2+1],
			}
			m.Intermediate[i*2] = temp.Real
			m.Intermediate[i*2+1] = temp.Imag
		}
	}

	// Forward substitution
	for i := int64(1); i <= size; i++ {
		temp := &Element{
			Real: m.Intermediate[i*2],
			Imag: m.Intermediate[i*2+1],
		}

		if temp.Real != 0.0 || temp.Imag != 0.0 {
			pivot := m.Diags[i]
			m.complexMultAssign(temp, pivot)

			m.Intermediate[i*2] = temp.Real
			m.Intermediate[i*2+1] = temp.Imag

			for element := pivot.NextInCol; element != nil; element = element.NextInCol {
				interm := &Element{
					Real: m.Intermediate[element.Row*2],
					Imag: m.Intermediate[element.Row*2+1],
				}
				m.complexMultSubtAssign(interm, temp, element)
				m.Intermediate[element.Row*2] = interm.Real
				m.Intermediate[element.Row*2+1] = interm.Imag
			}
		}
	}

	// Backward substitution
	for i := size; i > 0; i-- {
		temp := &Element{
			Real: m.Intermediate[i*2],
			Imag: m.Intermediate[i*2+1],
		}

		for element := m.Diags[i].NextInRow; element != nil; element = element.NextInRow {
			interm := &Element{
				Real: m.Intermediate[element.Col*2],
				Imag: m.Intermediate[element.Col*2+1],
			}
			m.complexMultSubtAssign(temp, element, interm)
		}

		m.Intermediate[i*2] = temp.Real
		m.Intermediate[i*2+1] = temp.Imag
	}

	if m.Config.SeparatedComplexVectors {
		solReal := make([]float64, matrixSize)
		solImag := make([]float64, matrixSize)
		for i := size; i > 0; i-- {
			extIdx := m.IntToExtColMap[i]
			solReal[extIdx] = m.Intermediate[i*2]
			solImag[extIdx] = m.Intermediate[i*2+1]
		}
		return solReal, solImag, nil
	} else {
		for i := size; i > 0; i-- {
			extIdx := m.IntToExtColMap[i]
			temp := Element{
				Real: m.Intermediate[i*2],
				Imag: m.Intermediate[i*2+1],
			}
			solution[extIdx*2] = temp.Real
			solution[extIdx*2+1] = temp.Imag
		}
		return solution, nil, nil
	}
}

func (m *Matrix) SolveComplexTransposed(rhs, irhs []float64) ([]float64, []float64, error) {
	size := m.Size
	matrixSize := m.Size + 1 // 1-based indexing

	if !m.Factored {
		return nil, nil, fmt.Errorf("matrix is not factored")
	}
	if len(rhs) < int(size) || len(irhs) < int(size) {
		return nil, nil, fmt.Errorf("rhs or irhs array size(%d,%d) is smaller than matrix size(%d)",
			len(rhs), len(irhs), size)
	}
	if !m.Complex {
		return nil, nil, fmt.Errorf("matrix must be complex")
	}
	if m.Intermediate == nil {
		return nil, nil, fmt.Errorf("intermediate vector not allocated")
	}

	// Initialize vectors
	if !m.Config.SeparatedComplexVectors {
		for i := int64(1); i <= size; i++ {
			extIdx := m.IntToExtColMap[i]
			temp := &Element{
				Real: rhs[extIdx*2],
				Imag: rhs[extIdx*2+1],
			}
			m.Intermediate[i*2] = temp.Real
			m.Intermediate[i*2+1] = temp.Imag
		}
	} else {
		for i := int64(1); i <= size; i++ {
			extIdx := m.IntToExtColMap[i]
			m.Intermediate[i*2] = rhs[extIdx]
			m.Intermediate[i*2+1] = irhs[extIdx]
		}
	}

	// Forward elimination
	for i := int64(1); i <= size; i++ {
		temp := &Element{
			Real: m.Intermediate[i*2],
			Imag: m.Intermediate[i*2+1],
		}

		if temp.Real != 0.0 || temp.Imag != 0.0 {
			for element := m.Diags[i].NextInRow; element != nil; element = element.NextInRow {
				interm := &Element{
					Real: m.Intermediate[element.Col*2],
					Imag: m.Intermediate[element.Col*2+1],
				}
				m.complexMultSubtAssign(interm, temp, element)
				m.Intermediate[element.Col*2] = interm.Real
				m.Intermediate[element.Col*2+1] = interm.Imag
			}
		}
	}

	// Backward substitution
	for i := size; i > 0; i-- {
		pivot := m.Diags[i]
		temp := &Element{
			Real: m.Intermediate[i*2],
			Imag: m.Intermediate[i*2+1],
		}

		for element := pivot.NextInCol; element != nil; element = element.NextInCol {
			interm := &Element{
				Real: m.Intermediate[element.Row*2],
				Imag: m.Intermediate[element.Row*2+1],
			}
			m.complexMultSubtAssign(temp, interm, element)
		}

		m.complexMultAssign(temp, pivot)
		m.Intermediate[i*2] = temp.Real
		m.Intermediate[i*2+1] = temp.Imag
	}

	if m.Config.SeparatedComplexVectors {
		solReal := make([]float64, matrixSize)
		solImag := make([]float64, matrixSize)
		for i := size; i > 0; i-- {
			extIdx := m.IntToExtRowMap[i]
			solReal[extIdx] = m.Intermediate[i*2]
			solImag[extIdx] = m.Intermediate[i*2+1]
		}
		return solReal, solImag, nil
	} else {
		solution := make([]float64, 2*(matrixSize))
		for i := size; i > 0; i-- {
			extIdx := m.IntToExtRowMap[i]
			solution[extIdx*2] = m.Intermediate[i*2]
			solution[extIdx*2+1] = m.Intermediate[i*2+1]
		}
		return solution, nil, nil
	}
}
