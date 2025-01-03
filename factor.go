package sparse

import (
	"fmt"
)

func (m *Matrix) OrderAndFactor(rhs []float64, relThreshold, absThreshold float64, diagPivoting bool) error {
	if relThreshold <= 0.0 || relThreshold > 1.0 {
		relThreshold = m.RelThreshold
	}
	m.RelThreshold = relThreshold
	if absThreshold < 0.0 {
		absThreshold = m.AbsThreshold
	}
	m.AbsThreshold = absThreshold

	size := m.Size
	var step int64 = 1

	if !m.NeedsOrdering {
		for step = 1; step <= size; step++ {
			pivot := m.Diags[step]
			if pivot == nil {
				m.NeedsOrdering = true
				break
			}

			largestInCol := m.FindBiggestInCol(pivot.NextInCol)
			if largestInCol*relThreshold < m.elementMag(pivot) {
				if m.Complex {
					m.ComplexRowColElimination(pivot)
				} else {
					m.RealRowColElimination(pivot)
				}
			} else {
				m.NeedsOrdering = true
				break
			}
		}

		if !m.NeedsOrdering {
			m.Factored = true
			return nil
		}
	} else {
		step = 1
		if !m.RowsLinked {
			m.LinkRows()
		}
	}

	m.CountMarkowitz(rhs, step)
	m.MarkowitzProducts(step)
	m.MaxRowCountInLowerTri = -1

	for ; step <= size; step++ {
		pivot := m.SearchForPivot(step, diagPivoting)
		if pivot == nil {
			m.SingularRow = step
			m.SingularCol = step
			return fmt.Errorf("matrix is singular at step %d", step)
		}

		m.ExchangeRowsAndCols(pivot, step)

		if m.Complex {
			m.ComplexRowColElimination(pivot)
		} else {
			m.RealRowColElimination(pivot)
		}

		m.UpdateMarkowitzNumbers(pivot)

		if m.Config.Annotate > 0 {
			m.WriteStatus(step)
		}
	}

	m.NeedsOrdering = false
	m.Reordered = true
	m.Factored = true
	return nil
}

func (m *Matrix) Factor() error {
	if m.NeedsOrdering {
		return m.OrderAndFactor(nil, 0.0, 0.0, true)
	}

	if !m.Partitioned {
		if err := m.Partition(DEFAULT_PARTITION); err != nil {
			return err
		}
	}

	if m.Complex {
		return m.FactorComplex()
	}

	if m.Diags[1] == nil || m.Diags[1].Real == 0.0 {
		m.SingularRow = 1
		m.SingularCol = 1
		return fmt.Errorf("zero pivot at step 1")
	}

	m.Diags[1].Real = 1.0 / m.Diags[1].Real

	for step := int64(2); step <= m.Size; step++ {
		if m.DoRealDirect[step] {
			// factorization - Direct
			for element := m.FirstInCol[step]; element != nil; element = element.NextInCol {
				m.Intermediate[element.Row] = element.Real
			}

			pColumn := m.FirstInCol[step]
			for pColumn != nil && pColumn.Row < step {
				element := m.Diags[pColumn.Row]
				pColumn.Real = m.Intermediate[pColumn.Row] * element.Real
				for element = element.NextInCol; element != nil; element = element.NextInCol {
					m.Intermediate[element.Row] -= pColumn.Real * element.Real
				}
				pColumn = pColumn.NextInCol
			}

			for element := m.Diags[step].NextInCol; element != nil; element = element.NextInCol {
				element.Real = m.Intermediate[element.Row]
			}

			if m.Intermediate[step] == 0.0 {
				return fmt.Errorf("zero pivot at step %d", step)
			}
			m.Diags[step].Real = 1.0 / m.Intermediate[step]
		} else {
			// factorization - Indirect
			dest := make([]*float64, m.Size+1)

			for element := m.FirstInCol[step]; element != nil; element = element.NextInCol {
				dest[element.Row] = &element.Real
			}

			for column := m.FirstInCol[step]; column != nil && column.Row < step; column = column.NextInCol {
				diag := m.Diags[column.Row]
				if diag == nil {
					continue
				}

				pColReal := *dest[column.Row] * diag.Real
				*dest[column.Row] = pColReal
				for element := diag.NextInCol; element != nil; element = element.NextInCol {
					*dest[element.Row] -= pColReal * element.Real
				}
			}

			diag := m.Diags[step]
			if diag == nil || diag.Real == 0.0 {
				m.SingularRow = step
				m.SingularCol = step
				return fmt.Errorf("zero pivot at step %d", step)
			}
			diag.Real = 1.0 / diag.Real
		}
	}

	m.Factored = true
	return nil
}

func (m *Matrix) FactorComplex() error {
	if m.Diags[1] == nil || (m.Diags[1].Real*m.Diags[1].Real+m.Diags[1].Imag*m.Diags[1].Imag) == 0 {
		m.SingularRow = 1
		m.SingularCol = 1
		return fmt.Errorf("zero pivot at step 1")
	}

	m.complexReciprocal(m.Diags[1])

	for step := int64(2); step <= m.Size; step++ {
		if m.DoComplexDirect[step] {
			dest := make([]*Element, m.Size+1)
			for i := range dest {
				dest[i] = &Element{}
			}
			for element := m.FirstInCol[step]; element != nil; element = element.NextInCol {
				dest[element.Row].Real = element.Real
				dest[element.Row].Imag = element.Imag
			}

			for column := m.FirstInCol[step]; column.Row < step; column = column.NextInCol {
				element := m.Diags[column.Row]
				m.complexMultAssign(dest[column.Row], element)
				column.Real = dest[column.Row].Real
				column.Imag = dest[column.Row].Imag

				for element = element.NextInCol; element != nil; element = element.NextInCol {
					m.complexMultSubtAssign(dest[element.Row], dest[column.Row], element)
				}
			}

			for element := m.FirstInCol[step]; element != nil; element = element.NextInCol {
				element.Real = dest[element.Row].Real
				element.Imag = dest[element.Row].Imag
			}

			if dest[step].Real*dest[step].Real+dest[step].Imag*dest[step].Imag == 0.0 {
				m.SingularRow = step
				m.SingularCol = step
				return fmt.Errorf("zero pivot at step %d", step)
			}

			m.complexReciprocal(dest[step])
			m.Diags[step].Real = dest[step].Real
			m.Diags[step].Imag = dest[step].Imag
		} else {
			dest := make([]*Element, m.Size+1)
			for element := m.FirstInCol[step]; element != nil; element = element.NextInCol {
				dest[element.Row] = element
			}

			for column := m.FirstInCol[step]; column.Row < step; column = column.NextInCol {
				element := m.Diags[column.Row]
				m.complexMultAssign(dest[column.Row], element)

				for element = element.NextInCol; element != nil; element = element.NextInCol {
					m.complexMultSubtAssign(dest[element.Row], dest[column.Row], element)
				}
			}

			if m.Diags[step].Real*m.Diags[step].Real+m.Diags[step].Imag*m.Diags[step].Imag == 0.0 {
				m.SingularRow = step
				m.SingularCol = step
				return fmt.Errorf("zero pivot at step %d", step)
			}
			m.complexReciprocal(m.Diags[step])
		}
	}

	m.Factored = true
	return nil
}

func (m *Matrix) Partition(mode int) error {
	if m.Partitioned {
		return nil
	}

	size := m.Size
	m.Partitioned = true

	if mode == DEFAULT_PARTITION {
		mode = m.Config.DefaultPartition
	}

	switch mode {
	case DIRECT_PARTITION:
		for step := int64(1); step <= size; step++ {
			if m.Config.Real {
				m.DoRealDirect[step] = true
			}
			if m.Config.Complex {
				m.DoComplexDirect[step] = true
			}
		}
		return nil
	case INDIRECT_PARTITION:
		for step := int64(1); step <= size; step++ {
			if m.Config.Real {
				m.DoRealDirect[step] = false
			}
			if m.Config.Complex {
				m.DoComplexDirect[step] = false
			}
		}
		return nil
	default:
		if mode != AUTO_PARTITION {
			return fmt.Errorf("unknown partition mode: %d", mode)
		}
	}

	nc := m.MarkowitzRow
	no := m.MarkowitzCol
	nm := m.MarkowitzProd

	for step := int64(1); step <= size; step++ {
		nc[step] = 0
		no[step] = 0
		nm[step] = 0

		for pElement := m.FirstInCol[step]; pElement != nil; pElement = pElement.NextInCol {
			nc[step]++
		}

		pColumn := m.FirstInCol[step]
		for pColumn != nil && pColumn.Row < step {
			pElement := m.Diags[pColumn.Row]
			nm[step]++
			for pElement = pElement.NextInCol; pElement != nil; pElement = pElement.NextInCol {
				no[step]++
			}
			pColumn = pColumn.NextInCol
		}
	}

	for step := int64(1); step <= size; step++ {
		if m.Config.Real {
			m.DoRealDirect[step] = (nm[step]+no[step] > 3*nc[step]-2*nm[step])
		}
		if m.Config.Complex {
			m.DoComplexDirect[step] = (nm[step]+no[step] > 7*nc[step]-4*nm[step])
		}
	}

	if m.Config.Annotate == 2 {
		var ops int64
		for step := int64(1); step <= size; step++ {
			ops += no[step]
		}
		m.OperationCount = int(ops)
	}

	return nil
}
