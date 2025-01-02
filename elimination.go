package sparse

import (
	"fmt"
)

func (m *Matrix) ComplexRowColElimination(pivot *Element) error {
	if m.elementMag(pivot) == 0.0 {
		m.SingularRow = pivot.Row
		return fmt.Errorf("matrix is singular at row %d", pivot.Row)
	}

	m.complexReciprocal(pivot)
	pUpper := pivot.NextInRow

	for pUpper != nil {
		m.complexMultAssign(pUpper, pivot)
		pSub := pUpper.NextInCol
		pLower := pivot.NextInCol
		ppAbove := &pUpper.NextInCol

		for pLower != nil {
			row := pLower.Row
			for pSub != nil && pSub.Row < row {
				ppAbove = &pSub.NextInCol
				pSub = pSub.NextInCol
			}

			if pSub == nil || pSub.Row > row {
				pSub = m.createElement(row, pUpper.Col, &pLower.NextInRow, ppAbove, true)
				if pSub == nil {
					return fmt.Errorf("memory allocation failed")
				}
			}

			m.complexMultSubtAssign(pSub, pUpper, pLower)
			pSub = pSub.NextInCol
			pLower = pLower.NextInCol
		}
		pUpper = pUpper.NextInRow
	}
	return nil
}

func (m *Matrix) RealRowColElimination(pivot *Element) error {
	if m.elementMag(pivot) == 0.0 {
		m.SingularRow = pivot.Row
		return fmt.Errorf("matrix is singular at row %d", pivot.Row)
	}

	pivot.Real = 1.0 / pivot.Real
	pUpper := pivot.NextInRow

	for pUpper != nil {
		pUpper.Real *= pivot.Real

		pSub := pUpper.NextInCol
		pLower := pivot.NextInCol
		ppAbove := &pUpper.NextInCol
		for pLower != nil {
			row := pLower.Row

			for pSub != nil && pSub.Row < row {
				ppAbove = &pSub.NextInCol
				pSub = pSub.NextInCol
			}

			if pSub == nil || pSub.Row > row {
				pSub = m.createElement(row, pUpper.Col, &pLower.NextInRow, ppAbove, true)
				if pSub == nil {
					return fmt.Errorf("memory allocation failed")
				}
			}

			pSub.Real -= pUpper.Real * pLower.Real
			pSub = pSub.NextInCol
			pLower = pLower.NextInCol
		}
		pUpper = pUpper.NextInRow
	}

	return nil
}
