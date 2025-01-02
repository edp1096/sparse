package sparse

import (
	"math"
)

func (m *Matrix) SearchForPivot(step int64, diagPivoting bool) *Element {
	if m.Singletons > 0 {
		if pivot := m.SearchForSingleton(step); pivot != nil {
			m.PivotSelectionMethod = 's'
			return pivot
		}
	}

	if diagPivoting {
		if pivot := m.QuicklySearchDiagonal(step); pivot != nil {
			m.PivotSelectionMethod = 'q'
			return pivot
		}

		if pivot := m.SearchDiagonal(step); pivot != nil {
			m.PivotSelectionMethod = 'd'
			return pivot
		}
	}

	pivot := m.SearchEntireMatrix(step)
	m.PivotSelectionMethod = 'e'

	return pivot
}

func (m *Matrix) FindBiggestInColExclude(elem *Element, step int64) float64 {
	current := m.FirstInCol[elem.Col]

	for current != nil && current.Row < step {
		current = current.NextInCol
	}

	if current == nil {
		return 0.0
	}

	var largest float64
	if current.Row != elem.Row {
		largest = m.elementMag(current)
	} else {
		largest = 0.0
	}

	for current = current.NextInCol; current != nil; current = current.NextInCol {
		magnitude := m.elementMag(current)
		if magnitude > largest && current.Row != elem.Row {
			largest = magnitude
		}
	}

	return largest
}

func (m *Matrix) FindBiggestInCol(element *Element) float64 {
	largest := 0.0

	current := element
	for current != nil {
		magnitude := m.elementMag(current)
		if magnitude > largest {
			largest = magnitude
		}
		current = current.NextInCol
	}

	return largest
}

func (m *Matrix) SearchForSingleton(step int64) *Element {
	m.MarkowitzProd[m.Size+1] = m.MarkowitzProd[step]
	m.MarkowitzProd[step-1] = 0

	singletons := m.Singletons
	m.Singletons--

	pMarkowitzProduct := m.Size + 1

	for singletons > 0 {
		for pMarkowitzProduct >= step && m.MarkowitzProd[pMarkowitzProduct] != 0 {
			pMarkowitzProduct--
		}

		i := pMarkowitzProduct
		if i < step {
			break
		}
		if i > m.Size {
			i = step
		}

		var chosenPivot *Element
		if pivot := m.Diags[i]; pivot != nil {
			pivotMag := m.elementMag(pivot)
			if pivotMag > m.AbsThreshold &&
				pivotMag > m.RelThreshold*m.FindBiggestInColExclude(pivot, step) {
				return pivot
			}
		} else {
			if m.MarkowitzCol[i] == 0 {
				pivot := m.FirstInCol[i]
				for pivot != nil && pivot.Row < step {
					pivot = pivot.NextInCol
				}
				if pivot != nil {
					chosenPivot = pivot
				}
			}
			if chosenPivot == nil && m.MarkowitzRow[i] == 0 {
				pivot := m.FirstInRow[i]
				for pivot != nil && pivot.Col < step {
					pivot = pivot.NextInRow
				}
				if pivot != nil {
					chosenPivot = pivot
				}
			}

			if chosenPivot != nil {
				pivotMag := m.elementMag(chosenPivot)
				if pivotMag > m.AbsThreshold && pivotMag > m.RelThreshold*m.FindBiggestInColExclude(chosenPivot, step) {
					return chosenPivot
				}
			}
		}

		singletons--
		pMarkowitzProduct--
	}

	m.Singletons++
	return nil
}

func (m *Matrix) QuicklySearchDiagonal(step int64) *Element {
	var chosenPivot *Element

	minMarkowitzProduct := int64(math.MaxInt64)
	m.MarkowitzProd[m.Size+1] = m.MarkowitzProd[step]
	m.MarkowitzProd[step-1] = -1

	index := m.Size + 2
	for {
		index--
		for m.MarkowitzProd[index] >= minMarkowitzProduct {
			index--
		}

		i := index

		if i < step {
			break
		}
		if i > m.Size {
			i = step
		}

		diag := m.Diags[i]
		if diag == nil {
			continue
		}
		magnitude := m.elementMag(diag)
		if magnitude <= m.AbsThreshold {
			continue
		}

		if m.MarkowitzProd[i] == 1 {
			otherInRow := diag.NextInRow
			otherInCol := diag.NextInCol

			if otherInRow == nil && otherInCol == nil {
				otherInRow := m.FirstInRow[i]
				for otherInRow != nil {
					if otherInRow.Col >= step && otherInRow.Col != i {
						break
					}
					otherInRow = otherInRow.NextInRow
				}

				otherInCol := m.FirstInCol[i]
				for otherInCol != nil {
					if otherInCol.Row >= step && otherInCol.Row != i {
						break
					}
					otherInCol = otherInCol.NextInCol
				}
			}

			if otherInRow != nil && otherInCol != nil {
				if otherInRow.Col == otherInCol.Row {
					largestOffDiag := math.Max(m.elementMag(otherInRow), m.elementMag(otherInCol))
					if magnitude >= largestOffDiag {
						return diag
					}
				}
			}
		}

		minMarkowitzProduct = m.MarkowitzProd[i]
		chosenPivot = diag
	}

	if chosenPivot != nil {
		largestInCol := m.FindBiggestInColExclude(chosenPivot, step)

		if m.elementMag(chosenPivot) <= m.RelThreshold*largestInCol {
			chosenPivot = nil
		}
	}

	return chosenPivot
}

func (m *Matrix) SearchDiagonal(step int64) *Element {
	var chosenPivot *Element
	minMarkowitzProduct := int64(math.MaxInt64)
	numberOfTies := int64(0)
	var ratioOfAccepted float64

	m.MarkowitzProd[m.Size+1] = m.MarkowitzProd[step]

	for i := m.Size; i >= step; i-- {
		if m.MarkowitzProd[i] > minMarkowitzProduct {
			continue
		}

		diag := m.Diags[i]
		if diag == nil {
			continue
		}

		magnitude := m.elementMag(diag)
		if magnitude <= m.AbsThreshold {
			continue
		}

		largestInCol := m.FindBiggestInColExclude(diag, step)
		if magnitude <= m.RelThreshold*largestInCol {
			continue
		}

		if m.MarkowitzProd[i] < minMarkowitzProduct {
			chosenPivot = diag
			minMarkowitzProduct = m.MarkowitzProd[i]
			ratioOfAccepted = largestInCol / magnitude
			numberOfTies = 0
		} else {
			numberOfTies++
			ratio := largestInCol / magnitude
			if ratio < ratioOfAccepted {
				chosenPivot = diag
				ratioOfAccepted = ratio
			}
			if numberOfTies >= minMarkowitzProduct*int64(m.Config.TiesMultiplier) {
				return chosenPivot
			}
		}
	}

	return chosenPivot
}

func (m *Matrix) SearchEntireMatrix(step int64) *Element {
	var chosenPivot *Element
	var pLargestElement *Element
	minMarkowitzProduct := int64(math.MaxInt64)
	largestElementMag := 0.0
	numberOfTies := int64(0)
	var ratioOfAccepted float64

	for i := step; i <= m.Size; i++ {
		current := m.FirstInCol[i]

		for current != nil && current.Row < step {
			current = current.NextInCol
		}

		largestInCol := m.FindBiggestInCol(current)
		if largestInCol == 0.0 {
			continue
		}

		for current != nil {
			magnitude := m.elementMag(current)
			if magnitude > largestElementMag {
				largestElementMag = magnitude
				pLargestElement = current
			}

			product := m.calculateMarkowitzProduct(m.MarkowitzRow[current.Row], m.MarkowitzCol[current.Col])

			if product <= minMarkowitzProduct && magnitude > m.RelThreshold*largestInCol && magnitude > m.AbsThreshold {
				if product < minMarkowitzProduct {
					chosenPivot = current
					minMarkowitzProduct = product
					ratioOfAccepted = largestInCol / magnitude
					numberOfTies = 0
				} else {
					numberOfTies++
					ratio := largestInCol / magnitude
					if ratio < ratioOfAccepted {
						chosenPivot = current
						ratioOfAccepted = ratio
					}
					if numberOfTies >= minMarkowitzProduct*int64(m.Config.TiesMultiplier) {
						return chosenPivot
					}
				}
			}
			current = current.NextInCol
		}
	}

	if chosenPivot != nil {
		return chosenPivot
	}

	if largestElementMag == 0.0 {
		// TODO: append error
		return nil
	}

	// TODO: append small pivot error
	return pLargestElement
}
