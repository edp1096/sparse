package sparse

// modified node admittance matrix - remove the diagonal zeros

func (m *Matrix) MNAPreorder() {
	if m.RowsLinked {
		return
	}

	m.Reordered = true
	size := m.Size
	startAt := int64(1)

	for {
		anotherPassNeeded := false
		swapped := false

		for j := startAt; j <= size; j++ {
			if m.Diags[j] == nil {
				twins, pTwin1, pTwin2 := m.CountTwins(j)

				if twins == 1 {
					m.SwapCols(pTwin1, pTwin2)
					swapped = true
				} else if twins > 1 && !anotherPassNeeded {
					anotherPassNeeded = true
					startAt = j
				}
			}
		}

		if anotherPassNeeded {
			for j := startAt; !swapped && j <= size; j++ {
				if m.Diags[j] == nil {
					_, pTwin1, pTwin2 := m.CountTwins(j)
					m.SwapCols(pTwin1, pTwin2)
					swapped = true
				}
			}
		}

		if !anotherPassNeeded {
			break
		}
	}
}

func (m *Matrix) CountTwins(col int64) (twins int, retTwin1, retTwin2 *Element) {
	pTwin1 := m.FirstInCol[col]
	for pTwin1 != nil {
		if m.elementMag(pTwin1) == 1.0 {
			row := pTwin1.Row
			pTwin2 := m.FirstInCol[row]
			for pTwin2 != nil && pTwin2.Row != col {
				pTwin2 = pTwin2.NextInCol
			}
			if pTwin2 != nil && m.elementMag(pTwin2) == 1.0 {
				twins++
				if twins >= 2 {
					return
				}
				retTwin1 = pTwin1
				retTwin2 = pTwin2
				retTwin1.Col = col
				retTwin2.Col = row
			}
		}
		pTwin1 = pTwin1.NextInCol
	}

	return
}

func (m *Matrix) SwapCols(pTwin1, pTwin2 *Element) {
	col1 := pTwin1.Col
	col2 := pTwin2.Col

	m.FirstInCol[col1], m.FirstInCol[col2] = m.FirstInCol[col2], m.FirstInCol[col1]
	m.IntToExtColMap[col1], m.IntToExtColMap[col2] = m.IntToExtColMap[col2], m.IntToExtColMap[col1]

	if m.Config.Translate {
		m.ExtToIntColMap[m.IntToExtColMap[col2]] = col2
		m.ExtToIntColMap[m.IntToExtColMap[col1]] = col1
	}

	m.Diags[col1] = pTwin2
	m.Diags[col2] = pTwin1

	m.NumberOfInterchangesIsOdd = !m.NumberOfInterchangesIsOdd
}
