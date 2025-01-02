package sparse

func (m *Matrix) findDiag(index int64) *Element {
	element := m.FirstInCol[index]

	for element != nil && element.Row < index {
		element = element.NextInCol
	}

	if element != nil && element.Row == index {
		return element
	}

	return nil
}

func (m *Matrix) ExchangeRowsAndCols(pivot *Element, step int64) {
	row := pivot.Row
	col := pivot.Col

	m.PivotsOriginalRow = row
	m.PivotsOriginalCol = col

	if row == step && col == step {
		return
	}

	if row == col {
		m.rowExchange(step, row)
		m.colExchange(step, col)

		m.MarkowitzProd[step], m.MarkowitzProd[row] = m.MarkowitzProd[row], m.MarkowitzProd[step]
		m.Diags[row], m.Diags[step] = m.Diags[step], m.Diags[row]
	} else {
		oldMarkowitzStep := m.MarkowitzProd[step]
		oldMarkowitzRow := m.MarkowitzProd[row]
		oldMarkowitzCol := m.MarkowitzProd[col]

		if row != step {
			m.rowExchange(step, row)
			m.MarkowitzProd[row] = m.calculateMarkowitzProduct(m.MarkowitzRow[row], m.MarkowitzCol[row])

			if (m.MarkowitzProd[row] == 0) != (oldMarkowitzRow == 0) {
				if oldMarkowitzRow == 0 {
					m.Singletons--
				} else {
					m.Singletons++
				}
			}
		}

		if col != step {
			m.colExchange(step, col)
			m.MarkowitzProd[col] = m.calculateMarkowitzProduct(m.MarkowitzCol[col], m.MarkowitzRow[col])

			if (m.MarkowitzProd[col] == 0) != (oldMarkowitzCol == 0) {
				if oldMarkowitzCol == 0 {
					m.Singletons--
				} else {
					m.Singletons++
				}
			}
			m.Diags[col] = m.findDiag(col)
		}

		if row != step {
			m.Diags[row] = m.findDiag(row)
		}
		m.Diags[step] = m.findDiag(step)

		m.MarkowitzProd[step] = int64(m.MarkowitzRow[step]) * int64(m.MarkowitzCol[step])
		if (m.MarkowitzProd[step] == 0) != (oldMarkowitzStep == 0) {
			if oldMarkowitzStep == 0 {
				m.Singletons--
			} else {
				m.Singletons++
			}
		}
	}
}

func (m *Matrix) rowExchange(row1, row2 int64) {
	if row1 > row2 {
		row1, row2 = row2, row1
	}

	row1Ptr := m.FirstInRow[row1]
	row2Ptr := m.FirstInRow[row2]

	for row1Ptr != nil || row2Ptr != nil {
		var column int64
		var element1, element2 *Element

		switch {
		case row1Ptr == nil:
			column = row2Ptr.Col
			element2 = row2Ptr
			row2Ptr = row2Ptr.NextInRow
		case row2Ptr == nil:
			column = row1Ptr.Col
			element1 = row1Ptr
			row1Ptr = row1Ptr.NextInRow
		case row1Ptr.Col < row2Ptr.Col:
			column = row1Ptr.Col
			element1 = row1Ptr
			row1Ptr = row1Ptr.NextInRow
		case row1Ptr.Col > row2Ptr.Col:
			column = row2Ptr.Col
			element2 = row2Ptr
			row2Ptr = row2Ptr.NextInRow
		default:
			column = row1Ptr.Col
			element1 = row1Ptr
			element2 = row2Ptr
			row1Ptr = row1Ptr.NextInRow
			row2Ptr = row2Ptr.NextInRow
		}

		m.exchangeColElements(row1, element1, row2, element2, column)
	}

	m.MarkowitzRow[row1], m.MarkowitzRow[row2] = m.MarkowitzRow[row2], m.MarkowitzRow[row1]
	m.FirstInRow[row1], m.FirstInRow[row2] = m.FirstInRow[row2], m.FirstInRow[row1]
	m.IntToExtRowMap[row1], m.IntToExtRowMap[row2] = m.IntToExtRowMap[row2], m.IntToExtRowMap[row1]

	if m.Config.Translate {
		m.ExtToIntRowMap[m.IntToExtRowMap[row1]] = row1
		m.ExtToIntRowMap[m.IntToExtRowMap[row2]] = row2
	}
}

func (m *Matrix) colExchange(col1, col2 int64) {
	if col1 > col2 {
		col1, col2 = col2, col1
	}

	col1Ptr := m.FirstInCol[col1]
	col2Ptr := m.FirstInCol[col2]

	for col1Ptr != nil || col2Ptr != nil {
		var row int64
		var element1, element2 *Element

		switch {
		case col1Ptr == nil:
			row = col2Ptr.Row
			element2 = col2Ptr
			col2Ptr = col2Ptr.NextInCol
		case col2Ptr == nil:
			row = col1Ptr.Row
			element1 = col1Ptr
			col1Ptr = col1Ptr.NextInCol
		case col1Ptr.Row < col2Ptr.Row:
			row = col1Ptr.Row
			element1 = col1Ptr
			col1Ptr = col1Ptr.NextInCol
		case col1Ptr.Row > col2Ptr.Row:
			row = col2Ptr.Row
			element2 = col2Ptr
			col2Ptr = col2Ptr.NextInCol
		default:
			row = col1Ptr.Row
			element1 = col1Ptr
			element2 = col2Ptr
			col1Ptr = col1Ptr.NextInCol
			col2Ptr = col2Ptr.NextInCol
		}

		m.exchangeRowElements(col1, element1, col2, element2, row)
	}

	m.MarkowitzCol[col1], m.MarkowitzCol[col2] = m.MarkowitzCol[col2], m.MarkowitzCol[col1]
	m.FirstInCol[col1], m.FirstInCol[col2] = m.FirstInCol[col2], m.FirstInCol[col1]
	m.IntToExtColMap[col1], m.IntToExtColMap[col2] = m.IntToExtColMap[col2], m.IntToExtColMap[col1]

	if m.Config.Translate {
		m.ExtToIntColMap[m.IntToExtColMap[col1]] = col1
		m.ExtToIntColMap[m.IntToExtColMap[col2]] = col2
	}
}

func (m *Matrix) exchangeColElements(row1 int64, element1 *Element, row2 int64, element2 *Element, column int64) {
	var elementAboveRow1, elementAboveRow2 **Element
	var elementBelowRow1, elementBelowRow2 *Element

	elementAboveRow1 = &m.FirstInCol[column]
	pElement := *elementAboveRow1
	for pElement.Row < row1 {
		elementAboveRow1 = &pElement.NextInCol
		pElement = *elementAboveRow1
	}

	if element1 != nil {
		elementBelowRow1 = element1.NextInCol
		if element2 == nil {
			if elementBelowRow1 != nil && elementBelowRow1.Row < row2 {
				*elementAboveRow1 = elementBelowRow1

				pElement = elementBelowRow1
				for pElement != nil && pElement.Row < row2 {
					elementAboveRow2 = &pElement.NextInCol
					pElement = *elementAboveRow2
				}

				*elementAboveRow2 = element1
				element1.NextInCol = pElement
			}
			element1.Row = row2
		} else {
			if elementBelowRow1.Row == row2 {
				element1.NextInCol = element2.NextInCol
				element2.NextInCol = element1
				*elementAboveRow1 = element2
			} else {
				pElement = elementBelowRow1
				for pElement.Row < row2 {
					elementAboveRow2 = &pElement.NextInCol
					pElement = *elementAboveRow2
				}

				elementBelowRow2 = element2.NextInCol

				*elementAboveRow1 = element2
				element2.NextInCol = elementBelowRow1
				*elementAboveRow2 = element1
				element1.NextInCol = elementBelowRow2
			}
			element1.Row = row2
			element2.Row = row1
		}
	} else {
		elementBelowRow1 = pElement

		if elementBelowRow1.Row != row2 {
			for pElement.Row < row2 {
				elementAboveRow2 = &pElement.NextInCol
				pElement = *elementAboveRow2
			}

			elementBelowRow2 = element2.NextInCol

			*elementAboveRow2 = elementBelowRow2
			*elementAboveRow1 = element2
			element2.NextInCol = elementBelowRow1
		}
		element2.Row = row1
	}
}

func (m *Matrix) exchangeRowElements(col1 int64, element1 *Element, col2 int64, element2 *Element, row int64) {
	elementLeftOfCol1 := &m.FirstInRow[row]
	pElement := *elementLeftOfCol1
	for pElement.Col < col1 {
		elementLeftOfCol1 = &pElement.NextInRow
		pElement = *elementLeftOfCol1
	}

	if element1 != nil {
		elementRightOfCol1 := element1.NextInRow
		if element2 == nil {
			if elementRightOfCol1 != nil && elementRightOfCol1.Col < col2 {
				*elementLeftOfCol1 = elementRightOfCol1

				pElement = elementRightOfCol1
				var elementLeftOfCol2 **Element
				for pElement != nil && pElement.Col < col2 {
					elementLeftOfCol2 = &pElement.NextInRow
					pElement = *elementLeftOfCol2
				}

				*elementLeftOfCol2 = element1
				element1.NextInRow = pElement
			}
			element1.Col = col2
		} else {
			if elementRightOfCol1.Col == col2 {
				element1.NextInRow = element2.NextInRow
				element2.NextInRow = element1
				*elementLeftOfCol1 = element2
			} else {
				pElement = elementRightOfCol1
				var elementLeftOfCol2 **Element
				for pElement.Col < col2 {
					elementLeftOfCol2 = &pElement.NextInRow
					pElement = *elementLeftOfCol2
				}

				elementRightOfCol2 := element2.NextInRow

				*elementLeftOfCol1 = element2
				element2.NextInRow = elementRightOfCol1
				*elementLeftOfCol2 = element1
				element1.NextInRow = elementRightOfCol2
			}
			element1.Col = col2
			element2.Col = col1
		}
	} else {
		elementRightOfCol1 := pElement

		if elementRightOfCol1.Col != col2 {
			var elementLeftOfCol2 **Element
			for pElement.Col < col2 {
				elementLeftOfCol2 = &pElement.NextInRow
				pElement = *elementLeftOfCol2
			}

			elementRightOfCol2 := element2.NextInRow

			*elementLeftOfCol2 = elementRightOfCol2
			*elementLeftOfCol1 = element2
			element2.NextInRow = elementRightOfCol1
		}
		element2.Col = col1
	}
}
