package sparse

func (m *Matrix) CountMarkowitz(rhs []float64, step int64) {
	/* Generate MarkowitzRow Count for each row. */
	for i := step; i <= m.Size; i++ {
		count := int64(-1)
		element := m.FirstInRow[i]

		for element != nil && element.Col < step {
			element = element.NextInRow
		}
		for element != nil {
			count++
			element = element.NextInRow
		}

		if rhs != nil {
			extRow := m.IntToExtRowMap[i]
			switch {
			case m.Config.SeparatedComplexVectors:
				if rhs[extRow] != 0.0 {
					count++
				}
			case m.Complex:
				if rhs[2*extRow] != 0.0 || rhs[2*extRow+1] != 0.0 {
					count++
				}
			default:
				if rhs[i] != 0.0 {
					count++
				}
			}
		}

		m.MarkowitzRow[i] = count
	}

	/* Generate the MarkowitzCol count for each column. */
	for i := step; i <= m.Size; i++ {
		count := int64(-1)
		element := m.FirstInCol[i]

		for element != nil && element.Row < step {
			element = element.NextInCol
		}
		for element != nil {
			count++
			element = element.NextInCol
		}

		m.MarkowitzCol[i] = count
	}
}

func (m *Matrix) MarkowitzProducts(step int64) {
	m.Singletons = 0

	for i := step; i <= m.Size; i++ {
		rowCount := m.MarkowitzRow[i]
		colCount := m.MarkowitzCol[i]

		m.MarkowitzProd[i] = markowitzProduct(rowCount, colCount)

		if m.MarkowitzProd[i] == 0 {
			m.Singletons++
		}
	}
}

func markowitzProduct(op1, op2 int64) int64 {
	LargestShortInteger := int64(32767)
	LargestLongInteger := int64(2147483647)

	if (op1 > LargestShortInteger && op2 != 0) || (op2 > LargestShortInteger && op1 != 0) {
		fProduct := float64(op1) * float64(op2)
		if fProduct >= float64(LargestLongInteger) {
			return LargestLongInteger
		}
		return int64(fProduct)
	}
	return op1 * op2
}

func (m *Matrix) calculateMarkowitzProduct(row, col int64) int64 {
	const (
		LargestShortInteger = 32767
		LargestLongInteger  = 2147483647
	)

	if (row > LargestShortInteger && col != 0) || (col > LargestShortInteger && row != 0) {
		product := float64(row) * float64(col)
		if product >= float64(LargestLongInteger) {
			return LargestLongInteger
		}
		return int64(product)
	}
	return row * col
}

func (m *Matrix) UpdateMarkowitzNumbers(pivot *Element) {
	const (
		LargestShortInteger = 32767
		LargestLongInteger  = 2147483647
	)

	markowitzRow := m.MarkowitzRow
	markowitzCol := m.MarkowitzCol
	markowitzProd := m.MarkowitzProd

	// Column traversal
	colPtr := pivot.NextInCol
	for colPtr != nil {
		row := colPtr.Row
		markowitzRow[row]--

		if (markowitzRow[row] > LargestShortInteger && markowitzCol[row] != 0) || (markowitzCol[row] > LargestShortInteger && markowitzRow[row] != 0) {
			product := float64(markowitzCol[row]) * float64(markowitzRow[row])
			if product >= float64(LargestLongInteger) {
				markowitzProd[row] = LargestLongInteger
			} else {
				markowitzProd[row] = int64(product)
			}
		} else {
			markowitzProd[row] = markowitzRow[row] * markowitzCol[row]
		}

		if markowitzRow[row] == 0 {
			m.Singletons++
		}

		colPtr = colPtr.NextInCol
	}

	// Row traversal
	rowPtr := pivot.NextInRow
	for rowPtr != nil {
		col := rowPtr.Col
		markowitzCol[col]--

		if (markowitzRow[col] > LargestShortInteger && markowitzCol[col] != 0) || (markowitzCol[col] > LargestShortInteger && markowitzRow[col] != 0) {
			product := float64(markowitzCol[col]) * float64(markowitzRow[col])
			if product >= float64(LargestLongInteger) {
				markowitzProd[col] = LargestLongInteger
			} else {
				markowitzProd[col] = int64(product)
			}
		} else {
			markowitzProd[col] = markowitzRow[col] * markowitzCol[col]
		}

		if markowitzCol[col] == 0 && markowitzRow[col] != 0 {
			m.Singletons++
		}

		rowPtr = rowPtr.NextInRow
	}

	m.MarkowitzRow = markowitzRow
	m.MarkowitzCol = markowitzCol
	m.MarkowitzProd = markowitzProd
}
