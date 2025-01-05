package sparse

import (
	"fmt"
	"math"
)

func (m *Matrix) WriteStatus(step int64) {
	fmt.Printf("Step = %d   ", step)
	fmt.Printf("Pivot found at %d,%d using ", m.PivotsOriginalRow, m.PivotsOriginalCol)

	switch m.PivotSelectionMethod {
	case 's':
		fmt.Println("SearchForSingleton")
	case 'q':
		fmt.Println("QuicklySearchDiagonal")
	case 'd':
		fmt.Println("SearchDiagonal")
	case 'e':
		fmt.Println("SearchEntireMatrix")
	}

	// Markowitz information
	fmt.Printf("MarkowitzRow     = ")
	for i := int64(1); i <= m.Size; i++ {
		fmt.Printf("%2d  ", m.MarkowitzRow[i])
	}
	fmt.Println()

	fmt.Printf("MarkowitzCol     = ")
	for i := int64(1); i <= m.Size; i++ {
		fmt.Printf("%2d  ", m.MarkowitzCol[i])
	}
	fmt.Println()

	fmt.Printf("MarkowitzProduct = ")
	for i := int64(1); i <= m.Size; i++ {
		fmt.Printf("%2d  ", m.MarkowitzProd[i])
	}
	fmt.Println()

	fmt.Printf("Singletons = %2d\n", m.Singletons)

	// Mapping information
	fmt.Printf("IntToExtRowMap     = ")
	for i := int64(1); i <= m.Size; i++ {
		fmt.Printf("%2d  ", m.IntToExtRowMap[i])
	}
	fmt.Println()

	fmt.Printf("IntToExtColMap     = ")
	for i := int64(1); i <= m.Size; i++ {
		fmt.Printf("%2d  ", m.IntToExtColMap[i])
	}
	fmt.Println()

	fmt.Printf("ExtToIntRowMap     = ")
	for i := int64(1); i <= m.Size; i++ {
		fmt.Printf("%2d  ", m.ExtToIntRowMap[i])
	}
	fmt.Println()

	fmt.Printf("ExtToIntColMap     = ")
	for i := int64(1); i <= m.Size; i++ {
		fmt.Printf("%2d  ", m.ExtToIntColMap[i])
	}
	fmt.Printf("\n\n")
}

func (m *Matrix) Print(printReordered bool, data bool, header bool) {
	if m == nil {
		return
	}

	top := m.Size
	if m.Config.Translate {
		top = m.ExtSize
	}
	printOrdToIntRowMap := make([]int64, top+1)
	printOrdToIntColMap := make([]int64, top+1)

	for i := int64(1); i <= m.Size; i++ {
		printOrdToIntRowMap[m.IntToExtRowMap[i]] = i
		printOrdToIntColMap[m.IntToExtColMap[i]] = i
	}

	compressMap := func(cmap []int64) []int64 {
		compressed := make([]int64, 1, top)
		for i := int64(1); i <= top; i++ {
			if cmap[i] != 0 {
				compressed = append(compressed, cmap[i])
			}
		}
		return compressed
	}

	printOrdToIntRowMap = compressMap(printOrdToIntRowMap)
	printOrdToIntColMap = compressMap(printOrdToIntColMap)

	if header {
		fmt.Printf("MATRIX SUMMARY\n\n")
		fmt.Printf("Size of matrix = %d x %d.\n", m.Size, m.Size)
		if m.Reordered && printReordered {
			fmt.Printf("Matrix has been reordered.\n")
		}
		fmt.Println()

		if m.Factored {
			fmt.Printf("Matrix after factorization:\n")
		} else {
			fmt.Printf("Matrix before factorization:\n")
		}
	}

	if m.Size == 0 {
		return
	}

	columns := m.Config.PrinterWidth
	if header {
		columns -= 5
	}
	if data {
		columns = (columns + 1) / 10
	}

	getRow := func(i int64) int64 {
		if printReordered {
			return i
		}
		return printOrdToIntRowMap[i]
	}

	getCol := func(j int64) int64 {
		if printReordered {
			return j
		}
		return printOrdToIntColMap[j]
	}

	startCol := int64(1)
	for startCol <= m.Size {
		stopCol := startCol + int64(columns) - 1
		if stopCol > m.Size {
			stopCol = m.Size
		}

		if header {
			if data {
				fmt.Printf("    ")
				for col := startCol; col <= stopCol; col++ {
					actualCol := getCol(col)
					fmt.Printf(" %9d", m.IntToExtColMap[actualCol])
				}
				fmt.Printf("\n\n")
			} else {
				firstCol := getCol(startCol)
				lastCol := getCol(stopCol)
				if printReordered {
					fmt.Printf("Columns %d to %d.\n", startCol, stopCol)
				} else {
					fmt.Printf("Columns %d to %d.\n",
						m.IntToExtColMap[firstCol],
						m.IntToExtColMap[lastCol])
				}
			}
		}

		for i := int64(1); i <= m.Size; i++ {
			row := getRow(i)

			if header {
				if printReordered && !data {
					fmt.Printf("%4d", i)
				} else {
					fmt.Printf("%4d", m.IntToExtRowMap[row])
				}
				if !data {
					fmt.Printf(" ")
				}
			}

			imagElements := make([]*Element, stopCol-startCol+1)

			for colIndex := startCol; colIndex <= stopCol; colIndex++ {
				col := getCol(colIndex)

				var element *Element
				for e := m.FirstInCol[col]; e != nil; e = e.NextInCol {
					if e.Row == row {
						element = e
						break
					}
				}

				imagElements[colIndex-startCol] = element

				if element != nil {
					if data {
						fmt.Printf(" %9.3g", element.Real)
					} else {
						fmt.Printf("x")
					}
				} else {
					if data {
						fmt.Printf("       ...")
					} else {
						fmt.Printf(".")
					}
				}
			}
			fmt.Println()

			if data && m.Complex {
				if header {
					fmt.Printf("    ")
				}
				for _, element := range imagElements {
					if element != nil {
						fmt.Printf(" %8.2gj", element.Imag)
					} else {
						fmt.Printf("          ")
					}
				}
				fmt.Println()
			}
		}

		fmt.Println()
		startCol = stopCol + 1
	}

	if header {
		stats := m.calculateStatistics(top)
		fmt.Printf("\nLargest element in matrix = %-1.4g.\n", stats.largestElement)
		fmt.Printf("Smallest element in matrix = %-1.4g.\n", stats.smallestElement)

		if m.Factored {
			fmt.Printf("\nLargest diagonal element = %-1.4g.\n", stats.largestDiag)
			fmt.Printf("Smallest diagonal element = %-1.4g.\n", stats.smallestDiag)
		} else {
			fmt.Printf("\nLargest pivot element = %-1.4g.\n", stats.largestDiag)
			fmt.Printf("Smallest pivot element = %-1.4g.\n", stats.smallestDiag)
		}

		density := float64(stats.elementCount) * 100.0 / float64(m.Size*m.Size)
		fmt.Printf("\nDensity = %.2f%%.\n", density)
		if !m.NeedsOrdering {
			fmt.Printf("Number of fill-ins = %d.\n", m.Fillins)
		}
		fmt.Println()
	}
}

type matrixStats struct {
	largestElement  float64
	smallestElement float64
	largestDiag     float64
	smallestDiag    float64
	elementCount    int64
}

func (m *Matrix) calculateStatistics(top int64) matrixStats {
	stats := matrixStats{
		smallestElement: math.MaxFloat64,
		smallestDiag:    math.MaxFloat64,
	}

	elementCount := int64(0)

	if m.Size == 0 {
		return stats
	}

	printOrdToIntRowMap := make([]int64, top+1)
	printOrdToIntColMap := make([]int64, top+1)

	for i := int64(1); i <= m.Size; i++ {
		printOrdToIntRowMap[m.IntToExtRowMap[i]] = i
		printOrdToIntColMap[m.IntToExtColMap[i]] = i
	}

	for i := int64(1); i <= top; i++ {
		row := printOrdToIntRowMap[i]

		for col := int64(1); col <= top; col++ {
			actualCol := printOrdToIntColMap[col]

			element := m.FirstInCol[actualCol]
			for element != nil && element.Row != row {
				element = element.NextInCol
			}

			if element != nil {
				elementCount++
				magnitude := m.elementMag(element)
				if m.Complex {
					magnitude = m.elementMag(element)
				}

				if magnitude > stats.largestElement {
					stats.largestElement = magnitude
				}
				if magnitude < stats.smallestElement && magnitude != 0 {
					stats.smallestElement = magnitude
				}

				if row == actualCol {
					if magnitude > stats.largestDiag {
						stats.largestDiag = magnitude
					}
					if magnitude < stats.smallestDiag && magnitude != 0 {
						stats.smallestDiag = magnitude
					}
				}
			}
		}
	}

	if elementCount == 0 {
		stats.smallestElement = 0
		stats.largestElement = 0
		stats.smallestDiag = 0
		stats.largestDiag = 0
	}

	stats.elementCount = elementCount
	return stats
}
