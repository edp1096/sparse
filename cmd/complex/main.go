package main

import (
	"fmt"
	"log"
	"math"

	"github.com/edp1096/sparse"
)

type CircuitMatrix struct {
	Size   int
	matrix *sparse.Matrix
}

func (m *CircuitMatrix) AddComplexElement(i, j int, real, imag float64) {
	fmt.Printf("Adding element at (i=%d, j=%d) with real=%f, imag=%f\n", i, j, real, imag)
	if i < 1 || j < 1 || i > m.Size || j > m.Size {
		fmt.Printf("Warning: Matrix index out of bounds (i=%d, j=%d, size=%d)\n", i, j, m.Size)
		return
	}

	element := m.matrix.GetElement(int64(i), int64(j))
	if element == nil {
		log.Fatalf("Error: Element at (i=%d, j=%d) not found\n", i, j)
	}
	element.Real += real
	element.Imag += imag
}

func (m *CircuitMatrix) PrintMatrixWithPhase() {
	fmt.Println("Matrix with phase angles:")
	for i := 1; i <= m.Size; i++ {
		for j := 1; j <= m.Size; j++ {
			element := m.matrix.GetElement(int64(i), int64(j))
			if element != nil {
				phase := math.Atan2(element.Imag, element.Real) * (180 / math.Pi)
				fmt.Printf("Element (%d, %d): %v + %vj, Phase: %v degrees\n", i, j, element.Real, element.Imag, phase)
			} else {
				fmt.Printf("Element (%d, %d): nil\n", i, j)
			}
		}
	}
}

func main() {
	matrixSize := 5
	config := &sparse.Configuration{
		Real:                    true,
		Complex:                 true,
		SeparatedComplexVectors: false,
		Expandable:              true,
		Translate:               false,
		ModifiedNodal:           true,
		TiesMultiplier:          5,
		PrinterWidth:            140,
		Annotate:                0,
	}

	matrix, err := sparse.Create(int64(matrixSize), config)
	if err != nil {
		log.Fatalf("Failed to create matrix: %v", err)
	}
	defer matrix.Destroy()

	matrix.Clear()

	circuitMatrix := &CircuitMatrix{Size: matrixSize, matrix: matrix}

	circuitMatrix.AddComplexElement(1, 1, 1.0, 0.5)
	circuitMatrix.AddComplexElement(2, 2, -1.0, -0.5)

	circuitMatrix.PrintMatrixWithPhase()
}
