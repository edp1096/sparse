package main

import (
	"fmt"
	"log"

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

func main() {
	// 매트릭스 크기 및 구성 설정
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

	// 매트릭스 생성
	matrix, err := sparse.Create(int64(matrixSize), config)
	if err != nil {
		log.Fatalf("Failed to create matrix: %v", err)
	}
	defer matrix.Destroy()

	// 매트릭스 초기화
	matrix.Clear()

	// 회로 매트릭스 초기화
	circuitMatrix := &CircuitMatrix{Size: matrixSize, matrix: matrix}

	// 예제 요소 추가
	circuitMatrix.AddComplexElement(1, 1, 1.0, 0.5)
	circuitMatrix.AddComplexElement(2, 2, -1.0, -0.5)

	// 매트릭스 상태 출력
	fmt.Println("Matrix after adding elements:")
	for i := 1; i <= matrixSize; i++ {
		for j := 1; j <= matrixSize; j++ {
			element := circuitMatrix.matrix.GetElement(int64(i), int64(j))
			if element != nil {
				fmt.Printf("Element (%d, %d): %v + %vj\n", i, j, element.Real, element.Imag)
			} else {
				fmt.Printf("Element (%d, %d): nil\n", i, j)
			}
		}
	}
}
