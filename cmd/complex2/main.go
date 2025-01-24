package main

import (
	"fmt"
	"math"

	"github.com/edp1096/sparse"
)

func main() {
	var err error

	config := &sparse.Configuration{
		Real:                    true,
		Complex:                 true,
		SeparatedComplexVectors: true,
		Expandable:              true,
		Translate:               false,
		ModifiedNodal:           true,
		TiesMultiplier:          5,
		PrinterWidth:            140,
		Annotate:                0,
	}

	A, err := sparse.Create(2, config)
	if err != nil {
		panic(err)
	}
	defer A.Destroy()
	A.Clear()

	// 2x2 복소수 테스트 매트릭스 설정
	fmt.Println("Setting up complex test matrix...")

	// // Y11
	// A.GetElement(1, 1).Real = 1.0
	// A.GetElement(1, 1).Imag = 0.0

	// // Y12
	// A.GetElement(1, 2).Real = 0.0
	// A.GetElement(1, 2).Imag = 1.0

	// // Y21
	// A.GetElement(2, 1).Real = 0.0
	// A.GetElement(2, 1).Imag = -1.0

	// // Y22
	// A.GetElement(2, 2).Real = 1.0
	// A.GetElement(2, 2).Imag = 0.0

	// Y11
	A.GetElement(1, 1).Real = 2.0
	A.GetElement(1, 1).Imag = 1.0

	// Y12
	A.GetElement(1, 2).Real = -0.5
	A.GetElement(1, 2).Imag = 1.0

	// Y21
	A.GetElement(2, 1).Real = -0.5
	A.GetElement(2, 1).Imag = -1.0

	// Y22
	A.GetElement(2, 2).Real = 2.0
	A.GetElement(2, 2).Imag = -1.0

	fmt.Println("\nOriginal matrix:")
	A.Print(false, true, true)

	err = A.Factor()
	if err != nil {
		panic(err)
	}

	fmt.Println("\nFactored matrix:")
	A.Print(true, true, true)

	// RHS
	// SeparatedComplexVectors = true
	b := make([]float64, 3)
	ib := make([]float64, 3)
	b[1] = 1.0
	ib[1] = 0.0
	b[2] = 0.0
	ib[2] = 0.0

	// SeparatedComplexVectors = false
	sb := make([]float64, 6) // 2 nodes * (Real + Imag)
	sb[2] = 1.0              // node 1 real
	sb[3] = 0.0              // node 1 imag
	sb[4] = 0.0              // node 2 real
	sb[5] = 0.0              // node 2 imag

	fmt.Println("\nRHS vector:")
	if config.SeparatedComplexVectors {
		for i := 1; i <= 2; i++ {
			fmt.Printf("Node %d: %.6f + j%.6f\n", i+1, b[i], ib[i])
		}
	} else {
		for i := 1; i <= 2; i++ {
			fmt.Printf("Node %d: %.6f + j%.6f\n", i+1, sb[2*i], sb[2*i+1])
		}
	}

	var xReal, xImag []float64
	if config.SeparatedComplexVectors {
		xReal, xImag, err = A.SolveComplex(b, ib)
		if err != nil {
			panic(err)
		}
	} else {
		// x, err := A.Solve(sb)
		xReal, xImag, err = A.SolveComplex(sb, nil)
		if err != nil {
			panic(err)
		}
	}

	fmt.Println("\nSolution vector:")
	for i := 1; i <= 2; i++ {
		phase := math.Atan2(xImag[i], xReal[i]) * 180 / math.Pi
		fmt.Printf("Node %d: %.6f + j%.6f, Phase: %.2f degrees\n", i, xReal[i], xImag[i], phase)
	}

}
