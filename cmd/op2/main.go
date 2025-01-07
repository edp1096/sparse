package main

import (
	"fmt"
	"log"

	"github.com/edp1096/sparse"
)

/* Modified nodal analysis */

func main() {
	const (
		R1  = 1000.0
		R2  = 1000.0
		Vin = 5.0
	)
	G1 := 1.0 / R1
	G2 := 1.0 / R2

	annotate := 0
	separatedComplexVectors := false

	config := &sparse.Configuration{
		Real:                    true,
		Complex:                 false,
		SeparatedComplexVectors: separatedComplexVectors,
		Expandable:              true,
		Translate:               false,
		ModifiedNodal:           true,
		TiesMultiplier:          5,
		PrinterWidth:            140,
		Annotate:                annotate,
	}

	A, err := sparse.Create(3, config) // 3x3 matrix (2 nodes + 1 branch)
	if err != nil {
		log.Fatalf("Failed to create matrix: %v", err)
	}
	defer A.Destroy()

	A.Clear()
	A.GetElement(1, 1).Real += G1
	A.GetElement(1, 2).Real += -G1
	A.GetElement(2, 1).Real += -G1
	A.GetElement(2, 2).Real += G1 + G2

	// Branch stamp
	A.GetElement(3, 1).Real = 1.0 // v1 coefficient
	A.GetElement(1, 3).Real = 1.0 // KCL current term

	fmt.Println("Matrix before factorization:")
	A.Print(false, true, true)

	b := make([]float64, A.Size+1)
	b[1] = 0
	b[2] = 0
	b[3] = Vin

	err = A.Factor()
	if err != nil {
		log.Fatalf("Failed to factor matrix: %v", err)
	}

	fmt.Println("\nMatrix after factorization:")
	A.Print(false, true, true)

	x, err := A.Solve(b)
	if err != nil {
		log.Fatalf("Failed to solve matrix: %v", err)
	}

	fmt.Printf("\nResults:\n")
	fmt.Printf("V1 = %.4f V\n", x[1])
	fmt.Printf("V2 = %.4f V\n", x[2])
	fmt.Printf("I_vs = %.4f mA\n", x[3]*1000) // Branch current
}
