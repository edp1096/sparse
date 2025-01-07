package main // import "op"

import (
	"fmt"
	"log"

	"github.com/edp1096/sparse"
)

func main() {
	const (
		R1  = 1000.0
		R2  = 2000.0
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

	A, err := sparse.Create(2, config)
	if err != nil {
		log.Fatalf("Failed to create matrix: %v", err)
	}
	defer A.Destroy()

	A.Clear()
	A.GetElement(1, 1).Real += G1
	A.GetElement(1, 2).Real += -G1
	A.GetElement(2, 1).Real += -G1
	A.GetElement(2, 2).Real += G1 + G2

	fmt.Println("Matrix before factorization:")
	A.Print(false, true, true)

	b := make([]float64, A.Size+1)
	b[1] = (1 / (R1 + R2)) * Vin
	b[2] = 0.0

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

	fmt.Printf("Voltage at node 1 (V1): %.4f V\n", x[1])
	fmt.Printf("Voltage at node 2 (V2, output): %.4f V\n", x[2])
	fmt.Printf("Division ratio (R2 / (R1 + R2)): %.4f\n", x[2]/Vin)
}
