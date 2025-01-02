package main

import (
	"fmt"
	"sparse"
)

func main() {
	var err error

	annotate := 0
	separatedComplexVectors := false

	config := &sparse.Configuration{
		Real:                    true,
		Complex:                 false,
		SeparatedComplexVectors: separatedComplexVectors,
		Expandable:              true,
		Translate:               true,
		ModifiedNodal:           true,
		TiesMultiplier:          5,
		PrinterWidth:            140,
		Annotate:                annotate,
	}

	A, err := sparse.Create(5, config)
	if err != nil {
		panic(err)
	}

	A.Clear()
	A.GetElement(1, 1).Real += 4
	A.GetElement(1, 2).Real += -2
	A.GetElement(1, 3).Real += 2
	A.GetElement(1, 4).Real += 1
	A.GetElement(1, 5).Real += 5

	A.GetElement(2, 1).Real += 2
	A.GetElement(2, 2).Real += 3
	A.GetElement(2, 3).Real += -1
	A.GetElement(2, 4).Real += 2
	A.GetElement(2, 5).Real += 3

	A.GetElement(3, 2).Real += 1
	A.GetElement(3, 3).Real += 5
	A.GetElement(3, 4).Real += 7
	A.GetElement(3, 5).Real += 2

	A.GetElement(4, 1).Real += 1
	A.GetElement(4, 2).Real += 2
	A.GetElement(4, 4).Real += 4
	A.GetElement(4, 5).Real += 1

	A.GetElement(5, 1).Real += 3
	A.GetElement(5, 2).Real += 1
	A.GetElement(5, 3).Real += 4
	A.GetElement(5, 4).Real += 2
	A.GetElement(5, 5).Real += 2

	A.Print(false, true, true)

	err = A.Factor()
	if err != nil {
		panic(err)
	}

	A.Print(true, true, true)

	b := make([]float64, 6)
	b[1] = 5.0
	b[2] = 0.0
	b[3] = 0.0
	b[4] = 0.0
	b[5] = 0.0

	fmt.Println("RHS b:")
	for i := 1; i <= 5; i++ {
		fmt.Printf("b[%d] = %.4f\n", i, b[i])
	}

	x, err := A.Solve(b)
	if err != nil {
		panic(err)
	}

	fmt.Println("Solution x:")
	for i := 1; i <= 5; i++ {
		fmt.Printf("x[%d] = %.4f\n", i, x[i])
	}

	Vin := 5.0
	fmt.Printf("Voltage at node 1 (V1): %.4f V\n", x[1])
	fmt.Printf("Voltage at node 2 (V2, output): %.4f V\n", x[2])
	fmt.Printf("Division ratio (R2 / (R1 + R2)): %.4f\n", x[2]/Vin)

	A.Destroy()
}
