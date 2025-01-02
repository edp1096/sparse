package main

import (
	"fmt"
	"log"
	"sparse"
)

func main() {
	var err error

	annotate := 0
	separatedComplexVectors := false

	config := &sparse.Configuration{
		Real:                    false,
		Complex:                 true,
		SeparatedComplexVectors: separatedComplexVectors,
		Expandable:              true,
		Translate:               true,
		ModifiedNodal:           true,
		TiesMultiplier:          5,
		PrinterWidth:            140,
		Annotate:                annotate,
	}

	A, err := sparse.Create(3, config)
	if err != nil {
		panic(err)
	}

	A.Clear()

	A.GetElement(1, 1).Real = 2
	A.GetElement(1, 1).Imag = 1
	A.GetElement(1, 2).Real = 1
	A.GetElement(1, 2).Imag = -1
	A.GetElement(1, 3).Real = 3
	A.GetElement(1, 3).Imag = 2

	A.GetElement(2, 1).Real = -1
	A.GetElement(2, 1).Imag = 1
	A.GetElement(2, 2).Real = 3
	A.GetElement(2, 2).Imag = 0
	A.GetElement(2, 3).Real = 2
	A.GetElement(2, 3).Imag = -1

	A.GetElement(3, 1).Real = 1
	A.GetElement(3, 1).Imag = -2
	A.GetElement(3, 2).Real = 2
	A.GetElement(3, 2).Imag = 1
	A.GetElement(3, 3).Real = 4
	A.GetElement(3, 3).Imag = 0

	fmt.Println("Original Complex Matrix:")
	A.Print(false, true, true)

	err = A.Factor()
	if err != nil {
		panic(err)
	}

	fmt.Println("\nFactored Complex Matrix:")
	A.Print(true, true, true)

	fmt.Println("\nSolution complex vector x:")
	if A.Config.SeparatedComplexVectors {
		b := make([]float64, 5)
		b[0] = 0.0
		b[1] = 3.0
		b[2] = 1.0
		b[3] = 2.0

		ib := make([]float64, 5)
		ib[0] = 0.0
		ib[1] = 2.0
		ib[2] = -1.0
		ib[3] = 1.0

		x, ix, err := A.SolveComplex(b, ib)
		if err != nil {
			panic(err)
		}

		fmt.Println("\nRHS complex vector b:")
		for i := 1; i <= 3; i++ {
			fmt.Printf("b[%d] = %.4f + %.4fj\n", i, b[i], ib[i])
		}

		log.Println(x)
		fmt.Println("\nSolution complex vector x:")
		for i := 1; i <= 3; i++ {
			fmt.Printf("x[%d] = %.4f + %.4fj\n", i, x[i], ix[i])
		}
	} else {
		b := make([]float64, 8) // 2*(size+1) for complex values
		b[2] = 3.0              // real part
		b[3] = 2.0              // imaginary part
		b[4] = 1.0
		b[5] = -1.0
		b[6] = 2.0
		b[7] = 1.0

		fmt.Println("\nRHS complex vector b:")
		for i := 1; i <= 3; i++ {
			fmt.Printf("b[%d] = %.4f + %.4fj\n", i, b[2*i], b[2*i+1])
		}

		x, _, err := A.SolveComplex(b, nil)
		if err != nil {
			panic(err)
		}

		fmt.Println("\nSolution complex vector x:")
		for i := 1; i <= 3; i++ {
			fmt.Printf("x[%d] = %.4f + %.4fj\n", i, x[2*i], x[2*i+1])
		}
	}

	A.Destroy()
}
