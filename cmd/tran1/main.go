package main

import (
	"fmt"
	"math"
	"sparse"
)

func main() {
	const (
		R1    = 1000.0 // First resistor: 1k
		R2    = 1000.0 // Second resistor: 1k
		Vpeak = 5.0    // Peak voltage: 5V
		freq  = 1000.0 // Frequency: 1kHz
	)

	// 3x3 matrix for MNA (2 nodes + 1 voltage source)
	config := &sparse.Configuration{
		Real:          true,
		Complex:       false,
		Expandable:    true,
		Translate:     true,
		ModifiedNodal: true,
	}

	A, err := sparse.Create(3, config)
	if err != nil {
		panic(fmt.Sprintf("Failed to create matrix: %v", err))
	}
	defer A.Destroy()

	startTime := 0.0
	endTime := 0.002    // 2ms
	timeStep := 0.00001 // 10us

	fmt.Println("Time (s) | Vin | V(1) | V(2) | I_source")
	fmt.Println("------------------------------------------")

	for t := startTime; t <= endTime; t += timeStep {
		A.Clear()

		G1 := 1.0 / R1
		G2 := 1.0 / R2

		A.GetElement(1, 1).Real += G1
		A.GetElement(1, 2).Real -= G1
		A.GetElement(2, 1).Real -= G1
		A.GetElement(2, 2).Real += G1

		A.GetElement(2, 2).Real += G2

		A.GetElement(1, 3).Real = 1.0 // KCL for node 1
		A.GetElement(3, 1).Real = 1.0 // KVL

		radian := 2.0 * math.Pi * freq * t
		vin := Vpeak * math.Sin(radian)

		// RHS
		b := make([]float64, A.Size+1)
		b[3] = vin // Voltage source

		err = A.Factor()
		if err != nil {
			fmt.Printf("Time %.6f: Factorization failed - %v\n", t, err)
			continue
		}

		x, err := A.Solve(b)
		if err != nil {
			fmt.Printf("Time %.6f: Solve failed - %v\n", t, err)
			continue
		}

		fmt.Printf("%.6f | %7.3f | %7.3f | %7.3f | %7.3f\n", t, vin, x[1], x[2], x[3])
	}

	fmt.Println("\nCircuit Parameters:")
	fmt.Printf("R1: %.0f Ohm\n", R1)
	fmt.Printf("R2: %.0f Ohm\n", R2)
	fmt.Printf("Vpeak: %.1f V\n", Vpeak)
	fmt.Printf("Frequency: %.0f Hz\n", freq)
	fmt.Printf("Expected voltage division ratio: %.3f\n", R2/(R1+R2))
}
