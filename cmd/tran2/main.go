package main

import (
	"fmt"
	"math"
	"sparse"
)

func main() {
	const (
		R     = 1000.0 // Resistor: 1k
		C     = 100e-6 // Capacitor: 100uF
		Vpeak = 5.0    // Peak voltage: 5V
		freq  = 1000.0 // Frequency: 1kHz
	)

	config := &sparse.Configuration{
		Real:          true,
		Complex:       false,
		Expandable:    true,
		Translate:     true,
		ModifiedNodal: true,
	}

	// 3x3 matrix for MNA (2 nodes + 1 voltage source)
	A, err := sparse.Create(3, config)
	if err != nil {
		panic(fmt.Sprintf("Failed to create matrix: %v", err))
	}
	defer A.Destroy()

	startTime := 0.0
	endTime := 1.0
	timeStep := 0.0001

	vcap := 0.0 // Initial condition - capacitor voltage

	fmt.Println("Time (s) | Vin | V(1) | V(2) | I_source")
	fmt.Println("------------------------------------------")

	for t := startTime; t <= endTime; t += timeStep {
		A.Clear()

		G := 1.0 / R
		Gc := C / timeStep

		A.GetElement(1, 1).Real += G
		A.GetElement(1, 2).Real -= G
		A.GetElement(2, 1).Real -= G
		A.GetElement(2, 2).Real += G

		A.GetElement(2, 2).Real += Gc

		A.GetElement(1, 3).Real = 1.0
		A.GetElement(3, 1).Real = 1.0

		radian := 2.0 * math.Pi * freq * t
		vin := Vpeak * math.Abs(math.Sin(radian))

		// RHS
		b := make([]float64, A.Size+1)
		b[2] = Gc * vcap // Capacitor history term
		b[3] = vin       // Voltage source

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

		vcap = x[2] // Store capacitor voltage for next iteration

		fmt.Printf("%.6f | %7.3f | %7.3f | %7.3f | %7.3f\n", t, vin, x[1], x[2], x[3])
	}

	fmt.Println("\nCircuit Parameters:")
	fmt.Printf("R: %.0f Ohm\n", R)
	fmt.Printf("C: %.0f uF\n", C*1e6)
	fmt.Printf("Vpeak: %.1f V\n", Vpeak)
	fmt.Printf("Frequency: %.0f Hz\n", freq)
	fmt.Printf("Time constant (RC): %.3f ms\n", R*C*1000)
}
