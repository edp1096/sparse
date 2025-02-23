package main

import (
	"fmt"
	"math"

	"github.com/edp1096/sparse"
)

func main() {
	const (
		R     = 100.0   // Resistor: 100 ohm
		G     = 1.0 / R // Conductance: 1/100 ohm
		L     = 1e-3    // Inductor: 1mH
		Vpeak = 5.0     // Peak voltage: 5V
		freq  = 1000.0  // Frequency: 1kHz
	)

	config := &sparse.Configuration{
		Real:          true,
		Complex:       false,
		Expandable:    true,
		Translate:     true,
		ModifiedNodal: true,
	}

	A, err := sparse.Create(4, config)
	if err != nil {
		panic(fmt.Sprintf("Failed to create matrix: %v", err))
	}
	defer A.Destroy()

	startTime := 0.0
	endTime := 0.002
	timeStep := 0.00001

	iL := 0.0
	vL := make([]float64, 0)

	b := make([]float64, A.Size+1)

	fmt.Println("Time (s) | Vin | V(1) | V(2) | I_source | I_L")
	fmt.Println("---------------------------------------------------")

	for t := startTime; t <= endTime; t += timeStep {
		A.Clear()

		// Resistor
		A.GetElement(1, 1).Real += G
		A.GetElement(2, 2).Real += G
		A.GetElement(1, 2).Real += -G
		A.GetElement(2, 1).Real += -G

		// Voltage source
		A.GetElement(1, 3).Real += 1.0
		A.GetElement(3, 1).Real += 1.0

		// Inductor
		A.GetElement(2, 4).Real += 1.0
		A.GetElement(4, 2).Real += 1.0
		A.GetElement(4, 4).Real += -L / timeStep

		vin := Vpeak * math.Sin(2.0*math.Pi*freq*t)

		b[1] = 0.0
		b[2] = 0.0
		b[3] = vin
		b[4] = -(L / timeStep) * iL

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

		iL = x[4]
		vL = append(vL, x[2])
		fmt.Printf("%.6f | %7.5f | %7.5f | %7.5f | %7.5f | %7.5f\n", t, vin, x[1], x[2], x[3], x[4])
	}

	fmt.Println("\nCircuit Parameters:")
	fmt.Printf("R: %.0f Ohm\n", R)
	fmt.Printf("L: %.0f mH\n", L*1e3)
	fmt.Printf("Vpeak: %.1f V\n", Vpeak)
	fmt.Printf("Frequency: %.0f Hz\n", freq)
	fmt.Printf("Time constant (L/R): %.3f ms\n", (L/R)*1000)

	maxVL := 0.0
	for _, v := range vL {
		if math.Abs(v) > maxVL {
			maxVL = math.Abs(v)
		}
	}
	maxDiDt := Vpeak * 2 * math.Pi * freq / math.Sqrt(R*R+math.Pow(2*math.Pi*freq*L, 2))
	theoryMaxVL := L * maxDiDt

	fmt.Println()
	fmt.Printf("Theory max VL: %.6f V\n", theoryMaxVL)
	fmt.Printf("Max VL: %.6f V\n", maxVL)
	errPct := 100 * math.Abs(maxVL-theoryMaxVL) / theoryMaxVL
	fmt.Printf("Err: %.6f%%\n", errPct)
}
