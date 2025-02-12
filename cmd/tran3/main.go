/*
Simple RL Circuit
*/

package main

import (
	"fmt"
	"math"

	"github.com/edp1096/sparse"
)

type IntegrationMethod int

type BackwardDifferentialFormula struct {
	coefficients []float64 // BDF coefficients
	beta         float64   // BDF f(y^(n+1)) coefficient
}

const (
	R     = 100.0
	L     = 1e-3
	Vpeak = 5.0
	freq  = 1000.0

	G = 1.0 / R

	tstop    = 0.002
	timestep = 1e-5

	integrationMethod = TrapezoidalMethod
	// integrationMethod = GearMethod
	methodOrder = 6
)

const (
	GearMethod IntegrationMethod = iota
	TrapezoidalMethod
)

// BDF coefficients
var bdfCoefficients = [6]BackwardDifferentialFormula{
	{[]float64{1.0}, 1.0},
	{[]float64{4.0 / 3.0, -1.0 / 3.0}, 2.0 / 3.0},
	{[]float64{18.0 / 11.0, -9.0 / 11.0, 2.0 / 11.0}, 6.0 / 11.0},
	{[]float64{48.0 / 25.0, -36.0 / 25.0, 16.0 / 25.0, -3.0 / 25.0}, 12.0 / 25.0},
	{[]float64{300.0 / 137.0, -300.0 / 137.0, 200.0 / 137.0, -75.0 / 137.0, 12.0 / 137.0}, 60.0 / 137.0},
	{[]float64{360.0 / 147.0, -450.0 / 147.0, 400.0 / 147.0, -225.0 / 147.0, 72.0 / 147.0, -10.0 / 147.0}, 60.0 / 147.0},
}

func GetBDFcoeffs(order int, dt float64) []float64 {
	bdf := bdfCoefficients[order-1]
	coeffs := make([]float64, order+1)
	scale := 1.0 / (bdf.beta * dt)
	coeffs[0] = scale

	for i := 1; i <= order; i++ {
		coeffs[i] = -bdf.coefficients[i-1] * scale
	}

	return coeffs
}

func GetTrapezoidalCoeffs(order int, dt float64) []float64 {
	coeffs := make([]float64, 1)

	coeffs[0] = 2.0 / dt
	if order == 1 {
		coeffs[0] = 1.0 / dt
	}

	return coeffs
}

func GetIntegratorCoeffs(order int, dt float64) []float64 {
	switch integrationMethod {
	case TrapezoidalMethod:
		return GetTrapezoidalCoeffs(order, dt)
	default:
		return GetBDFcoeffs(order, dt)
	}
}

func main() {
	t0 := 0.0
	dt := timestep
	N := int(tstop/dt) + 1 // data count

	config := &sparse.Configuration{
		Real:          true,
		Complex:       false,
		ModifiedNodal: true,
		Translate:     true,
	}

	iL := make([]float64, N)
	vL := make([]float64, N)

	fmt.Println("Time (s) |   iL (A)   |  vL (V)")
	fmt.Println("-----------------------------------")

	A, err := sparse.Create(4, config)
	if err != nil {
		panic(err)
	}
	defer A.Destroy()

	b := make([]float64, A.Size+1) // RHS

	// DC operating point
	A.Clear()

	A.GetElement(1, 1).Real += G
	A.GetElement(1, 3).Real += 1.0
	A.GetElement(2, 3).Real += 1.0
	A.GetElement(2, 4).Real += 1.0
	A.GetElement(3, 2).Real += 1.0
	A.GetElement(4, 1).Real += -1.0
	A.GetElement(4, 2).Real += 1.0

	b[1] = 0.0
	b[2] = 0.0
	b[3] = Vpeak * math.Sin(2.0*math.Pi*freq*t0)
	b[4] = 0.0

	A.MNAPreorder()
	A.Factor()

	x, err := A.Solve(b)
	if err != nil {
		panic(err)
	}

	iL[0] = x[3]
	vL[0] = x[2] - x[1]

	fmt.Printf("%.6f | %.8f | %.8f\n", t0, iL[0], vL[0])

	// Transient
	for i := 0; i < N-1; i++ {
		order := methodOrder
		if i+1 < methodOrder {
			order = i + 1
		}

		coeffs := GetIntegratorCoeffs(order, dt)
		tNext := float64(i+1) * dt

		A.Clear()

		A.GetElement(1, 1).Real += G
		A.GetElement(1, 3).Real += 1.0
		A.GetElement(2, 3).Real += 1.0
		A.GetElement(2, 4).Real += 1.0
		A.GetElement(3, 2).Real += 1.0
		A.GetElement(4, 1).Real += -1.0
		A.GetElement(4, 2).Real += 1.0
		A.GetElement(4, 3).Real += coeffs[0] * L

		b[1] = 0.0
		b[2] = 0.0
		b[3] = Vpeak * math.Sin(2.0*math.Pi*freq*tNext)
		b[4] = 0.0
		switch integrationMethod {
		case TrapezoidalMethod:
			if i < 2 || order == 1 {
				b[4] = coeffs[0] * iL[i] * L // order1 (Backward Euler)
			} else {
				b[4] = coeffs[0]*iL[i]*L - vL[i] // order2 (현재 잘 동작하는 코드)
			}

		default:
			for j := 1; j <= order; j++ {
				b[4] -= coeffs[j] * iL[i+1-j]
			}
			b[4] *= L
		}

		A.MNAPreorder()
		A.Factor()

		x, err = A.Solve(b)
		if err != nil {
			panic(err)
		}

		iL[i+1] = x[3]
		vL[i+1] = x[2] - x[1]

		if (i+1)%(N/20) == 0 {
			fmt.Printf("%.6f | %.7f | %.7f\n", tNext, iL[i+1], vL[i+1])
		}
	}

	maxVL := 0.0
	for _, v := range vL {
		if math.Abs(v) > maxVL {
			maxVL = math.Abs(v)
		}
	}
	maxDiDt := Vpeak * 2 * math.Pi * freq / math.Sqrt(R*R+math.Pow(2*math.Pi*freq*L, 2))
	theoryMaxVL := L * maxDiDt

	fmt.Println()
	if integrationMethod == TrapezoidalMethod {
		fmt.Printf("Trapezoidal\n")
	} else {
		fmt.Printf("Gear%d\n", methodOrder)
	}
	fmt.Printf("Theory max VL: %.6f V\n", theoryMaxVL)
	fmt.Printf("Max VL: %.6f V\n", maxVL)
	errPct := 100 * math.Abs(maxVL-theoryMaxVL) / theoryMaxVL
	fmt.Printf("Err: %.6f%%\n", errPct)
}
