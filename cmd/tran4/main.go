package main

import (
	"fmt"
	"math"

	"github.com/edp1096/sparse"
)

type IntegrationMethod int

const (
	GearMethod IntegrationMethod = iota
	TrapezoidalMethod
)

type BackwardDifferentialFormula struct {
	coefficients []float64
	beta         float64
}

const (
	R     = 100.0
	L     = 1e-3
	Vpeak = 5.0
	freq  = 1000.0
	G     = 1.0 / R

	tstop    = 0.002
	timestep = 1e-5

	targetLTE    = 1e-6
	minTimeStep  = 1e-9
	safetyFactor = 0.9
)

var (
	integrationMethod = GearMethod // TrapezoidalMethod or GearMethod
	methodOrder       = 6          // Trapezoidal: 1 or 2, Gear: 1 ~ 6

	// BDF coefficients
	bdfCoefficients = [6]BackwardDifferentialFormula{
		{[]float64{1.0}, 1.0},                                                                                                // 1st order
		{[]float64{4.0 / 3.0, -1.0 / 3.0}, 2.0 / 3.0},                                                                        // 2nd order
		{[]float64{18.0 / 11.0, -9.0 / 11.0, 2.0 / 11.0}, 6.0 / 11.0},                                                        // 3rd order
		{[]float64{48.0 / 25.0, -36.0 / 25.0, 16.0 / 25.0, -3.0 / 25.0}, 12.0 / 25.0},                                        // 4th order
		{[]float64{300.0 / 137.0, -300.0 / 137.0, 200.0 / 137.0, -75.0 / 137.0, 12.0 / 137.0}, 60.0 / 137.0},                 // 5th order
		{[]float64{360.0 / 147.0, -450.0 / 147.0, 400.0 / 147.0, -225.0 / 147.0, 72.0 / 147.0, -10.0 / 147.0}, 60.0 / 147.0}, // 6th order
	}

	lteTrapCoeffs = []float64{0.5, 1.0 / 12.0}
	lteBdfCoeffs  = []float64{0.5, 2.0 / 3.0, 6.0 / 11.0, 12.0 / 25.0, 60.0 / 137.0, 60.0 / 147.0}

	stabilityFactors = []float64{2.00, 2.00, 1.98, 1.92, 1.76, 1.56} // stability factor
)

func GetBDFcoeffs(order int, dt float64) []float64 {
	if order < 1 || order > 6 {
		panic(fmt.Sprintf("Invalid BDF order: %d", order))
	}

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
	if order < 1 || order > 2 {
		panic(fmt.Sprintf("Invalid Trapezoidal order: %d", order))
	}

	coeffs := make([]float64, 1)
	coeffs[0] = 2.0 / dt
	if order == 1 {
		coeffs[0] = 1.0 / dt // Backward Euler
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

func calculateLTE(values []float64, dt float64, order int) float64 {
	if len(values) < order+1 {
		return 0.0
	}

	// Derivatives values
	n := len(values)
	der := make([]float64, order+1)
	der[0] = values[n-1]

	for i := 1; i <= order; i++ {
		for j := 0; j < len(der)-i; j++ {
			der[j] = (der[j] - der[j+1]) / dt
		}
	}

	var lteCoeff float64
	switch integrationMethod {
	case TrapezoidalMethod:
		lteCoeff = lteTrapCoeffs[order-1]
	default:
		lteCoeff = lteBdfCoeffs[order-1]
	}

	return math.Abs(lteCoeff * math.Pow(dt, float64(order+1)) * der[0])
}

func calculateNewTimeStep(currentStep, lte float64, order int) float64 {
	if lte < 1e-15 {
		return currentStep
	}

	factor := math.Pow(targetLTE/lte, 1.0/float64(order+1)) * safetyFactor
	factor = math.Max(0.1, math.Min(factor, 10.0))

	newStep := currentStep * factor
	maxStep := stabilityFactors[order-1] * timestep

	newStep = math.Min(newStep, maxStep)
	newStep = math.Max(newStep, minTimeStep)

	return newStep
}

func checkMethodOrder() {
	if integrationMethod == TrapezoidalMethod && methodOrder > 2 {
		fmt.Println("Trapezoidal method order must be 1 or 2. Run as 2")
		methodOrder = 2
	}
	if integrationMethod == GearMethod && (methodOrder < 1 || methodOrder > 6) {
		fmt.Println("BDF/Gear method order must be 1 ~ 6. Run as 6")
		methodOrder = 6
	}
}

func main() {
	checkMethodOrder()

	config := &sparse.Configuration{
		Real:          true,
		Complex:       false,
		ModifiedNodal: true,
		Translate:     true,
	}

	t := 0.0
	dt := timestep
	timePoints := []float64{t}
	currents := []float64{0.0}
	voltages := []float64{0.0}

	currentMethod := map[IntegrationMethod]string{GearMethod: "Gear/BDF", TrapezoidalMethod: "Trapezoidal"}[integrationMethod]
	fmt.Printf("Integration Method: %s (Order %d)\n\n", currentMethod, methodOrder)

	fmt.Printf("\n%8s | %12s | %12s | %10s | %5s | %6s\n", "Time (s)", "iL (A)", "vL (V)", "Step (s)", "Order", "LTE")
	fmt.Println("--------------------------------------------------------------------------")

	A, err := sparse.Create(4, config)
	if err != nil {
		panic(err)
	}
	defer A.Destroy()

	// DC Operating point
	b := make([]float64, A.Size+1)
	A.Clear()

	A.GetElement(1, 1).Real += G
	A.GetElement(1, 3).Real += 1.0
	A.GetElement(2, 3).Real += 1.0
	A.GetElement(2, 4).Real += 1.0
	A.GetElement(3, 2).Real += 1.0
	A.GetElement(4, 1).Real += -1.0
	A.GetElement(4, 2).Real += 1.0

	b[3] = Vpeak * math.Sin(2.0*math.Pi*freq*t)

	A.MNAPreorder()
	A.Factor()

	x, err := A.Solve(b)
	if err != nil {
		panic(err)
	}

	currents[0] = x[3]
	voltages[0] = x[2] - x[1]

	fmt.Printf("%8.6f | %12.8f | %12.8f | %10.2e | %5d | %10.2e\n", t, currents[0], voltages[0], dt, 1, 0.0)

	// Transient
	currentOrder := 1
	for t < tstop {
		currentOrder = len(currents)
		if len(currents) > methodOrder {
			currentOrder = methodOrder
		}

		coeffs := GetIntegratorCoeffs(currentOrder, dt)

		tNext := t + dt
		if tNext > tstop {
			tNext = tstop
		}

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
			if currentOrder == 1 {
				b[4] = coeffs[0] * currents[len(currents)-1] * L // Backward Euler
			} else {
				b[4] = coeffs[0]*currents[len(currents)-1]*L - voltages[len(voltages)-1] // Trapezoidal
			}
		default:
			for j := 1; j <= currentOrder; j++ {
				b[4] -= coeffs[j] * currents[len(currents)-j]
			}
			b[4] *= L
		}

		A.MNAPreorder()
		A.Factor()

		x, err = A.Solve(b)
		if err != nil {
			panic(err)
		}

		current := x[3]
		voltage := x[2] - x[1]

		// Error control
		lte := calculateLTE(currents, dt, currentOrder)
		if lte > targetLTE && len(currents) >= currentOrder+1 {
			dt = calculateNewTimeStep(dt, lte, currentOrder)
			continue
		}

		t = tNext
		timePoints = append(timePoints, t)
		currents = append(currents, current)
		voltages = append(voltages, voltage)

		dt = calculateNewTimeStep(dt, lte, currentOrder)

		if len(timePoints)%10 == 0 || t == tstop {
			fmt.Printf("%8.6f | %12.8f | %12.8f | %10.2e | %5d | %10.2e\n", t, current, voltage, dt, currentOrder, lte)
		}
	}

	// Results compare with thorytical max voltage
	maxVL := 0.0
	for _, v := range voltages {
		if math.Abs(v) > maxVL {
			maxVL = math.Abs(v)
		}
	}

	maxDiDt := Vpeak * 2 * math.Pi * freq / math.Sqrt(R*R+math.Pow(2*math.Pi*freq*L, 2))
	theoryMaxVL := L * maxDiDt

	fmt.Println("\nResults Summary:")
	fmt.Printf("Theoretical max VL: %.6f V\n", theoryMaxVL)
	fmt.Printf("Computed max VL: %.6f V\n", maxVL)
	errPct := 100 * math.Abs(maxVL-theoryMaxVL) / theoryMaxVL
	fmt.Printf("Error: %.6f%%\n", errPct)
	fmt.Printf("Final timestep: %.8e s\n", dt)
	fmt.Printf("Final order: %d\n", currentOrder)
}
