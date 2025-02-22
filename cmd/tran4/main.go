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
	R     = 100.0   // Resistor: 100 ohm
	G     = 1.0 / R // Conductance
	L     = 1e-3    // Inductor: 1mH
	Vpeak = 5.0     // Peak voltage: 5V
	freq  = 1000.0  // Frequency: 1kHz

	targetLTE    = 1e-6
	minTimeStep  = 1e-9
	maxStepSize  = 5e-5
	safetyFactor = 0.9

	StartTime = 0.0
	EndTime   = 0.002
	TimeStep  = 1e-5
)

var (
	integrationMethod = GearMethod // TrapezoidalMethod or GearMethod
	methodOrder       = 6          // Trapezoidal: 1 or 2, Gear: 1 ~ 6

	BdfCoefficients = [6]BackwardDifferentialFormula{
		{[]float64{1.0}, 1.0},                                                                                                // 1st order
		{[]float64{4.0 / 3.0, -1.0 / 3.0}, 2.0 / 3.0},                                                                        // 2nd order
		{[]float64{18.0 / 11.0, -9.0 / 11.0, 2.0 / 11.0}, 6.0 / 11.0},                                                        // 3rd order
		{[]float64{48.0 / 25.0, -36.0 / 25.0, 16.0 / 25.0, -3.0 / 25.0}, 12.0 / 25.0},                                        // 4th order
		{[]float64{300.0 / 137.0, -300.0 / 137.0, 200.0 / 137.0, -75.0 / 137.0, 12.0 / 137.0}, 60.0 / 137.0},                 // 5th order
		{[]float64{360.0 / 147.0, -450.0 / 147.0, 400.0 / 147.0, -225.0 / 147.0, 72.0 / 147.0, -10.0 / 147.0}, 60.0 / 147.0}, // 6th order
	}

	LteTrapCoeffs = []float64{0.5, 1.0 / 12.0}
	LteBdfCoeffs  = []float64{0.5, 2.0 / 3.0, 6.0 / 11.0, 12.0 / 25.0, 60.0 / 137.0, 60.0 / 147.0}

	StabilityFactors = []float64{2.00, 2.00, 1.98, 1.92, 1.76, 1.56} // stability factor
)

func GetBDFcoeffs(order int, dt float64) []float64 {
	if order < 1 || order > 6 {
		panic(fmt.Sprintf("Invalid BDF order: %d", order))
	}

	bdf := BdfCoefficients[order-1]
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
		lteCoeff = LteTrapCoeffs[order-1]
	default:
		lteCoeff = LteBdfCoeffs[order-1]
	}

	return math.Abs(lteCoeff * math.Pow(dt, float64(order+1)) * der[0])
}

func calculateNewTimeStep(currentStep, lte float64, order int, timestep float64) float64 {
	if lte < 1e-15 {
		return currentStep
	}

	factor := math.Pow(targetLTE/lte, 1.0/float64(order+1)) * safetyFactor
	factor = math.Max(0.1, math.Min(factor, 10.0))

	newStep := currentStep * factor
	maxStep := StabilityFactors[order-1] * timestep

	newStep = math.Min(newStep, maxStep)
	newStep = math.Max(newStep, minTimeStep)

	return newStep
}

func main() {
	if integrationMethod == TrapezoidalMethod && methodOrder > 2 {
		fmt.Println("Trapezoidal method order must be 1 or 2. Run as 2")
		methodOrder = 2
	}
	if integrationMethod == GearMethod && (methodOrder < 1 || methodOrder > 6) {
		fmt.Println("BDF/Gear method order must be 1 ~ 6. Run as 6")
		methodOrder = 6
	}

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

	t := StartTime
	endTime := EndTime
	timestep := TimeStep
	dt := timestep

	currents := []float64{0.0}
	voltages := []float64{0.0}
	vL := make([]float64, 0)

	b := make([]float64, A.Size+1)
	lte := 0.0

	currentMethod := map[IntegrationMethod]string{GearMethod: "Gear/BDF", TrapezoidalMethod: "Trapezoidal"}[integrationMethod]
	fmt.Printf("Integration Method: %s (Order %d)\n\n", currentMethod, methodOrder)

	fmt.Printf("Time (s) | Vin | V(1) | V(2) | I_source | I_L | dt | Order | LTE\n")
	fmt.Println("-------------------------------------------------------------------------")

	for t <= endTime {
		currentOrder := min(len(currents), methodOrder)
		coeffs := GetIntegratorCoeffs(currentOrder, dt)

		A.Clear()

		A.GetElement(1, 1).Real += G
		A.GetElement(1, 2).Real += -G
		A.GetElement(1, 3).Real += 1.0

		A.GetElement(2, 1).Real += -G
		A.GetElement(2, 2).Real += G
		A.GetElement(2, 4).Real += 1.0

		A.GetElement(3, 1).Real += 1.0

		A.GetElement(4, 2).Real += 1.0
		A.GetElement(4, 4).Real = -coeffs[0] * L

		vin := Vpeak * math.Sin(2.0*math.Pi*freq*t)

		b[1] = 0.0
		b[2] = 0.0
		b[3] = vin
		b[4] = 0.0

		switch integrationMethod {
		case TrapezoidalMethod:
			if currentOrder == 1 {
				b[4] = -coeffs[0] * L * currents[len(currents)-1]
			} else {
				b[4] = -coeffs[0]*L*(currents[len(currents)-1]) - voltages[len(voltages)-1]
			}

		default:
			maxOrder := math.Min(float64(currentOrder), float64(len(currents)-1))
			for j := 1; j <= int(maxOrder); j++ {
				b[4] += coeffs[j] * currents[len(currents)-j]
			}
			b[4] *= L
		}

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

		if currentOrder > 1 {
			lte = calculateLTE(currents, dt, currentOrder)
			if lte > targetLTE && len(currents) >= currentOrder+1 {
				dt = calculateNewTimeStep(dt, lte, currentOrder, timestep)
				continue
			}
		}

		currents = append(currents, x[4])
		voltages = append(voltages, x[2])
		vL = append(vL, x[2])

		fmt.Printf("%.6f | %.5f | %.5f | %.5f | %.5f | %.5f | %.2e | %d | %.2e\n", t, vin, x[1], x[2], x[3], x[4], dt, currentOrder, lte)

		tPrev := t
		dtPrev := dt
		t += dt
		dt = calculateNewTimeStep(dt, lte, currentOrder, timestep)
		if t > endTime && tPrev < endTime {
			t = endTime
			dt = dtPrev
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

	fmt.Println("\nCircuit Parameters:")
	fmt.Printf("R: %.0f Ohm\n", R)
	fmt.Printf("L: %.0f mH\n", L*1e3)
	fmt.Printf("Vpeak: %.1f V\n", Vpeak)
	fmt.Printf("Frequency: %.0f Hz\n", freq)
	fmt.Printf("Time constant (L/R): %.3f ms\n", (L/R)*1000)
	fmt.Printf("Theory max VL: %.6f V\n", theoryMaxVL)
	fmt.Printf("Max VL: %.6f V\n", maxVL)
	fmt.Printf("Err: %.6f%%\n", 100*math.Abs(maxVL-theoryMaxVL)/theoryMaxVL)
}
