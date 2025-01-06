package main

import (
	"fmt"
	"math"

	"github.com/edp1096/sparse"
)

type Complex struct {
	Re float64
	Im float64
}

func main() {
	var err error

	annotate := 0
	separatedComplexVectors := false
	defaultPartition := sparse.AUTO_PARTITION

	config := &sparse.Configuration{
		Real:                    true,
		Complex:                 true,
		SeparatedComplexVectors: separatedComplexVectors,
		Expandable:              true,
		Translate:               true,
		ModifiedNodal:           true,
		TiesMultiplier:          5,
		DefaultPartition:        defaultPartition,
		PrinterWidth:            140,
		Annotate:                annotate,
	}

	A, err := sparse.Create(2, config)
	if err != nil {
		panic(err)
	}

	stamps := make([]sparse.Template, 3)
	A.GetAdmittance(1, 0, &stamps[0])
	A.GetAdmittance(1, 2, &stamps[1])
	A.GetAdmittance(2, 0, &stamps[2])

	/* Drive the circuit at node 1. */
	// SeparatedComplexVectors = true
	b := make([]float64, 3)
	ib := make([]float64, 3)
	b[1] = 1.0
	ib[1] = 0.0
	b[2] = 0.0
	ib[2] = 0.0

	// SeparatedComplexVectors = false
	sb := make([]float64, 6)
	sb[2] = 1.0 // node1 real
	sb[3] = 0.0 // node1 imag
	sb[4] = 0.0 // node2 real
	sb[5] = 0.0 // node2 imag

	/* Perform AC analysis over a range of frequencies. */
	// for f := 0.0; f <= 100000.0; f += 1000.0 {
	for f := 0.0; f <= 2000.0; f += 100.0 {
		omega := 2.0 * math.Pi * f

		A.Clear()
		stamps[0].AddComplexQuad(1.0/50.0, 10e-6*omega)
		stamps[1].AddRealQuad(1.0 / 200.0)
		stamps[2].AddComplexQuad(1.0/50.0, 10e-6*omega)

		err = A.Factor()
		if err != nil {
			panic(err)
		}

		var magnitude float64
		if config.SeparatedComplexVectors {
			xReal, xImag, err := A.SolveComplex(b, ib)
			if err != nil {
				panic(err)
			}
			magnitude = math.Sqrt(xReal[2]*xReal[2] + xImag[2]*xImag[2])
		} else {
			x, _, err := A.SolveComplex(sb, nil)
			if err != nil {
				panic(err)
			}
			magnitude = math.Sqrt(x[4]*x[4] + x[5]*x[5])
		}

		db := 20 * math.Log10(magnitude)

		fmt.Printf("f = %04.0f Hz, h = %.6f, Gain = %.2f dB\n", f, magnitude, db)
	}

	A.Destroy()
}
