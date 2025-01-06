package main

import (
	"github.com/edp1096/sparse"
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

	// stamps := make([]*sparse.Element, 11)
	// stamps[0] = A.GetElement(1, 1)
	// stamps[1] = A.GetElement(1, 4)
	// stamps[2] = A.GetElement(2, 2)
	// stamps[3] = A.GetElement(2, 3)
	// stamps[4] = A.GetElement(3, 2)
	// stamps[5] = A.GetElement(3, 3)
	// stamps[6] = A.GetElement(4, 1)
	// stamps[7] = A.GetElement(4, 4)
	// stamps[8] = A.GetElement(4, 5)
	// stamps[9] = A.GetElement(5, 4)
	// stamps[10] = A.GetElement(5, 5)

	// stamps[0].Real += 10
	// stamps[1].Real += 4
	// stamps[2].Real += 20
	// stamps[3].Real += 5
	// stamps[4].Real += 2
	// stamps[5].Real += 30
	// stamps[6].Real += 4
	// stamps[7].Real += 40
	// stamps[8].Real += 6
	// stamps[9].Real += 6
	// stamps[10].Real += 50

	A.GetElement(1, 1).Real += 10
	A.GetElement(1, 4).Real += 4
	A.GetElement(2, 2).Real += 20
	A.GetElement(2, 3).Real += 5
	A.GetElement(3, 2).Real += 2
	A.GetElement(3, 3).Real += 30
	A.GetElement(4, 1).Real += 4
	A.GetElement(4, 4).Real += 40
	A.GetElement(4, 5).Real += 6
	A.GetElement(5, 4).Real += 6
	A.GetElement(5, 5).Real += 50

	A.Print(false, true, true)

	A.Factor()

	A.Print(false, true, true)
}
