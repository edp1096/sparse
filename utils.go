package sparse

import (
	"fmt"
	"math"
)

func (m *Matrix) ElementCount() int {
	return m.Elements
}

func (m *Matrix) FillinCount() int {
	return m.Fillins
}

func (m *Matrix) GetSize(external bool) int64 {
	if m.Config.Translate && external {
		return m.ExtSize
	}
	return m.Size
}

// Calculating determinant after LU factorization. Only use after Factor and before Clear
func (m *Matrix) Determinant() (determinant float64, exponent int, imagDeterminant *float64) {
	if m == nil || !m.Factored {
		return 0.0, 0, nil
	}

	if m.SingularRow > 0 || m.SingularCol > 0 {
		var singularZero float64 = 0.0
		if m.Complex {
			return 0.0, 0, &singularZero
		}
		return 0.0, 0, nil
	}

	exponent = 0
	if m.Complex {
		detReal, detImag := 1.0, 0.0

		for i := int64(1); i <= m.Size; i++ {
			pivotReal := m.Diags[i].Real
			pivotImag := m.Diags[i].Imag
			denominator := pivotReal*pivotReal + pivotImag*pivotImag

			tempReal := (detReal*pivotReal + detImag*pivotImag) / denominator
			tempImag := (detImag*pivotReal - detReal*pivotImag) / denominator
			detReal, detImag = tempReal, tempImag

			// Scaling
			norm := math.Max(math.Abs(detReal), math.Abs(detImag))
			if norm != 0.0 {
				for norm >= 1.0e12 {
					detReal *= 1.0e-12
					detImag *= 1.0e-12
					exponent += 12
					norm = math.Max(math.Abs(detReal), math.Abs(detImag))
				}
				for norm < 1.0e-12 {
					detReal *= 1.0e12
					detImag *= 1.0e12
					exponent -= 12
					norm = math.Max(math.Abs(detReal), math.Abs(detImag))
				}
			}
		}

		// Scale to 1.0 <= x < 10.0
		norm := math.Max(math.Abs(detReal), math.Abs(detImag))
		if norm != 0.0 {
			for norm >= 10.0 {
				detReal *= 0.1
				detImag *= 0.1
				exponent++
				norm = math.Max(math.Abs(detReal), math.Abs(detImag))
			}
			for norm < 1.0 {
				detReal *= 10.0
				detImag *= 10.0
				exponent--
				norm = math.Max(math.Abs(detReal), math.Abs(detImag))
			}
		}

		if m.NumberOfInterchangesIsOdd {
			detReal = -detReal
			detImag = -detImag
		}

		return detReal, exponent, &detImag
	} else {
		det := 1.0

		for i := int64(1); i <= m.Size; i++ {
			det /= m.Diags[i].Real

			// Scaling
			if det != 0.0 {
				for math.Abs(det) >= 1.0e12 {
					det *= 1.0e-12
					exponent += 12
				}
				for math.Abs(det) < 1.0e-12 {
					det *= 1.0e12
					exponent -= 12
				}
			}
		}

		// Scale to 1.0 <= x < 10.0
		if det != 0.0 {
			for math.Abs(det) >= 10.0 {
				det *= 0.1
				exponent++
			}
			for math.Abs(det) < 1.0 {
				det *= 10.0
				exponent--
			}
		}

		if m.NumberOfInterchangesIsOdd {
			det = -det
		}

		return det, exponent, nil
	}
}

func (m *Matrix) LargestElement() float64 {
	if m == nil {
		return 0.0
	}

	if !m.Factored {
		max := 0.0
		for i := int64(1); i <= m.Size; i++ {
			for element := m.FirstInCol[i]; element != nil; element = element.NextInCol {
				mag := 0.0
				if m.Complex {
					mag = complexInfNorm(element.Real, element.Imag)
				} else {
					mag = math.Abs(element.Real)
				}
				if mag > max {
					max = mag
				}
			}
		}

		return max
	}

	// Factored case
	if m.SingularRow > 0 || m.SingularCol > 0 {
		return 0.0
	}

	maxRow := 0.0
	maxCol := 0.0

	for i := int64(1); i <= m.Size; i++ {
		diag := m.Diags[i]
		if m.Complex {
			// Lower triangular matrix - Complex
			recipReal := diag.Real
			recipImag := diag.Imag
			m.complexReciprocal(diag)
			pivotNorm := complexInfNorm(diag.Real, diag.Imag)
			diag.Real = recipReal
			diag.Imag = recipImag
			if pivotNorm > maxRow {
				maxRow = pivotNorm
			}

			// Row elements
			for element := m.FirstInRow[i]; element != diag; element = element.NextInRow {
				mag := complexInfNorm(element.Real, element.Imag)
				if mag > maxRow {
					maxRow = mag
				}
			}

			// Upper triangular matrix
			absColSum := 1.0
			for element := m.FirstInCol[i]; element != diag; element = element.NextInCol {
				absColSum += complexInfNorm(element.Real, element.Imag)
			}
			if absColSum > maxCol {
				maxCol = absColSum
			}
		} else {
			// Lower triangular matrix - Real
			pivot := 1.0 / diag.Real
			mag := math.Abs(pivot)
			if mag > maxRow {
				maxRow = mag
			}

			// Row elements
			for element := m.FirstInRow[i]; element != diag; element = element.NextInRow {
				mag := math.Abs(element.Real)
				if mag > maxRow {
					maxRow = mag
				}
			}

			// Upper triangular matrix
			absColSum := 1.0
			for element := m.FirstInCol[i]; element != diag; element = element.NextInCol {
				absColSum += math.Abs(element.Real)
			}
			if absColSum > maxCol {
				maxCol = absColSum
			}
		}
	}

	return maxRow * maxCol
}

// Returns a bound on the magnitude of the largest element in E = A - LU
func (m *Matrix) Roundoff(rho float64) float64 {
	if m == nil || !m.Factored {
		return 0.0
	}

	// Compute Barlow's bound if not given
	if rho < 0.0 {
		rho = m.LargestElement()
	}

	// Find the maximum number of off-diagonals in L
	maxCount := int64(0)
	if m.MaxRowCountInLowerTri < 0 {
		for i := m.Size; i > 0; i-- {
			count := int64(0)
			for element := m.FirstInRow[i]; element != nil && element.Col < i; element = element.NextInRow {
				count++
			}
			if count > maxCount {
				maxCount = count
			}
		}
		m.MaxRowCountInLowerTri = maxCount
	} else {
		maxCount = m.MaxRowCountInLowerTri
	}

	const machineResolution = 2.2204460492503131e-016 // DBL_EPSILON

	// Compute error bound
	gear := 1.01 * ((float64(maxCount)+1)*m.RelThreshold + 1.0) * float64(maxCount*maxCount)
	reid := 3.01 * float64(m.Size)

	if gear < reid {
		return machineResolution * rho * gear
	}
	return machineResolution * rho * reid
}

// Calculates the infinity norm of the matrix
func (m *Matrix) Norm() float64 {
	if m == nil || m.Factored {
		return 0.0
	}

	if !m.RowsLinked {
		m.LinkRows()
	}

	max := 0.0

	if m.Complex {
		for i := m.Size; i > 0; i-- {
			absRowSum := 0.0
			for element := m.FirstInRow[i]; element != nil; element = element.NextInRow {
				absRowSum += complex1Norm(element.Real, element.Imag)
			}
			if max < absRowSum {
				max = absRowSum
			}
		}
	} else {
		for i := m.Size; i > 0; i-- {
			absRowSum := 0.0
			for element := m.FirstInRow[i]; element != nil; element = element.NextInRow {
				absRowSum += math.Abs(element.Real)
			}
			if max < absRowSum {
				max = absRowSum
			}
		}
	}

	return max
}

// Condition returns reciprocal of the condition number
func (m *Matrix) Condition(normOfMatrix float64) (float64, error) {
	if m == nil || !m.Factored {
		return 0.0, fmt.Errorf("matrix not valid or not factored")
	}
	if normOfMatrix == 0.0 {
		return 0.0, fmt.Errorf("singular")
	}

	if m.Complex {
		return m.complexCondition(normOfMatrix)
	}

	size := m.Size
	matrixSize := m.Size + 1 // 1-based indexing

	t := make([]float64, matrixSize)
	tm := make([]float64, matrixSize)

	// Part 1. Ay = e
	e := 1.0
	for i := int64(1); i <= size; i++ {
		pivot := m.Diags[i]
		var em float64
		if t[i] < 0.0 {
			em = -e
		} else {
			em = e
		}
		wm := (em + t[i]) * pivot.Real

		if math.Abs(wm) > SLACK {
			scaleFactor := 1.0 / math.Max(SLACK*SLACK, math.Abs(wm))
			for k := size; k > 0; k-- {
				t[k] *= scaleFactor
			}
			e *= scaleFactor
			em *= scaleFactor
			wm = (em + t[i]) * pivot.Real
		}

		wp := (t[i] - em) * pivot.Real
		asp := math.Abs(t[i] - em)
		asm := math.Abs(em + t[i])

		// Update T for both values of W
		for element := pivot.NextInCol; element != nil; element = element.NextInCol {
			row := element.Row
			tm[row] = t[row] - (wm * element.Real)
			t[row] -= (wp * element.Real)
			asp += math.Abs(t[row])
			asm += math.Abs(tm[row])
		}

		// Choose better value
		if asm > asp {
			t[i] = wm
			for element := pivot.NextInCol; element != nil; element = element.NextInCol {
				t[element.Row] = tm[element.Row]
			}
		} else {
			t[i] = wp
		}
	}

	// Compute 1-norm of T and scale
	asw := 0.0
	for i := size; i > 0; i-- {
		asw += math.Abs(t[i])
	}
	scaleFactor := 1.0 / (SLACK * asw)
	if scaleFactor < 0.5 {
		for i := size; i > 0; i-- {
			t[i] *= scaleFactor
		}
		e *= scaleFactor
	}

	// Backward Substitution
	for i := size; i >= 1; i-- {
		for element := m.Diags[i].NextInRow; element != nil; element = element.NextInRow {
			t[i] -= element.Real * t[element.Col]
		}
		if math.Abs(t[i]) > SLACK {
			scaleFactor := 1.0 / math.Max(SLACK*SLACK, math.Abs(t[i]))
			for k := size; k > 0; k-- {
				t[k] *= scaleFactor
			}
			e *= scaleFactor
		}
	}

	// More computation and scaling...
	asy := 0.0
	for i := size; i > 0; i-- {
		asy += math.Abs(t[i])
	}
	scaleFactor = 1.0 / (SLACK * asy)
	if scaleFactor < 0.5 {
		for i := size; i > 0; i-- {
			t[i] *= scaleFactor
		}
		asy = 1.0 / SLACK
		e *= scaleFactor
	}

	// Compute infinity-norm
	maxY := 0.0
	for i := size; i > 0; i-- {
		if maxY < math.Abs(t[i]) {
			maxY = math.Abs(t[i])
		}
	}

	// Part 2: A* z = y
	// Forward elimination
	for i := int64(1); i <= size; i++ {
		for element := m.Diags[i].NextInRow; element != nil; element = element.NextInRow {
			t[element.Col] -= t[i] * element.Real
		}
		if math.Abs(t[i]) > SLACK {
			scaleFactor := 1.0 / math.Max(SLACK*SLACK, math.Abs(t[i]))
			for k := size; k > 0; k-- {
				t[k] *= scaleFactor
			}
			asy *= scaleFactor
		}
	}

	// Compute 1-norm and scale
	asv := 0.0
	for i := size; i > 0; i-- {
		asv += math.Abs(t[i])
	}
	scaleFactor = 1.0 / (SLACK * asv)
	if scaleFactor < 0.5 {
		for i := size; i > 0; i-- {
			t[i] *= scaleFactor
		}
		asy *= scaleFactor
	}

	// Backward substitution
	for i := size; i >= 1; i-- {
		pivot := m.Diags[i]
		for element := pivot.NextInCol; element != nil; element = element.NextInCol {
			t[i] -= element.Real * t[element.Row]
		}
		t[i] *= pivot.Real
		if math.Abs(t[i]) > SLACK {
			scaleFactor := 1.0 / math.Max(SLACK*SLACK, math.Abs(t[i]))
			for k := size; k > 0; k-- {
				t[k] *= scaleFactor
			}
			asy *= scaleFactor
		}
	}

	// Compute final 1-norm
	asz := 0.0
	for i := size; i > 0; i-- {
		asz += math.Abs(t[i])
	}

	// Compute condition estimates
	linpack := asy / asz
	oleary := e / maxY
	invNormOfInverse := math.Min(linpack, oleary)

	return invNormOfInverse / normOfMatrix, nil
}

type ComplexTemp struct {
	Real, Imag float64
}

// complexCondition returns reciprocal of the condition number for complex matrices
func (m *Matrix) complexCondition(normOfMatrix float64) (float64, error) {
	size := m.Size
	matrixSize := m.Size + 1 // 1-based indexing

	t := make([]ComplexTemp, matrixSize)
	tm := make([]ComplexTemp, matrixSize)

	// Initialize T to zeros
	for i := size; i > 0; i-- {
		t[i].Real = 0.0
		t[i].Imag = 0.0
	}

	// Part 1. Ay = e
	e := 1.0
	for i := int64(1); i <= size; i++ {
		pivot := m.Diags[i]

		var em float64
		if t[i].Real < 0.0 {
			em = -e
		} else {
			em = e
		}

		// Wm = T[i] + em
		wm := &Element{Real: t[i].Real + em, Imag: t[i].Imag}
		asm := complex1Norm(wm.Real, wm.Imag)
		m.complexMultAssign(wm, pivot)

		if complex1Norm(wm.Real, wm.Imag) > SLACK {
			scaleFactor := 1.0 / math.Max(SLACK*SLACK, complex1Norm(wm.Real, wm.Imag))
			for k := size; k > 0; k-- {
				t[k].Real *= scaleFactor
				t[k].Imag *= scaleFactor
			}
			e *= scaleFactor
			em *= scaleFactor
			asm *= scaleFactor
			wm.scalarMultAssign(scaleFactor)
		}

		// Wp = T[i] - em
		wp := &Element{Real: t[i].Real - em, Imag: t[i].Imag}
		asp := complex1Norm(wp.Real, wp.Imag)
		m.complexMultAssign(wp, pivot)

		// Update T for both values of W
		for element := pivot.NextInCol; element != nil; element = element.NextInCol {
			row := element.Row

			// Tm[Row] = T[Row] - (Wm * element)
			tmElement := &Element{Real: t[row].Real, Imag: t[row].Imag}
			m.complexMultSubtAssign(tmElement, wm, element)
			tm[row] = ComplexTemp{Real: tmElement.Real, Imag: tmElement.Imag}

			// T[Row] -= Wp * element
			tElement := &Element{Real: t[row].Real, Imag: t[row].Imag}
			m.complexMultSubtAssign(tElement, wp, element)
			t[row] = ComplexTemp{Real: tElement.Real, Imag: tElement.Imag}

			asp += complex1Norm(t[row].Real, t[row].Imag)
			asm += complex1Norm(tm[row].Real, tm[row].Imag)
		}

		// Choose better value
		if asm > asp {
			t[i] = ComplexTemp{Real: wm.Real, Imag: wm.Imag}
			for element := pivot.NextInCol; element != nil; element = element.NextInCol {
				t[element.Row] = tm[element.Row]
			}
		} else {
			t[i] = ComplexTemp{Real: wp.Real, Imag: wp.Imag}
		}
	}

	// Compute 1-norm of T and scale
	asw := 0.0
	for i := size; i > 0; i-- {
		asw += complex1Norm(t[i].Real, t[i].Imag)
	}
	scaleFactor := 1.0 / (SLACK * asw)
	if scaleFactor < 0.5 {
		for i := size; i > 0; i-- {
			t[i].Real *= scaleFactor
			t[i].Imag *= scaleFactor
		}
		e *= scaleFactor
	}

	// Backward Substitution
	for i := size; i >= 1; i-- {
		for element := m.Diags[i].NextInRow; element != nil; element = element.NextInRow {
			tElement := &Element{Real: t[element.Col].Real, Imag: t[element.Col].Imag}
			result := &Element{Real: t[i].Real, Imag: t[i].Imag}
			m.complexMultSubtAssign(result, tElement, element)
			t[i] = ComplexTemp{Real: result.Real, Imag: result.Imag}
		}
		if complex1Norm(t[i].Real, t[i].Imag) > SLACK {
			scaleFactor := 1.0 / math.Max(SLACK*SLACK, complex1Norm(t[i].Real, t[i].Imag))
			for k := size; k > 0; k-- {
				t[k].Real *= scaleFactor
				t[k].Imag *= scaleFactor
			}
			e *= scaleFactor
		}
	}

	// Compute 1-norm of T and scale
	asy := 0.0
	for i := size; i > 0; i-- {
		asy += complex1Norm(t[i].Real, t[i].Imag)
	}
	scaleFactor = 1.0 / (SLACK * asy)
	if scaleFactor < 0.5 {
		for i := size; i > 0; i-- {
			t[i].Real *= scaleFactor
			t[i].Imag *= scaleFactor
		}
		asy = 1.0 / SLACK
		e *= scaleFactor
	}

	// Compute infinity-norm
	maxY := 0.0
	for i := size; i > 0; i-- {
		norm := complex1Norm(t[i].Real, t[i].Imag)
		if maxY < norm {
			maxY = norm
		}
	}

	// Forward elimination
	for i := int64(1); i <= size; i++ {
		for element := m.Diags[i].NextInRow; element != nil; element = element.NextInRow {
			tElement := &Element{Real: t[element.Col].Real, Imag: t[element.Col].Imag}
			result := &Element{Real: t[i].Real, Imag: t[i].Imag}
			m.complexMultSubtAssign(result, tElement, element)
			t[i] = ComplexTemp{Real: result.Real, Imag: result.Imag}
		}
		if complex1Norm(t[i].Real, t[i].Imag) > SLACK {
			scaleFactor := 1.0 / math.Max(SLACK*SLACK, complex1Norm(t[i].Real, t[i].Imag))
			for k := size; k > 0; k-- {
				t[k].Real *= scaleFactor
				t[k].Imag *= scaleFactor
			}
			asy *= scaleFactor
		}
	}

	// Compute 1-norm and scale
	asv := 0.0
	for i := size; i > 0; i-- {
		asv += complex1Norm(t[i].Real, t[i].Imag)
	}
	scaleFactor = 1.0 / (SLACK * asv)
	if scaleFactor < 0.5 {
		for i := size; i > 0; i-- {
			t[i].Real *= scaleFactor
			t[i].Imag *= scaleFactor
		}
		asy *= scaleFactor
	}

	// Backward substitution
	for i := size; i >= 1; i-- {
		pivot := m.Diags[i]
		for element := pivot.NextInCol; element != nil; element = element.NextInCol {
			tElement := &Element{Real: t[element.Row].Real, Imag: t[element.Row].Imag}
			result := &Element{Real: t[i].Real, Imag: t[i].Imag}
			m.complexMultSubtAssign(result, tElement, element)
			t[i] = ComplexTemp{Real: result.Real, Imag: result.Imag}
		}

		// T[I] *= pivot
		tElement := &Element{Real: t[i].Real, Imag: t[i].Imag}
		m.complexMultAssign(tElement, pivot)
		t[i] = ComplexTemp{Real: tElement.Real, Imag: tElement.Imag}

		if complex1Norm(t[i].Real, t[i].Imag) > SLACK {
			scaleFactor := 1.0 / math.Max(SLACK*SLACK, complex1Norm(t[i].Real, t[i].Imag))
			for k := size; k > 0; k-- {
				t[k].Real *= scaleFactor
				t[k].Imag *= scaleFactor
			}
			asy *= scaleFactor
		}
	}

	// final 1-norm
	asz := 0.0
	for i := size; i > 0; i-- {
		asz += complex1Norm(t[i].Real, t[i].Imag)
	}

	linpack := asy / asz
	oleary := e / maxY
	invNormOfInverse := math.Min(linpack, oleary)

	return invNormOfInverse / normOfMatrix, nil
}

func (m *Matrix) PseudoCondition() float64 {
	if m == nil || !m.Factored || m.SingularRow > 0 || m.SingularCol > 0 {
		return 0.0
	}

	var mag float64
	if m.Complex {
		mag = complexInfNorm(m.Diags[1].Real, m.Diags[1].Imag)
	} else {
		mag = math.Abs(m.Diags[1].Real)
	}
	maxPivot := mag
	minPivot := mag

	for i := int64(2); i <= m.Size; i++ {
		if m.Complex {
			mag = complexInfNorm(m.Diags[i].Real, m.Diags[i].Imag)
		} else {
			mag = math.Abs(m.Diags[i].Real)
		}
		if mag > maxPivot {
			maxPivot = mag
		} else if mag < minPivot {
			minPivot = mag
		}
	}

	if maxPivot <= 0.0 {
		return 0.0
	}
	return maxPivot / minPivot
}

func (m *Matrix) Multiply(solution []float64, isolution []float64) ([]float64, []float64, error) {
	if !m.RowsLinked {
		m.LinkRows()
	}
	if !m.InternalVectorsAllocated {
		if err := m.CreateInternalVectors(); err != nil {
			return nil, nil, fmt.Errorf("failed to create internal vectors: %v", err)
		}
	}

	top := m.Size
	if m.Config.Translate {
		top = m.ExtSize
	}

	// rhs := make([]float64, m.Size+1)
	rhs := make([]float64, top+1)
	var irhs []float64
	if m.Complex && m.Config.SeparatedComplexVectors {
		// irhs = make([]float64, m.Size+1)
		irhs = make([]float64, top+1)
	}

	if m.Complex {
		return m.MultiplyComplexMatrix(solution, isolution)
	}

	for i := int64(1); i <= m.Size; i++ {
		m.Intermediate[i] = solution[m.IntToExtColMap[i]]
	}

	for i := int64(1); i <= m.Size; i++ {
		element := m.FirstInRow[i]
		sum := 0.0

		for element != nil {
			sum += element.Real * m.Intermediate[element.Col]
			element = element.NextInRow
		}
		rhs[m.IntToExtRowMap[i]] = sum
	}

	return rhs, irhs, nil
}

func (m *Matrix) MultiplyComplexMatrix(solution []float64, isolution []float64) ([]float64, []float64, error) {
	separated := m.Config.SeparatedComplexVectors
	matrixSize := m.Size + 1 // 1-based indexing

	rhs := make([]float64, matrixSize)
	var irhs []float64
	if separated {
		irhs = make([]float64, matrixSize)
	} else {
		rhs = make([]float64, 2*(matrixSize))
	}

	vector := make([]Element, matrixSize)
	for i := int64(1); i <= m.Size; i++ {
		extIdx := m.IntToExtColMap[i]
		if separated {
			vector[i].Real = solution[extIdx]
			vector[i].Imag = isolution[extIdx]
		} else {
			vector[i].Real = solution[2*extIdx]
			vector[i].Imag = solution[2*extIdx+1]
		}
	}

	for i := int64(1); i <= m.Size; i++ {
		element := m.FirstInRow[i]
		sum := Element{}

		for element != nil {
			m.complexMultAddAssign(&sum, element, &vector[element.Col])
			element = element.NextInRow
		}

		extIdx := m.IntToExtRowMap[i]
		if separated {
			rhs[extIdx] = sum.Real
			irhs[extIdx] = sum.Imag
		} else {
			rhs[2*extIdx] = sum.Real
			rhs[2*extIdx+1] = sum.Imag
		}
	}

	return rhs, irhs, nil
}

func (m *Matrix) MultplyTransposed(solution []float64, isolution []float64) ([]float64, []float64, error) {
	if !m.InternalVectorsAllocated {
		if err := m.CreateInternalVectors(); err != nil {
			return nil, nil, fmt.Errorf("failed to create internal vectors: %v", err)
		}
	}

	matrixSize := m.Size + 1 // 1-based indexing

	// Create result vectors with proper size
	rhs := make([]float64, matrixSize)
	var irhs []float64
	if m.Complex && m.Config.SeparatedComplexVectors {
		irhs = make([]float64, matrixSize)
	}

	if m.Complex {
		return m.MultiplyComplexTransposedMatrix(solution, isolution)
	}

	// Initialize Intermediate vector with reordered Solution vector
	for i := int64(1); i <= m.Size; i++ {
		m.Intermediate[i] = solution[m.IntToExtRowMap[i]]
	}

	// Multiply transposed matrix by intermediate vector
	for i := int64(1); i <= m.Size; i++ {
		element := m.FirstInCol[i]
		sum := 0.0

		for element != nil {
			sum += element.Real * m.Intermediate[element.Row]
			element = element.NextInCol
		}
		rhs[m.IntToExtColMap[i]] = sum
	}

	return rhs, irhs, nil
}

func (m *Matrix) MultiplyComplexTransposedMatrix(solution []float64, isolution []float64) ([]float64, []float64, error) {
	separated := m.Config.SeparatedComplexVectors
	matrixSize := m.Size + 1 // 1-based indexing

	// Create result vectors
	rhs := make([]float64, matrixSize)
	var irhs []float64
	if separated {
		irhs = make([]float64, matrixSize)
	} else {
		rhs = make([]float64, 2*(matrixSize))
	}

	vector := make([]Element, matrixSize)
	for i := int64(1); i <= m.Size; i++ {
		extIdx := m.IntToExtRowMap[i]
		if separated {
			vector[i].Real = solution[extIdx]
			vector[i].Imag = isolution[extIdx]
		} else {
			vector[i].Real = solution[2*extIdx-1]
			vector[i].Imag = solution[2*extIdx]
		}
	}

	for i := int64(1); i <= m.Size; i++ {
		element := m.FirstInCol[i]
		sum := Element{}

		for element != nil {
			m.complexMultAddAssign(&sum, element, &vector[element.Row])
			element = element.NextInCol
		}

		extIdx := m.IntToExtColMap[i]
		if separated {
			rhs[extIdx] = sum.Real
			irhs[extIdx] = sum.Imag
		} else {
			rhs[2*extIdx-1] = sum.Real
			rhs[2*extIdx] = sum.Imag
		}
	}

	return rhs, irhs, nil
}

func (m *Matrix) CalculateNormalizedResidual(rhs, solution, irhs, isolution []float64) (float64, float64, error) {
	size := m.Size
	isComplex := m.Complex

	m.Initialize()

	maxRHS := 0.0
	if isComplex {
		if m.Config.SeparatedComplexVectors {
			for i := int64(1); i <= size; i++ {
				maxRHS = max(maxRHS, math.Abs(rhs[i]))
				maxRHS = max(maxRHS, math.Abs(irhs[i]))
			}
		} else {
			for i := int64(1); i <= size; i++ {
				idx := 2*i - 1
				maxRHS = max(maxRHS, math.Abs(rhs[idx]))
				maxRHS = max(maxRHS, math.Abs(rhs[idx+1]))
			}
		}
	} else {
		for i := int64(1); i <= size; i++ {
			maxRHS = max(maxRHS, math.Abs(rhs[i]))
		}
	}

	var rhsVerif, irhsVerif []float64
	var err error
	if m.Config.Transpose {
		rhsVerif, irhsVerif, err = m.MultplyTransposed(solution, isolution)
	} else {
		rhsVerif, irhsVerif, err = m.Multiply(solution, isolution)
	}
	if err != nil {
		return 0, 0, fmt.Errorf("failed to multiply matrix: %v", err)
	}

	residual := 0.0
	if isComplex {
		if m.Config.SeparatedComplexVectors {
			for i := int64(1); i <= size; i++ {
				residual += math.Abs(rhs[i]-rhsVerif[i]) + math.Abs(irhs[i]-irhsVerif[i])
			}
		} else {
			for i := int64(1); i <= size; i++ {
				idx := 2*i - 1
				residual += math.Abs(rhs[idx]-rhsVerif[idx]) +
					math.Abs(rhs[idx+1]-rhsVerif[idx+1])
			}
		}
	} else {
		for i := int64(1); i <= size; i++ {
			residual += math.Abs(rhs[i] - rhsVerif[i])
		}
	}

	if maxRHS == 0.0 {
		return 0, 0, nil
	}
	return residual / maxRHS, maxRHS, nil
}
