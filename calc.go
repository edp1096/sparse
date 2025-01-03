package sparse

import "math"

func (m *Matrix) elementMag(e *Element) float64 {
	if m.Complex {
		return math.Abs(e.Real) + math.Abs(e.Imag)
	}
	return math.Abs(e.Real)
}

/* Macro function that multiply-assigns a complex number by a scalar. */
func (e *Element) scalarMultAssign(scalar float64) {
	e.Real *= scalar
	e.Imag *= scalar
}

// complex1Norm returns 1-norm of a complex number (|real| + |imag|)
func complex1Norm(real, imag float64) float64 {
	return math.Abs(real) + math.Abs(imag)
}

// complexInfNorm returns the infinity norm of a complex number
func complexInfNorm(real, imag float64) float64 {
	return math.Max(math.Abs(real), math.Abs(imag))
}

// complexMultAssign sets a = a * b for complex numbers
func (m *Matrix) complexMultAssign(a, b *Element) {
	aReal := a.Real
	a.Real = aReal*b.Real - a.Imag*b.Imag
	a.Imag = aReal*b.Imag + a.Imag*b.Real
}

// Use from nothing
// // complexMultAdd sets to = (mult_a * mult_b) + add for complex numbers
// func (m *Matrix) complexMultAdd(to, mult_a, mult_b, add *Element) {
// 	to.Real = mult_a.Real*mult_b.Real - mult_a.Imag*mult_b.Imag + add.Real
// 	to.Imag = mult_a.Real*mult_b.Imag + mult_a.Imag*mult_b.Real + add.Imag
// }

// complexMultAddAssign sets to += (from_a * from_b) for complex numbers
func (m *Matrix) complexMultAddAssign(to, from_a, from_b *Element) {
	to.Real += from_a.Real*from_b.Real - from_a.Imag*from_b.Imag
	to.Imag += from_a.Real*from_b.Imag + from_a.Imag*from_b.Real
}

// complexMultSubtAssign sets result = result - (a * b) for complex numbers
func (m *Matrix) complexMultSubtAssign(result, a, b *Element) {
	result.Real -= a.Real*b.Real - a.Imag*b.Imag
	result.Imag -= a.Real*b.Imag + a.Imag*b.Real
}

// complexReciprocal sets e = 1 / e
func (m *Matrix) complexReciprocal(e *Element) {
	if (e.Real >= e.Imag && e.Real > -e.Imag) || (e.Real < e.Imag && e.Real <= -e.Imag) {
		r := e.Imag / e.Real
		e.Real = 1.0 / (e.Real + r*e.Imag)
		e.Imag = -r * e.Real
	} else {
		r := e.Real / e.Imag
		e.Imag = -1.0 / (e.Imag + r*e.Real)
		e.Real = -r * e.Imag
	}
}

/* Element, Quad Template */

func (e *Element) AddComplexElement(real, imag float64) {
	e.Real += real
	e.Imag += imag
}

func (t *Template) AddRealQuad(real float64) {
	t.Element1.Real += real
	t.Element2.Real += real
	t.Element3Negated.Real -= real
	t.Element4Negated.Real -= real
}

func (t *Template) AddImagQuad(imag float64) {
	t.Element1.Imag += imag
	t.Element2.Imag += imag
	t.Element3Negated.Imag -= imag
	t.Element4Negated.Imag -= imag
}

func (t *Template) AddComplexQuad(real, imag float64) {
	t.AddRealQuad(real)
	t.AddImagQuad(imag)
}
