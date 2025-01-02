package sparse

import (
	"math"

	"golang.org/x/exp/constraints"
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

func (m *Matrix) elementMag(e *Element) float64 {
	if m.Complex {
		return math.Abs(e.Real) + math.Abs(e.Imag)
	}
	return math.Abs(e.Real)
}

// complexMultAssign sets a = a * b for complex numbers
func (m *Matrix) complexMultAssign(a, b *Element) {
	aReal := a.Real
	a.Real = aReal*b.Real - a.Imag*b.Imag
	a.Imag = aReal*b.Imag + a.Imag*b.Real
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

func min[T constraints.Ordered](a, b T) T {
	if a < b {
		return a
	}
	return b
}
