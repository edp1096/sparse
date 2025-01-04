package sparse

const (
	SLACK float64 = 1e4

	DEFAULT_PARTITION  int = 0
	DIRECT_PARTITION   int = 1
	INDIRECT_PARTITION int = 2
	AUTO_PARTITION     int = 3
)

// Replace spConfig
type Configuration struct {
	Real                    bool
	Complex                 bool
	SeparatedComplexVectors bool

	Expandable        bool
	Translate         bool
	Initialize        bool
	DiagonalPivoting  bool
	ArrayOffset       bool // Not use. array offset as 0 flag for fortran. default offset is 1
	ModifiedMarkowitz bool // Not use
	Delete            bool // Not use
	Strip             bool // Not use
	ModifiedNodal     bool
	QuadElement       bool // Not use, regardless getAdmittance, getQuad using
	Transpose         bool // Flag for transpose job
	Scaling           bool // Not use
	Documentation     bool // Not use. fortran
	Stability         bool
	Condition         bool
	PseudoCondition   bool
	Determinant       bool
	Multiplication    bool
	Fortran           bool // Not use. fortran
	Debug             bool // Not use

	DefaultThreshold      float64 // For relative threshold
	DiagPivotingAsDefault bool
	SpaceForElements      int
	SpaceForFillIns       int
	ElementsPerAllocation int
	MinimumAllocatedSize  int
	ExpansionFactor       float64
	MaxMarkowitzTies      int64
	TiesMultiplier        int
	DefaultPartition      int
	PrinterWidth          int // Default: 80
	Annotate              int // 0: None, 1: OnStrangeBehavior , 2: Full
}

type Matrix struct {
	Config Configuration

	Size        int64 // Matrix size
	ExtSize     int64 // Matrix external size
	CurrentSize int64 // Matrix size after translate
	Complex     bool  // Matrix is complex

	DoRealDirect    []bool // Address mode real - Partition function
	DoComplexDirect []bool // Address mode complex - Partition function
	OperationCount  int    // Operation count for inner loop of factorization - Partition function

	// Elements      *list.List // Element list. not use

	Diags                 []*Element // Diagonal elements (as reciprocal) [1...Size]
	FirstInRow            []*Element // First element in each row [1...Size]
	FirstInCol            []*Element // First element in each column [1...Size]
	Intermediate          []float64  // Temporary vector for rhs, solution, etc. [1...Size]
	MarkowitzRow          []int64    // Markowitz counts of each row [1...Size]
	MarkowitzCol          []int64    // Markowitz counts of each column [1...Size]
	MarkowitzProd         []int64    // Markowitz products [1...Size]
	MaxRowCountInLowerTri int64      // Maximum number of off-diagonals in L
	RelThreshold          float64    // Relative threshold
	AbsThreshold          float64    // Absolute threshold

	// Factoring status flags
	NeedsOrdering             bool // reorder neccessary
	NumberOfInterchangesIsOdd bool // when reorder, number of interchanges is odd
	Partitioned               bool // partitioned
	Factored                  bool // factor done
	Reordered                 bool // reorder done
	RowsLinked                bool // rows linked

	SingularRow int64 // Singular row number
	SingularCol int64 // Singular column number

	// Counts
	Elements   int // Element count
	Fillins    int // Fill-in count
	Singletons int // singleton count

	// Pivot
	PivotsOriginalRow    int64 // Original pivot row number
	PivotsOriginalCol    int64 // Original pivot column number
	PivotSelectionMethod byte  // pivot choose method ('s', 'q', 'd', 'e')

	// Internal vectors allocation state - not use
	InternalVectorsAllocated bool

	IntToExtRowMap []int64 // Internal->External rows map [1...Size]
	IntToExtColMap []int64 // Internal->External columns map [1...Size]
	ExtToIntRowMap []int64 // External->Internal rows map [1...Size]
	ExtToIntColMap []int64 // External->Internal columns map [1...Size]

	// TrashCan *Element // Not use yet. eat more memory
}

type ComplexNumber struct {
	Real float64
	Imag float64
}

type Element struct {
	Real      float64
	Imag      float64
	Row       int64
	Col       int64
	NextInRow *Element
	NextInCol *Element
	InitInfo  *ComplexNumber
}

type Template struct {
	Element1        *Element
	Element2        *Element
	Element3Negated *Element
	Element4Negated *Element
}
