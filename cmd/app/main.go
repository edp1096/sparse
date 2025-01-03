package main

import (
	"bufio"
	"flag"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"runtime"
	"strconv"
	"strings"
	"time"

	"sparse"
)

var (
	annotate                = 0
	separatedComplexVectors = false
	stability               = true
	condition               = true
	pseudoCondition         = true
	determinant             = true
	multiplication          = true
	defaultPartition        = sparse.AUTO_PARTITION
)

type App struct {
	matrix         *sparse.Matrix
	filename       string
	description    string
	solutionOnly   bool
	printLimit     int
	relThreshold   float64
	absThreshold   float64
	iterations     int
	useColumnAsRHS bool
	columnAsRHS    int64
	rhs            []float64
	irhs           []float64
	solution       []float64
	isolution      []float64

	largestBefore   float64
	largestAfter    float64
	roundoff        float64
	infNorm         float64
	conditionNumber float64
	psudoCondition  float64
	determinant     float64
	iDeterminant    *float64
	detExponent     int
	detTime         float64

	buildTime  float64
	factorTime float64
	solveTime  float64
	startTime  time.Time
}

func InitApp() *App {
	return &App{
		absThreshold: 0.0,
		relThreshold: 0.0,
		startTime:    time.Now(),
	}
}

func (a *App) readMatrixFromFile(filename string) error {
	file, err := os.Open(filename)
	if err != nil {
		return fmt.Errorf("error opening file: %v", err)
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)
	lineNumber := 1

	if !scanner.Scan() {
		return fmt.Errorf("empty file")
	}
	line := strings.TrimSpace(scanner.Text())

	if strings.HasPrefix(line, "Starting") {
		if strings.HasPrefix(line, "Starting complex") {
			a.matrix.Complex = true
		}
		if !scanner.Scan() {
			return fmt.Errorf("missing description")
		}
	}
	a.description = strings.TrimSpace(scanner.Text())

	// Print description at 1st line of matrix
	if !a.solutionOnly {
		fmt.Printf("\n%s\n\n", a.description)
	}
	lineNumber++

	if !scanner.Scan() {
		return fmt.Errorf("missing size information")
	}
	line = scanner.Text()
	if len(strings.TrimSpace(line)) == 0 {
		return fmt.Errorf("syntax error in file '%s' at line %d", a.filename, lineNumber)
	}

	fields := strings.Fields(line)
	if len(fields) < 1 {
		return fmt.Errorf("invalid size line")
	}

	size, err := strconv.ParseInt(fields[0], 10, 64)
	if err != nil {
		return fmt.Errorf("invalid size value: %v", err)
	}

	isComplex := false
	if len(fields) > 1 && strings.ToLower(fields[1]) == "complex" {
		isComplex = true
	}

	config := &sparse.Configuration{
		Real:                    true,
		Complex:                 isComplex,
		SeparatedComplexVectors: separatedComplexVectors,
		Expandable:              true,
		Translate:               true,
		Initialize:              true,
		ModifiedNodal:           true,
		Stability:               stability,
		Condition:               condition,
		PseudoCondition:         pseudoCondition,
		Determinant:             determinant,
		Multiplication:          multiplication,
		DefaultPartition:        defaultPartition,
		TiesMultiplier:          5,
		PrinterWidth:            140,
		Annotate:                annotate,
	}

	a.matrix, err = sparse.Create(size, config)
	if err != nil {
		return fmt.Errorf("failed to create matrix: %v", err)
	}

	a.matrix.Initialize()

	a.rhs = make([]float64, size+1)
	if isComplex {
		if config.SeparatedComplexVectors {
			a.irhs = make([]float64, size+1)
		} else {
			a.rhs = make([]float64, 2*(size+1))
		}
	}

	// RHS_Col
	rhsCol := int64(1)
	if a.useColumnAsRHS {
		rhsCol = min(a.matrix.Size, a.columnAsRHS)
	}

	matrixEnded := false
	rhsValues := make([]string, 0)

	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" {
			continue
		}

		fields := strings.Fields(line)
		if !matrixEnded {
			if len(fields) >= 3 {
				row, err := strconv.ParseInt(fields[0], 10, 64)
				if err != nil {
					continue
				}
				col, err := strconv.ParseInt(fields[1], 10, 64)
				if err != nil {
					continue
				}

				// 0 0 0 check
				if row == 0 && col == 0 {
					matrixEnded = true
					continue
				}

				real := float64(0.0)
				imag := float64(0.0)

				real, err = strconv.ParseFloat(fields[2], 64)
				if err != nil {
					continue
				}

				element := a.matrix.GetElement(row, col)
				if element != nil {
					element.Real += real
					if isComplex && len(fields) >= 4 {
						if imag, err = strconv.ParseFloat(fields[3], 64); err == nil {
							element.Imag += imag
						}
					}

					if element.InitInfo == nil {
						element.InitInfo = &sparse.ComplexNumber{}
						element.InitInfo.Real = real
						element.InitInfo.Imag = imag
					}
				}

				if col == rhsCol {
					if a.matrix.Complex {
						if a.matrix.Config.SeparatedComplexVectors {
							a.rhs[row] = real
							a.irhs[row] = imag
						} else {
							idx := 2 * row
							a.rhs[idx] = real
							a.rhs[idx+1] = imag
						}
					} else {
						a.rhs[row] = real
					}
				}
			}
		} else {
			rhsValues = append(rhsValues, line)
		}
	}

	// RHS vector
	if !a.useColumnAsRHS {
		if len(rhsValues) > 0 && strings.HasPrefix(rhsValues[0], "Beginning") {
			rhsValues = rhsValues[1:]
		}

		for i := int64(0); i < size && i < int64(len(rhsValues)); i++ {
			fields := strings.Fields(rhsValues[i])
			if len(fields) > 0 {
				if value, err := strconv.ParseFloat(fields[0], 64); err == nil {
					if isComplex && !config.SeparatedComplexVectors {
						a.rhs[2*i+2] = value
						if len(fields) > 1 {
							if imag, err := strconv.ParseFloat(fields[1], 64); err == nil {
								a.rhs[2*i+2+1] = imag
							}
						}
					} else {
						a.rhs[i+1] = value
						if isComplex && len(fields) > 1 {
							if imag, err := strconv.ParseFloat(fields[1], 64); err == nil {
								a.irhs[i+1] = imag
							}
						}
					}
				}
			}
		}
	}

	if err := scanner.Err(); err != nil {
		return fmt.Errorf("error reading file: %v", err)
	}

	if !a.solutionOnly {
		fmt.Printf("Matrix is %d x %d ", size, size)
		if isComplex {
			fmt.Printf("and complex.\n")
		} else {
			fmt.Printf("and real.\n")
		}
	}

	return nil
}

func (a *App) solve() error {
	if !a.solutionOnly && a.matrix.Config.Stability {
		a.largestBefore = a.matrix.LargestElement()
	}
	if !a.solutionOnly && a.matrix.Config.Condition {
		a.infNorm = a.matrix.Norm()
	}

	initialFactorStart := time.Now()
	if err := a.matrix.OrderAndFactor(a.rhs, a.relThreshold, a.absThreshold, true); err != nil {
		return fmt.Errorf("initial order and factor failed: %v", err)
	}

	initialFactorTime := time.Since(initialFactorStart).Seconds()

	var err error
	solveStart := time.Now()
	if a.matrix.Complex {
		if a.matrix.Config.Transpose {
			a.solution, a.isolution, err = a.matrix.SolveComplexTransposed(a.rhs, a.irhs)
		} else {
			a.solution, a.isolution, err = a.matrix.SolveComplex(a.rhs, a.irhs)
		}
		if err != nil {
			return fmt.Errorf("initial complex solve failed: %v", err)
		}
	} else {
		if a.matrix.Config.Transpose {
			a.solution, err = a.matrix.SolveTransposed(a.rhs)
		} else {
			a.solution, err = a.matrix.Solve(a.rhs)
		}
		if err != nil {
			return fmt.Errorf("initial solve failed: %v", err)
		}
	}
	a.solveTime += time.Since(solveStart).Seconds()

	if !a.solutionOnly && a.matrix.Config.Stability {
		a.largestAfter = a.matrix.LargestElement()
		a.roundoff = a.matrix.Roundoff(a.largestAfter)
	}

	var conditionTime float64
	if !a.solutionOnly && a.matrix.Config.Condition {
		conditionStart := time.Now()
		a.conditionNumber, err = a.matrix.Condition(a.infNorm)
		if err != nil {
			fmt.Printf("initial condition number failed: %v", err)
		}
		conditionTime = time.Since(conditionStart).Seconds()
	}

	if !a.solutionOnly && a.matrix.Config.PseudoCondition {
		a.psudoCondition = a.matrix.PseudoCondition()
	}

	var partitionTime float64
	if !a.solutionOnly {
		partitionStart := time.Now()
		a.matrix.Partition(sparse.DEFAULT_PARTITION)
		partitionTime = time.Since(partitionStart).Seconds()
	}

	for i := 1; i <= a.iterations; i++ {
		buildStart := time.Now()
		if err = a.matrix.Initialize(); err != nil {
			return fmt.Errorf("initialize failed: %v", err)
		}
		a.buildTime += time.Since(buildStart).Seconds()

		factorStart := time.Now()
		if err = a.matrix.Factor(); err != nil {
			return fmt.Errorf("factor failed: %v", err)
		}
		a.factorTime += time.Since(factorStart).Seconds()

		solveStart := time.Now()
		if a.matrix.Complex {
			if a.matrix.Config.Transpose {
				a.solution, a.isolution, err = a.matrix.SolveComplexTransposed(a.rhs, a.irhs)
			} else {
				a.solution, a.isolution, err = a.matrix.SolveComplex(a.rhs, a.irhs)
			}
		} else {
			if a.matrix.Config.Transpose {
				a.solution, err = a.matrix.SolveTransposed(a.rhs)
			} else {
				a.solution, err = a.matrix.Solve(a.rhs)
			}
		}
		a.solveTime += time.Since(solveStart).Seconds()
		if err != nil {
			return fmt.Errorf("solve failed: %v", err)
		}
	}

	limit := a.matrix.Size
	if !a.solutionOnly && a.printLimit > 0 && int64(a.printLimit) < limit {
		limit = int64(a.printLimit)
	}
	if !a.matrix.Complex {
		if !a.solutionOnly {
			fmt.Println("Solution:")
		}
		for i := int64(1); i <= limit; i++ {
			fmt.Printf("%-16.9g\n", a.solution[i])
		}
	} else {
		if !a.solutionOnly {
			fmt.Println("Complex solution:")
		}
		if a.matrix.Config.SeparatedComplexVectors {
			for i := int64(1); i <= limit; i++ {
				fmt.Printf("%-16.9g   %-.9g j\n", a.solution[i], a.isolution[i])
			}
		} else {
			for i := int64(1); i <= limit; i++ {
				fmt.Printf("%-16.9g   %-.9g j\n", a.solution[i*2], a.solution[i*2+1])
			}
		}
	}
	if !a.solutionOnly && limit < a.matrix.Size && limit != 0 {
		fmt.Printf("Solution list truncated.\n")
	}
	fmt.Println()

	var det float64
	if !a.solutionOnly && a.matrix.Config.Determinant {
		startTime := time.Now()
		a.determinant, a.detExponent, a.iDeterminant = a.matrix.Determinant()
		a.detTime = time.Since(startTime).Seconds()

		if a.matrix.Complex {
			det = math.Hypot(a.determinant, *a.iDeterminant)
			for det >= 10.0 {
				det *= 0.1
				a.detExponent++
			}
		} else {
			det = a.determinant
		}
	}

	var normalizedResidual, maxRHS float64
	if a.matrix.Config.Multiplication {
		normalizedResidual, maxRHS, err = a.matrix.CalculateNormalizedResidual(a.rhs, a.solution, a.irhs, a.isolution)
		if err != nil {
			fmt.Printf("initial normalized residual failed: %v", err)
		}
	}

	if !a.solutionOnly {
		additionalLines := ""

		fmt.Printf("Statistics:\n")
		fmt.Printf("Initial factor time = %.2f.\n", initialFactorTime)
		fmt.Printf("Partition time = %.2f.\n", partitionTime)
		if a.iterations > 0 {
			fmt.Printf("Build time = %.3f.\n", a.buildTime/float64(a.iterations))
			fmt.Printf("Factor time = %.3f.\n", a.factorTime/float64(a.iterations))
			fmt.Printf("Solve time = %.3f.\n", a.solveTime/float64(a.iterations))
			if a.matrix.Config.Condition {
				fmt.Printf("Condition time = %.3f.\n", conditionTime/float64(a.iterations))
			}
		}

		if a.matrix.Config.Stability {
			if a.largestBefore != 0.0 {
				additionalLines += fmt.Sprintf("Growth = %.2g\n", a.largestAfter/a.largestBefore)
			}
			additionalLines += fmt.Sprintf("Max error in matrix = %.2g\n", a.roundoff)
		}
		if a.matrix.Config.Condition {
			additionalLines += fmt.Sprintf("Condition number = %.2g\n", 1/a.conditionNumber)
		}
		if a.matrix.Config.Condition && a.matrix.Config.Stability {
			additionalLines += fmt.Sprintf("Estimated upper bound of error in solution = %.2g\n", a.roundoff/a.conditionNumber)
		}
		if a.matrix.Config.PseudoCondition {
			additionalLines += fmt.Sprintf("PseudoCondition = %.2g\n", a.psudoCondition)
		}
		if a.matrix.Config.Determinant {
			fmt.Printf("Determinant time = %.2f.\n", a.detTime)
			if det != 0.0 && a.detExponent != 0 {
				additionalLines += fmt.Sprintf("Determinant = %.3ge%d\n", det, a.detExponent)
			} else {
				additionalLines += fmt.Sprintf("Determinant = %.3g\n", det)
			}
		}
		if a.matrix.Config.Multiplication && maxRHS != 0.0 {
			additionalLines += fmt.Sprintf("Normalized residual = %.2g\n", normalizedResidual)
		}

		fmt.Printf("\nTotal number of elements = %d\n", a.matrix.ElementCount())
		fmt.Printf("Average number of elements per row initially = %.2f\n", float64(a.matrix.ElementCount()-a.matrix.FillinCount())/float64(a.matrix.GetSize(false)))
		fmt.Printf("Total number of fill-ins = %d\n", a.matrix.Fillins)
		fmt.Println()

		fmt.Print(additionalLines)
	}
	return nil
}

func (a *App) printResourceUsage() {
	var m runtime.MemStats
	runtime.ReadMemStats(&m)

	fmt.Printf("\nAggregate resource usage:\n")
	fmt.Printf("    Time required = %.4f seconds.\n", time.Since(a.startTime).Seconds())
	fmt.Printf("    Heap memory used = %d kBytes\n", m.HeapAlloc/1024)
	fmt.Printf("    Total memory from OS = %d kBytes\n\n", m.Sys/1024)
}

func main() {
	solutionOnly := flag.Bool("s", false, "Print solution rather than run statistics")
	relThreshold := flag.Float64("r", 0.001, "Use x as relative threshold")
	absThreshold := flag.Float64("a", 0.0, "Use x as absolute threshold")
	printLimit := flag.Int("n", 9, "Print first n terms of solution vector")
	iterations := flag.Int("i", 1, "Repeat build/factor/solve n times")
	columnAsRHS := flag.Int("b", -1, "Use n'th column of matrix as b in Ax=b")
	flag.Parse()

	args := flag.Args()
	if len(args) < 1 {
		fmt.Println("Error: Please provide a matrix file")
		os.Exit(1)
	}

	a := InitApp()
	a.filename = args[0]
	a.solutionOnly = *solutionOnly
	a.printLimit = *printLimit
	a.iterations = *iterations

	if *columnAsRHS > 0 {
		a.useColumnAsRHS = true
		a.columnAsRHS = int64(*columnAsRHS)
	}

	fmt.Printf("Sparse Go\nCopyright (c) 2025, Robert Sungwook Shin.\n\n")

	if err := a.readMatrixFromFile(args[0]); err != nil {
		fmt.Printf("%s: %v\n", filepath.Base(os.Args[0]), err)
		os.Exit(1)
	}

	if *relThreshold != 0.001 && *relThreshold > 0.0 {
		a.matrix.RelThreshold = *relThreshold
	}
	a.matrix.AbsThreshold = *absThreshold

	a.matrix.Initialize()

	if a.matrix.Config.ModifiedNodal {
		a.matrix.MNAPreorder()
	}

	if err := a.solve(); err != nil {
		fmt.Printf("%s: %v\n", filepath.Base(os.Args[0]), err)
		os.Exit(1)
	}

	if !a.solutionOnly {
		a.printResourceUsage()
	}

	a.matrix.Destroy()
}
