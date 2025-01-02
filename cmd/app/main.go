package main

import (
	"bufio"
	"flag"
	"fmt"
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
)

type App struct {
	matrix       *sparse.Matrix
	filename     string
	description  string
	solutionOnly bool
	printLimit   int
	relThreshold float64
	absThreshold float64
	iterations   int
	columnAsRHS  bool
	rhs          []float64
	irhs         []float64
	solution     []float64
	isolution    []float64
	buildTime    float64
	factorTime   float64
	solveTime    float64
	startTime    time.Time
}

func InitApp() *App {
	return &App{
		absThreshold: 0.0,
		relThreshold: 0.0,
		startTime:    time.Now(),
	}
}

func (t *App) readMatrixFromFile(filename string) error {
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
			t.matrix.Complex = true
		}
		if !scanner.Scan() {
			return fmt.Errorf("missing description")
		}
	}
	t.description = strings.TrimSpace(scanner.Text())

	// Print description at 1st line of matrix
	if !t.solutionOnly {
		fmt.Printf("\n%s\n\n", t.description)
	}
	lineNumber++

	if !scanner.Scan() {
		return fmt.Errorf("missing size information")
	}
	line = scanner.Text()
	if len(strings.TrimSpace(line)) == 0 {
		return fmt.Errorf("syntax error in file '%s' at line %d", t.filename, lineNumber)
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
		ModifiedNodal:           true,
		TiesMultiplier:          5,
		PrinterWidth:            140,
		Annotate:                annotate,
	}

	t.matrix, err = sparse.Create(size, config)
	if err != nil {
		return fmt.Errorf("failed to create matrix: %v", err)
	}

	t.matrix.Clear()

	t.rhs = make([]float64, size+1)
	if isComplex {
		if config.SeparatedComplexVectors {
			t.irhs = make([]float64, size+1)
		} else {
			t.rhs = make([]float64, 2*(size+1))
		}
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

				real, err := strconv.ParseFloat(fields[2], 64)
				if err != nil {
					continue
				}

				element := t.matrix.GetElement(row, col)
				if element != nil {
					element.Real += real

					if isComplex && len(fields) >= 4 {
						if imag, err := strconv.ParseFloat(fields[3], 64); err == nil {
							element.Imag += imag
						}
					}
				}

				if t.columnAsRHS && col == t.matrix.Size {
					t.rhs[row] = real
					if isComplex && len(fields) >= 4 {
						if imag, err := strconv.ParseFloat(fields[3], 64); err == nil {
							t.irhs[row] = imag
						}
					}
				}
			}
		} else {
			rhsValues = append(rhsValues, line)
		}
	}

	// RHS vector
	if !t.columnAsRHS {
		if len(rhsValues) > 0 && strings.HasPrefix(rhsValues[0], "Beginning") {
			rhsValues = rhsValues[1:]
		}

		t.rhs[0] = 0
		t.rhs[1] = 0

		for i := int64(0); i < size && i < int64(len(rhsValues)); i++ {
			fields := strings.Fields(rhsValues[i])
			if len(fields) > 0 {
				if value, err := strconv.ParseFloat(fields[0], 64); err == nil {
					if isComplex && !config.SeparatedComplexVectors {
						t.rhs[2*i+2] = value
						if len(fields) > 1 {
							if imag, err := strconv.ParseFloat(fields[1], 64); err == nil {
								t.rhs[2*i+2+1] = imag
							}
						}
					} else {
						t.rhs[i+1] = value
						if isComplex && len(fields) > 1 {
							if imag, err := strconv.ParseFloat(fields[1], 64); err == nil {
								t.irhs[i+1] = imag
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

	if !t.solutionOnly {
		fmt.Printf("Matrix is %d x %d ", size, size)
		if isComplex {
			fmt.Printf("and complex.\n")
		} else {
			fmt.Printf("and real.\n")
		}
	}

	return nil
}

func (t *App) solve() error {
	orderFactorStart := time.Now()
	if err := t.matrix.OrderAndFactor(t.rhs, t.relThreshold, t.absThreshold, true); err != nil {
		return fmt.Errorf("initial order and factor failed: %v", err)
	}

	initialFactorTime := time.Since(orderFactorStart).Seconds()

	partitionStart := time.Now()
	t.matrix.Partition(sparse.DEFAULT_PARTITION)
	partitionTime := time.Since(partitionStart).Seconds()

	var err error
	if t.matrix.Complex {
		if t.matrix.Config.Transpose {
			t.solution, t.isolution, err = t.matrix.SolveComplexTransposed(t.rhs, t.irhs)
		} else {
			t.solution, t.isolution, err = t.matrix.SolveComplex(t.rhs, t.irhs)
		}
		if err != nil {
			return fmt.Errorf("initial complex solve failed: %v", err)
		}
	} else {
		if t.matrix.Config.Transpose {
			t.solution, err = t.matrix.SolveTransposed(t.rhs)
		} else {
			t.solution, err = t.matrix.Solve(t.rhs)
		}
		if err != nil {
			return fmt.Errorf("initial solve failed: %v", err)
		}
	}

	limit := t.matrix.Size
	if t.printLimit > 0 && int64(t.printLimit) < limit {
		limit = int64(t.printLimit)
	}
	if !t.matrix.Complex {
		fmt.Println("\nSolution:")
		for i := int64(1); i <= limit; i++ {
			fmt.Printf("%-16.9g\n", t.solution[i])
		}
		if limit < t.matrix.Size && limit != 0 {
			fmt.Printf("Solution list truncated.\n")
		}
	} else {
		if t.matrix.Config.Transpose {
			t.solution, t.isolution, err = t.matrix.SolveComplexTransposed(t.rhs, t.irhs)
		} else {
			t.solution, t.isolution, err = t.matrix.SolveComplex(t.rhs, t.irhs)
		}
		if err != nil {
			return fmt.Errorf("initial complex solve failed: %v", err)
		}

		fmt.Println("\nComplex solution:")
		if t.matrix.Config.SeparatedComplexVectors {
			for i := int64(1); i <= limit; i++ {
				fmt.Printf("%-16.9g   %-.9g j\n", t.solution[i], t.isolution[i])
			}
		} else {
			for i := int64(1); i <= limit; i++ {
				fmt.Printf("%-16.9g   %-.9g j\n", t.solution[i*2], t.solution[i*2+1])
			}
		}
	}

	if t.solutionOnly {
		limit := t.matrix.Size
		if t.printLimit > 0 && int64(t.printLimit) < limit {
			limit = int64(t.printLimit)
		}
		fmt.Printf("\nSolution (first %d terms):\n", limit)
		if t.matrix.Complex {
			for i := int64(1); i <= limit; i++ {
				fmt.Printf("%-16.9g   %-.9g j\n", t.solution[i], t.isolution[i])
			}
		} else {
			for i := int64(1); i <= limit; i++ {
				fmt.Printf("%-.9g\n", t.solution[i])
			}
		}
		if limit < t.matrix.Size && limit != 0 {
			fmt.Printf("Solution list truncated.\n")
		}
		fmt.Println()
	} else {
		fmt.Printf("\nStatistics:\n")
		fmt.Printf("Initial factor time = %.2f seconds\n", initialFactorTime)
		fmt.Printf("Partition time = %.2f seconds\n", partitionTime)
		if t.iterations > 0 {
			fmt.Printf("Build time = %.3f seconds\n", t.buildTime/float64(t.iterations))
			fmt.Printf("Factor time = %.3f seconds\n", t.factorTime/float64(t.iterations))
			fmt.Printf("Solve time = %.3f seconds\n", t.solveTime/float64(t.iterations))
		}

		fmt.Printf("\nTotal number of elements = %d\n", t.matrix.ElementCount())
		fmt.Printf("Average number of elements per row initially = %.2f\n", float64(t.matrix.ElementCount()-t.matrix.FillinCount())/float64(t.matrix.GetSize(false)))
		fmt.Printf("Total number of fill-ins = %d\n", t.matrix.Fillins)
	}
	return nil
}

func (t *App) printResourceUsage() {
	var m runtime.MemStats
	runtime.ReadMemStats(&m)

	fmt.Printf("\nAggregate resource usage:\n")
	fmt.Printf("    Time required = %.4f seconds\n", time.Since(t.startTime).Seconds())
	fmt.Printf("    Heap memory used = %d kBytes\n", m.HeapAlloc/1024)
	fmt.Printf("    Total memory from OS = %d kBytes\n\n", m.Sys/1024)
}

func main() {
	solutionOnly := flag.Bool("s", false, "Print solution rather than run statistics")
	relThreshold := flag.Float64("r", 0.001, "Use x as relative threshold")
	absThreshold := flag.Float64("a", 0.0, "Use x as absolute threshold")
	printLimit := flag.Int("n", 9, "Print first n terms of solution vector")
	iterations := flag.Int("i", 1, "Repeat build/factor/solve n times")
	columnAsRHS := flag.Bool("b", false, "Use n'th column of matrix as b in Ax=b")
	flag.Parse()

	args := flag.Args()
	if len(args) < 1 {
		fmt.Println("Error: Please provide a matrix file")
		os.Exit(1)
	}

	app := InitApp()
	app.filename = args[0]
	app.solutionOnly = *solutionOnly
	app.printLimit = *printLimit
	app.iterations = *iterations
	app.columnAsRHS = bool(*columnAsRHS)

	if !app.solutionOnly {
		// fmt.Printf("Sparse1.4\nCopyright (c) 2003, Kenneth S. Kundert.\nAll rights reserved.\n")
		fmt.Printf("Sparse golang\nSungwook Shin\n\n")
	}

	if err := app.readMatrixFromFile(args[0]); err != nil {
		fmt.Printf("%s: %v\n", filepath.Base(os.Args[0]), err)
		os.Exit(1)
	}

	app.matrix.RelThreshold = *relThreshold
	app.matrix.AbsThreshold = *absThreshold

	if app.matrix.Config.ModifiedNodal {
		app.matrix.MNAPreorder()
	}

	if err := app.solve(); err != nil {
		fmt.Printf("%s: %v\n", filepath.Base(os.Args[0]), err)
		os.Exit(1)
	}

	if !app.solutionOnly {
		app.printResourceUsage()
	}

	app.matrix.Destroy()
}
