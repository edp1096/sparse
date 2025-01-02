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
	buildTime      float64
	factorTime     float64
	solveTime      float64
	startTime      time.Time
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
		ModifiedNodal:           true,
		TiesMultiplier:          5,
		PrinterWidth:            140,
		Annotate:                annotate,
	}

	a.matrix, err = sparse.Create(size, config)
	if err != nil {
		return fmt.Errorf("failed to create matrix: %v", err)
	}

	a.matrix.Clear()

	a.rhs = make([]float64, size+1)
	if isComplex {
		if config.SeparatedComplexVectors {
			a.irhs = make([]float64, size+1)
		} else {
			a.rhs = make([]float64, 2*(size+1))
		}
	}

	// RHS_Col 설정
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
	orderFactorStart := time.Now()
	if err := a.matrix.OrderAndFactor(a.rhs, a.relThreshold, a.absThreshold, true); err != nil {
		return fmt.Errorf("initial order and factor failed: %v", err)
	}

	initialFactorTime := time.Since(orderFactorStart).Seconds()

	partitionStart := time.Now()
	a.matrix.Partition(sparse.DEFAULT_PARTITION)
	partitionTime := time.Since(partitionStart).Seconds()

	var err error
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

	limit := a.matrix.Size
	if a.printLimit > 0 && int64(a.printLimit) < limit {
		limit = int64(a.printLimit)
	}
	if !a.matrix.Complex {
		fmt.Println("\nSolution:")
		for i := int64(1); i <= limit; i++ {
			fmt.Printf("%-16.9g\n", a.solution[i])
		}
		if limit < a.matrix.Size && limit != 0 {
			fmt.Printf("Solution list truncated.\n")
		}
	} else {
		if a.matrix.Config.Transpose {
			a.solution, a.isolution, err = a.matrix.SolveComplexTransposed(a.rhs, a.irhs)
		} else {
			a.solution, a.isolution, err = a.matrix.SolveComplex(a.rhs, a.irhs)
		}
		if err != nil {
			return fmt.Errorf("initial complex solve failed: %v", err)
		}

		fmt.Println("\nComplex solution:")
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

	if a.solutionOnly {
		limit := a.matrix.Size
		if a.printLimit > 0 && int64(a.printLimit) < limit {
			limit = int64(a.printLimit)
		}
		fmt.Printf("\nSolution (first %d terms):\n", limit)
		if a.matrix.Complex {
			for i := int64(1); i <= limit; i++ {
				fmt.Printf("%-16.9g   %-.9g j\n", a.solution[i], a.isolution[i])
			}
		} else {
			for i := int64(1); i <= limit; i++ {
				fmt.Printf("%-.9g\n", a.solution[i])
			}
		}
		if limit < a.matrix.Size && limit != 0 {
			fmt.Printf("Solution list truncated.\n")
		}
		fmt.Println()
	} else {
		fmt.Printf("\nStatistics:\n")
		fmt.Printf("Initial factor time = %.2f seconds\n", initialFactorTime)
		fmt.Printf("Partition time = %.2f seconds\n", partitionTime)
		if a.iterations > 0 {
			fmt.Printf("Build time = %.3f seconds\n", a.buildTime/float64(a.iterations))
			fmt.Printf("Factor time = %.3f seconds\n", a.factorTime/float64(a.iterations))
			fmt.Printf("Solve time = %.3f seconds\n", a.solveTime/float64(a.iterations))
		}

		fmt.Printf("\nTotal number of elements = %d\n", a.matrix.ElementCount())
		fmt.Printf("Average number of elements per row initially = %.2f\n", float64(a.matrix.ElementCount()-a.matrix.FillinCount())/float64(a.matrix.GetSize(false)))
		fmt.Printf("Total number of fill-ins = %d\n", a.matrix.Fillins)
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
	app.useColumnAsRHS = bool(*columnAsRHS)

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
