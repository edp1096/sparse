.PHONY: default all sparse factor1 solve1 solve2 op tran1 tran2 ac1 clean

default: all
all: sparse factor1 solve1 solve2 op tran1 tran2 ac1

BINARY_DIR := bin

sparse:
	go build -o $(BINARY_DIR)/ ./cmd/$@

factor1:
	go build -o $(BINARY_DIR)/ ./cmd/$@

solve1:  
	go build -o $(BINARY_DIR)/ ./cmd/$@

solve2:  
	go build -o $(BINARY_DIR)/ ./cmd/$@

op:
	go build -o $(BINARY_DIR)/ ./cmd/$@

tran1:
	go build -o $(BINARY_DIR)/ ./cmd/$@

tran2:
	go build -o $(BINARY_DIR)/ ./cmd/$@

ac1:
	go build -o $(BINARY_DIR)/ ./cmd/$@

clean:
	rm -rf $(BINARY_DIR)/*.exe
	rm -rf $(BINARY_DIR)/*.log
