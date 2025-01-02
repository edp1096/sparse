.PHONY: default all app factor0 solve1 solve2 ac op clean

default: all
all: app factor0 solve1 solve2 ac op

BINARY_DIR := bin

app:
	go build -o $(BINARY_DIR)/ ./cmd/$@

factor0:
	go build -o $(BINARY_DIR)/ ./cmd/$@

solve1:  
	go build -o $(BINARY_DIR)/ ./cmd/$@

solve2:  
	go build -o $(BINARY_DIR)/ ./cmd/$@

ac:
	go build -o $(BINARY_DIR)/ ./cmd/$@

op:
	go build -o $(BINARY_DIR)/ ./cmd/$@

clean:
	rm -rf $(BINARY_DIR)/*.exe
	rm -rf $(BINARY_DIR)/*.log
