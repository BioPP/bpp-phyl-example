BIOPP_PATH=$(HOME)/.local

all:
	g++ -I$(BIOPP_PATH)/include -L$(BIOPP_PATH)/lib ExTrees.cpp -lbpp-phyl -lbpp-seq -lbpp-core -o extrees

clean:
	rm extrees
