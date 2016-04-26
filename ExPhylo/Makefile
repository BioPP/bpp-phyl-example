BIOPP_PATH=$(HOME)/.local

all:
	g++ -std=c++11 -I$(BIOPP_PATH)/include -L$(BIOPP_PATH)/lib ExPhylo.cpp -lbpp-phyl -lbpp-seq -lbpp-core -o exphylo

clean:
	rm exphylo
	rm Simulations.fasta
	rm SSU_BioNJ.dnd
	rm SSU_ML.dnd
