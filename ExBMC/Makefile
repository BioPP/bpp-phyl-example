BIOPP_PATH=$(HOME)/.local

all:
	g++ -std=c++11 -I$(BIOPP_PATH)/include -L$(BIOPP_PATH)/lib njml.cpp -l bpp-phyl -lbpp-seq -lbpp-core -o njml

clean:
	rm njml
	rm output.dnd
	rm profile.txt
	rm messages.txt
  
