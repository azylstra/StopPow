
INCLUDE = -I../../src
DIR = ../../src/

UNAME = $(shell uname)
ifeq ($(UNAME), Darwin)
	compiler = clang++
	opts = -c -Wall -lib=std-c++11 -std=c++11 -O3
	rm = rm
	BIN_FILE = StopPow
else ifeq ($(UNAME), Linux)
	compiler = g++
	rm = rm
	BIN_FILE = StopPow
else # assume Windows
	compiler = g++
	rm = del
	BIN_FILE = StopPow.exe
endif

objects = main.o StopPow.o StopPow_LP.o StopPow_BetheBloch.o StopPow_SRIM.o PlotGen.o


$(BIN_FILE): $(objects)
	$(compiler) -o $(BIN_FILE) $(objects)
	$(rm) $(objects)

main.o: main.cpp
	$(compiler) $(opts) $(INCLUDE) main.cpp
	
StopPow.o: $(DIR)StopPow.cpp 
	$(compiler) $(opts) $(INCLUDE) $(DIR)StopPow.cpp

StopPow_LP.o: $(DIR)StopPow.cpp $(DIR)StopPow_LP.cpp
	$(compiler) $(opts) $(INCLUDE) $(DIR)StopPow.cpp $(DIR)StopPow_LP.cpp

StopPow_BetheBloch.o: $(DIR)StopPow.cpp $(DIR)StopPow_BetheBloch.cpp
	$(compiler) $(opts) $(INCLUDE) $(DIR)StopPow.cpp $(DIR)StopPow_BetheBloch.cpp

StopPow_SRIM.o: $(DIR)StopPow.cpp $(DIR)StopPow_SRIM.cpp
	$(compiler) $(opts) $(INCLUDE) $(DIR)StopPow.cpp $(DIR)StopPow_SRIM.cpp

PlotGen.o: $(DIR)PlotGen.cpp
	$(compiler) $(opts) $(INCLUDE) $(DIR)PlotGen.cpp

clean:
	rm $(BIN_FILE)