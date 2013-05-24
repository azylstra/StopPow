
INCLUDE = -I../src
DIR = ../src/

UNAME = $(shell uname)
ifeq ($(UNAME), Darwin)
	compiler = clang++
	opts = -c -Wall -lib=std-c++11 -std=c++11
	rm = rm
	BIN_FILE_1 = test1.out
	BIN_FILE_2 = test2.out
else ifeq ($(UNAME), Linux)
	compiler = g++
	opts = -c -Wall -fPIC -std=c++11
	rm = rm
	BIN_FILE_1 = test1.out
	BIN_FILE_2 = test2.out
else # assume Windows
	compiler = g++
	rm = del
	BIN_FILE_1 = test1.exe
	BIN_FILE_2 = test2.exe
endif

objects = StopPow.o StopPow_LP.o StopPow_BetheBloch.o StopPow_SRIM.o PlotGen.o
BIN_1_O = test1.o
BIN_2_O = test2.o

test: $(BIN_FILE_1) $(BIN_FILE_2)
	./$(BIN_FILE_1)
	./$(BIN_FILE_2)

test_verbose: $(BIN_FILE_1) $(BIN_FILE_2)
	./$(BIN_FILE_1) --verbose
	./$(BIN_FILE_2) --verbose
	
$(BIN_FILE_1): $(BIN_1_O) $(objects)
	$(compiler) -o $(BIN_FILE_1) $(BIN_1_O) $(objects)

$(BIN_FILE_2): $(BIN_2_O) $(objects)
	$(compiler) -o $(BIN_FILE_2) $(BIN_2_O) $(objects)

$(BIN_1_O): test1.cpp
	$(compiler) $(opts) $(INCLUDE) test1.cpp

$(BIN_2_O): test2.cpp
	$(compiler) $(opts) $(INCLUDE) test2.cpp
	
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
	$(rm) $(objects) $(BIN_1_O) $(BIN_2_O)
	rm $(BIN_FILE_1)
	rm $(BIN_FILE_2)