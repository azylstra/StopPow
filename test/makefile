
INCLUDE = -I../src
DIR = ../src/

UNAME = $(shell uname)
ifeq ($(UNAME), Darwin)
	compiler = clang++
	opts = -c -Wall -std=c++11 -stdlib=libc++ -O3
	linkopts = -std=c++11 -stdlib=libc++
	rm = rm
	BIN_FILE_1 = test1.out
	BIN_FILE_2 = test2.out
	BIN_FILE_3 = test3.out
else ifeq ($(UNAME), Linux)
	compiler = g++
	opts = -c -Wall -fPIC -std=c++11 -O3
	rm = rm
	BIN_FILE_1 = test1.out
	BIN_FILE_2 = test2.out
	BIN_FILE_3 = test3.out
else # assume Windows
	compiler = g++
	rm = del
	BIN_FILE_1 = test1.exe
	BIN_FILE_2 = test2.exe
	BIN_FILE_3 = test3.exe
endif

objects = StopPow.o StopPow_LP.o StopPow_BetheBloch.o StopPow_SRIM.o PlotGen.o AtomicData.o
BIN_1_O = test1.o
BIN_2_O = test2.o
BIN_3_O = test3.o

test: $(BIN_FILE_1) $(BIN_FILE_2) $(BIN_FILE_3)
	./$(BIN_FILE_1)
	./$(BIN_FILE_2)
	./$(BIN_FILE_3)

test_verbose: $(BIN_FILE_1) $(BIN_FILE_2) $(BIN_FILE_3)
	./$(BIN_FILE_1) --verbose
	./$(BIN_FILE_2) --verbose
	./$(BIN_FILE_3) --verbose
	
$(BIN_FILE_1): $(BIN_1_O) $(objects)
	$(compiler) $(linkopts) -o $(BIN_FILE_1) $(BIN_1_O) $(objects)

$(BIN_FILE_2): $(BIN_2_O) $(objects)
	$(compiler) $(linkopts) -o $(BIN_FILE_2) $(BIN_2_O) $(objects)

$(BIN_FILE_3): $(BIN_3_O) $(objects)
	$(compiler) $(linkopts) -o $(BIN_FILE_3) $(BIN_3_O) $(objects)

$(BIN_1_O): test1.cpp
	$(compiler) $(opts) $(INCLUDE) test1.cpp

$(BIN_2_O): test2.cpp
	$(compiler) $(opts) $(INCLUDE) test2.cpp

$(BIN_3_O): test3.cpp
	$(compiler) $(opts) $(INCLUDE) test3.cpp
	
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

AtomicData.o: $(DIR)AtomicData.cpp
	$(compiler) $(opts) $(INCLUDE) $(DIR)AtomicData.cpp

clean:
	$(rm) $(objects) $(BIN_1_O) $(BIN_2_O) $(BIN_3_O)
	rm $(BIN_FILE_1)
	rm $(BIN_FILE_2)
	rm $(BIN_FILE_3)