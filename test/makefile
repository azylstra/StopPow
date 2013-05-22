
INCLUDE = -I../src
DIR = ../src/

UNAME = $(shell uname)
ifeq ($(UNAME), Darwin)
	compiler = clang++
	opts = -c -Wall -lib=std-c++11 -std=c++11
	rm = rm
	BIN_FILE_1 = test1
else ifeq ($(UNAME), Linux)
	compiler = g++
	opts = -c -Wall -fPIC -std=c++11
	rm = rm
	BIN_FILE_1 = test1
else # assume Windows
	compiler = g++
	rm = del
	BIN_FILE_1 = test1.exe
endif

objects = test1.o StopPow.o StopPow_LP.o StopPow_BetheBloch.o StopPow_SRIM.o

test: $(BIN_FILE_1)
	./$(BIN_FILE_1)

test_verbose: $(BIN_FILE_1)
	./$(BIN_FILE_1) --verbose
	
$(BIN_FILE_1): $(objects)
	$(compiler) -o $(BIN_FILE_1) $(objects)
	$(rm) $(objects)

test1.o: test1.cpp
	$(compiler) $(opts) $(INCLUDE) test1.cpp
	
StopPow.o: $(DIR)StopPow.cpp 
	$(compiler) $(opts) $(INCLUDE) $(DIR)StopPow.cpp

StopPow_LP.o: $(DIR)StopPow.cpp $(DIR)StopPow_LP.cpp
	$(compiler) $(opts) $(INCLUDE) $(DIR)StopPow.cpp $(DIR)StopPow_LP.cpp

StopPow_BetheBloch.o: $(DIR)StopPow.cpp $(DIR)StopPow_BetheBloch.cpp
	$(compiler) $(opts) $(INCLUDE) $(DIR)StopPow.cpp $(DIR)StopPow_BetheBloch.cpp

StopPow_SRIM.o: $(DIR)StopPow.cpp $(DIR)StopPow_SRIM.cpp
	$(compiler) $(opts) $(INCLUDE) $(DIR)StopPow.cpp $(DIR)StopPow_SRIM.cpp

clean:
	rm $(BIN_FILE)