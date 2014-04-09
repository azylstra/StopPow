
INCLUDE = -I../src
DIR = ../src/

UNAME = $(shell uname)
ifeq ($(UNAME), Darwin)
	compiler = clang++
	opts = -c -Wall -std=c++11 -stdlib=libc++ -O3
	linkopts = -std=c++11 -stdlib=libc++ -lgsl
	rm = rm -f
	BIN_FILE_0 = test0.out
	BIN_FILE_1 = test1.out
	BIN_FILE_2 = test2.out
	BIN_FILE_3 = test3.out
	BIN_FILE_4 = test4.out
	BIN_FILE_5 = test5.out
else ifeq ($(UNAME), Linux)
	compiler = g++
	opts = -c -Wall -fPIC -std=c++11 -O3
	linkopts = -lgsl
	rm = rm -f
	BIN_FILE_0 = test0.out
	BIN_FILE_1 = test1.out
	BIN_FILE_2 = test2.out
	BIN_FILE_3 = test3.out
	BIN_FILE_4 = test4.out
	BIN_FILE_5 = test5.out
else # assume Windows
	compiler = g++
	rm = del
	BIN_FILE_0 = test0.exe
	BIN_FILE_1 = test1.exe
	BIN_FILE_2 = test2.exe
	BIN_FILE_3 = test3.exe
	BIN_FILE_4 = test4.exe
	BIN_FILE_5 = test5.exe
endif

objects = StopPow.o StopPow_Plasma.o StopPow_PartialIoniz.o StopPow_LP.o StopPow_BetheBloch.o StopPow_SRIM.o StopPow_Grabowski.o StopPow_AZ.o StopPow_Zimmerman.o StopPow_BPS.o StopPow_Mehlhorn.o PlotGen.o AtomicData.o
BIN_0_O = test0.o
BIN_1_O = test1.o
BIN_2_O = test2.o
BIN_3_O = test3.o
BIN_4_O = test4.o
BIN_5_O = test5.o

test: $(BIN_FILE_0) $(BIN_FILE_1) $(BIN_FILE_2) $(BIN_FILE_3) $(BIN_FILE_4) $(BIN_FILE_5)
	./$(BIN_FILE_0)
	./$(BIN_FILE_1)
	./$(BIN_FILE_2)
	./$(BIN_FILE_3)
	./$(BIN_FILE_4)
	./$(BIN_FILE_5)

test_verbose: $(BIN_FILE_0) $(BIN_FILE_1) $(BIN_FILE_2) $(BIN_FILE_3) $(BIN_FILE_4) $(BIN_FILE_5)
	./$(BIN_FILE_0) --verbose
	./$(BIN_FILE_1) --verbose
	./$(BIN_FILE_2) --verbose
	./$(BIN_FILE_3) --verbose
	./$(BIN_FILE_4) --verbose
	./$(BIN_FILE_5) --verbose
	
$(BIN_FILE_0): $(BIN_0_O) $(objects)
	$(compiler) $(linkopts) -o $(BIN_FILE_0) $(BIN_0_O) $(objects)
	
$(BIN_FILE_1): $(BIN_1_O) $(objects)
	$(compiler) $(linkopts) -o $(BIN_FILE_1) $(BIN_1_O) $(objects)

$(BIN_FILE_2): $(BIN_2_O) $(objects)
	$(compiler) $(linkopts) -o $(BIN_FILE_2) $(BIN_2_O) $(objects)

$(BIN_FILE_3): $(BIN_3_O) $(objects)
	$(compiler) $(linkopts) -o $(BIN_FILE_3) $(BIN_3_O) $(objects)

$(BIN_FILE_4): $(BIN_4_O) $(objects)
	$(compiler) $(linkopts) -o $(BIN_FILE_4) $(BIN_4_O) $(objects)

$(BIN_FILE_5): $(BIN_5_O) $(objects)
	$(compiler) $(linkopts) -o $(BIN_FILE_5) $(BIN_5_O) $(objects)

$(BIN_0_O): test0.cpp
	$(compiler) $(opts) $(INCLUDE) test0.cpp

$(BIN_1_O): test1.cpp
	$(compiler) $(opts) $(INCLUDE) test1.cpp

$(BIN_2_O): test2.cpp
	$(compiler) $(opts) $(INCLUDE) test2.cpp

$(BIN_3_O): test3.cpp
	$(compiler) $(opts) $(INCLUDE) test3.cpp

$(BIN_4_O): test4.cpp
	$(compiler) $(opts) $(INCLUDE) test4.cpp

$(BIN_5_O): test5.cpp
	$(compiler) $(opts) $(INCLUDE) test5.cpp
	
StopPow.o: $(DIR)StopPow.cpp 
	$(compiler) $(opts) $(INCLUDE) $(DIR)StopPow.cpp
	
StopPow_Plasma.o: $(DIR)StopPow_Plasma.cpp 
	$(compiler) $(opts) $(INCLUDE) $(DIR)StopPow_Plasma.cpp
	
StopPow_PartialIoniz.o: $(DIR)StopPow_PartialIoniz.cpp 
	$(compiler) $(opts) $(INCLUDE) $(DIR)StopPow_PartialIoniz.cpp

StopPow_LP.o: $(DIR)StopPow.cpp $(DIR)StopPow_LP.cpp
	$(compiler) $(opts) $(INCLUDE) $(DIR)StopPow.cpp $(DIR)StopPow_LP.cpp

StopPow_BetheBloch.o: $(DIR)StopPow.cpp $(DIR)StopPow_BetheBloch.cpp
	$(compiler) $(opts) $(INCLUDE) $(DIR)StopPow.cpp $(DIR)StopPow_BetheBloch.cpp

StopPow_SRIM.o: $(DIR)StopPow.cpp $(DIR)StopPow_SRIM.cpp
	$(compiler) $(opts) $(INCLUDE) $(DIR)StopPow.cpp $(DIR)StopPow_SRIM.cpp

StopPow_Grabowski.o: $(DIR)StopPow.cpp $(DIR)StopPow_Grabowski.cpp
	$(compiler) $(opts) $(INCLUDE) $(DIR)StopPow.cpp $(DIR)StopPow_Grabowski.cpp

StopPow_AZ.o: $(DIR)StopPow.cpp $(DIR)StopPow_AZ.cpp
	$(compiler) $(opts) $(INCLUDE) $(DIR)StopPow.cpp $(DIR)StopPow_AZ.cpp

StopPow_Zimmerman.o: $(DIR)StopPow.cpp $(DIR)StopPow_Zimmerman.cpp
	$(compiler) $(opts) $(INCLUDE) $(DIR)StopPow.cpp $(DIR)StopPow_Zimmerman.cpp

StopPow_BPS.o: $(DIR)StopPow.cpp $(DIR)StopPow_BPS.cpp
	$(compiler) $(opts) $(INCLUDE) $(DIR)StopPow.cpp $(DIR)StopPow_BPS.cpp

StopPow_Mehlhorn.o: $(DIR)StopPow.cpp $(DIR)StopPow_Mehlhorn.cpp
	$(compiler) $(opts) $(INCLUDE) $(DIR)StopPow.cpp $(DIR)StopPow_Mehlhorn.cpp

PlotGen.o: $(DIR)PlotGen.cpp
	$(compiler) $(opts) $(INCLUDE) $(DIR)PlotGen.cpp

AtomicData.o: $(DIR)AtomicData.cpp
	$(compiler) $(opts) $(INCLUDE) $(DIR)AtomicData.cpp

clean:
	$(rm) $(objects) $(BIN_0_O) $(BIN_1_O) $(BIN_2_O) $(BIN_3_O) $(BIN_4_O) $(BIN_5_O) $(BIN_FILE_0) $(BIN_FILE_1) $(BIN_FILE_2) $(BIN_FILE_3) $(BIN_FILE_4) $(BIN_FILE_5)
	