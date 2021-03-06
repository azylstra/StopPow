
INCLUDE = -I../src
DIR = ../src/

UNAME = $(shell uname)
ifeq ($(UNAME), Darwin)
	compiler = clang++
	opts = -c  -I/usr/local/include -Wall -std=c++11 -stdlib=libc++ -g
	linkopts = -L/usr/local/lib -std=c++11 -stdlib=libc++ -lgsl -lgslcblas
	rm = rm -f
	BIN_FILE_0 = test0.out
	BIN_FILE_1 = test1.out
	BIN_FILE_2 = test2.out
	BIN_FILE_3 = test3.out
	BIN_FILE_4 = test4.out
	BIN_FILE_5 = test5.out
	BIN_FILE_6 = test6.out
	BIN_FILE_7 = test7.out
	BIN_FILE_8 = test8.out
else ifeq ($(UNAME), Linux)
	compiler = g++
	opts = -c -Wall -fPIC -std=c++11 -O3
	linkopts = -lgsl -lgslcblas
	rm = rm -f
	BIN_FILE_0 = test0.out
	BIN_FILE_1 = test1.out
	BIN_FILE_2 = test2.out
	BIN_FILE_3 = test3.out
	BIN_FILE_4 = test4.out
	BIN_FILE_5 = test5.out
	BIN_FILE_6 = test6.out
	BIN_FILE_7 = test7.out
	BIN_FILE_8 = test8.out
else # assume Windows
	compiler = g++
	rm = del
	BIN_FILE_0 = test0.exe
	BIN_FILE_1 = test1.exe
	BIN_FILE_2 = test2.exe
	BIN_FILE_3 = test3.exe
	BIN_FILE_4 = test4.exe
	BIN_FILE_5 = test5.exe
	BIN_FILE_6 = test6.exe
	BIN_FILE_7 = test7.exe
	BIN_FILE_8 = test8.exe
endif

objects = StopPow.o StopPow_Plasma.o StopPow_PartialIoniz.o StopPow_LP.o StopPow_BetheBloch.o StopPow_SRIM.o StopPow_Grabowski.o StopPow_AZ.o StopPow_Zimmerman.o StopPow_BPS.o StopPow_Mehlhorn.o StopPow_Fit.o PlotGen.o AtomicData.o Spectrum.o Fit.o
BIN_0_O = test0.o
BIN_1_O = test1.o
BIN_2_O = test2.o
BIN_3_O = test3.o
BIN_4_O = test4.o
BIN_5_O = test5.o
BIN_6_O = test6.o
BIN_7_O = test7.o
BIN_8_O = test8.o

test: $(BIN_FILE_0) $(BIN_FILE_1) $(BIN_FILE_2) $(BIN_FILE_3) $(BIN_FILE_4) $(BIN_FILE_5) $(BIN_FILE_6) $(BIN_FILE_7)
	./$(BIN_FILE_0)
	./$(BIN_FILE_1)
	./$(BIN_FILE_2)
	./$(BIN_FILE_3)
	./$(BIN_FILE_4)
	./$(BIN_FILE_5)
	./$(BIN_FILE_6)
	./$(BIN_FILE_7)
	./$(BIN_FILE_8)

test_verbose: $(BIN_FILE_0) $(BIN_FILE_1) $(BIN_FILE_2) $(BIN_FILE_3) $(BIN_FILE_4) $(BIN_FILE_5) $(BIN_FILE_6) $(BIN_FILE_7)
	./$(BIN_FILE_0) --verbose
	./$(BIN_FILE_1) --verbose
	./$(BIN_FILE_2) --verbose
	./$(BIN_FILE_3) --verbose
	./$(BIN_FILE_4) --verbose
	./$(BIN_FILE_5) --verbose
	./$(BIN_FILE_6) --verbose
	./$(BIN_FILE_7) --verbose
	./$(BIN_FILE_8) --verbose
	
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

$(BIN_FILE_6): $(BIN_6_O) $(objects)
	$(compiler) $(linkopts) -o $(BIN_FILE_6) $(BIN_6_O) $(objects)

$(BIN_FILE_7): $(BIN_7_O) $(objects)
	$(compiler) $(linkopts) -o $(BIN_FILE_7) $(BIN_7_O) $(objects)

$(BIN_FILE_8): $(BIN_8_O) $(objects)
	$(compiler) $(linkopts) -o $(BIN_FILE_8) $(BIN_8_O) $(objects)

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

$(BIN_6_O): test6.cpp
	$(compiler) $(opts) $(INCLUDE) test6.cpp

$(BIN_7_O): test7.cpp
	$(compiler) $(opts) $(INCLUDE) test7.cpp

$(BIN_8_O): test8.cpp
	$(compiler) $(opts) $(INCLUDE) test8.cpp
	
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

StopPow_Fit.o: $(DIR)StopPow.cpp $(DIR)StopPow_Fit.cpp
	$(compiler) $(opts) $(INCLUDE) $(DIR)StopPow.cpp $(DIR)StopPow_Fit.cpp

PlotGen.o: $(DIR)PlotGen.cpp
	$(compiler) $(opts) $(INCLUDE) $(DIR)PlotGen.cpp

AtomicData.o: $(DIR)AtomicData.cpp
	$(compiler) $(opts) $(INCLUDE) $(DIR)AtomicData.cpp

Spectrum.o: $(DIR)Spectrum.cpp
	$(compiler) $(opts) $(INCLUDE) $(DIR)Spectrum.cpp

Fit.o: $(DIR)Fit.cpp
	$(compiler) $(opts) $(INCLUDE) $(DIR)Fit.cpp

clean:
	$(rm) $(objects) $(BIN_0_O) $(BIN_1_O) $(BIN_2_O) $(BIN_3_O) $(BIN_4_O) $(BIN_5_O) $(BIN_6_O) $(BIN_7_O) $(BIN_8_O) $(BIN_FILE_0) $(BIN_FILE_1) $(BIN_FILE_2) $(BIN_FILE_3) $(BIN_FILE_4) $(BIN_FILE_5) $(BIN_FILE_6) $(BIN_FILE_7) $(BIN_FILE_8)
	