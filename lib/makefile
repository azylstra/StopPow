
INCLUDE = -I../src
DIR = ../src/
INCLUDES_SRC = ../src/*.h

UNAME = $(shell uname)
ifeq ($(UNAME), Darwin)
	compiler = clang++
	opts = -c -Wall -std=c++11 -stdlib=libc++ -O3
	linker = libtool
	link_opts = -static
	LIB_FILE = StopPow.a
	cp = cp
	mv = mv
	mkdir = mkdir -p
	rm = rm -f
	rm_dir = rm -rf
else ifeq ($(UNAME), Linux)
	compiler = g++
	opts = -c -Wall -fPIC -std=c++11 -O3
	linker = ar
	link_opts = rvs
	LIB_FILE = StopPow.a
	cp = cp
	mv = mv
	mkdir = mkdir -p
	rm = rm -f
	rm_dir = rm -rf
else ifeq ($(UNAME), CYGWIN_NT-6.1-WOW64)
	compiler = g++
	opts = -c -Wall -fPIC -std=c++11 -O3
	linker = ar
	link_opts = rvs
	LIB_FILE = StopPow.a
	cp = cp
	mv = mv
	mkdir = mkdir -p
	rm = rm -f
	rm_dir = rm -rf
else # assume Windows
	compiler = g++
	opts = -c -Wall -fPIC -std=c++11 -O3
	linker = ar
	link_opts = rvs
	LIB_FILE = StopPow.a
	mkdir = mkdir
	cp = copy
	mv = move
	rm = del
	rm_dir = rmdir /s /q
	INCLUDES_SRC = ..\src\*.h
endif

objects = StopPow.o StopPow_LP.o StopPow_BetheBloch.o StopPow_SRIM.o StopPow_AZ.o StopPow_Mehlhorn.o PlotGen.o AtomicData.o Spectrum.o

INCLUDES_DST = include
LIB_DST = lib

$(LIB_FILE): $(objects)
	$(linker) $(link_opts) -o $(LIB_FILE) $(objects)
	$(mkdir) $(INCLUDES_DST)
	$(cp) $(INCLUDES_SRC) $(INCLUDES_DST)
	$(mkdir) $(LIB_DST)
	$(rm) $(objects)
	$(mv) $(LIB_FILE) $(LIB_DST)$(slash)
	
StopPow.o: $(DIR)StopPow.cpp 
	$(compiler) $(opts) $(INCLUDE) $(DIR)StopPow.cpp

StopPow_LP.o: $(DIR)StopPow.cpp $(DIR)StopPow_LP.cpp
	$(compiler) $(opts) $(INCLUDE) $(DIR)StopPow.cpp $(DIR)StopPow_LP.cpp

StopPow_BetheBloch.o: $(DIR)StopPow.cpp $(DIR)StopPow_BetheBloch.cpp
	$(compiler) $(opts) $(INCLUDE) $(DIR)StopPow.cpp $(DIR)StopPow_BetheBloch.cpp

StopPow_SRIM.o: $(DIR)StopPow.cpp $(DIR)StopPow_SRIM.cpp
	$(compiler) $(opts) $(INCLUDE) $(DIR)StopPow.cpp $(DIR)StopPow_SRIM.cpp

StopPow_AZ.o: $(DIR)StopPow.cpp $(DIR)StopPow_AZ.cpp
	$(compiler) $(opts) $(INCLUDE) $(DIR)StopPow.cpp $(DIR)StopPow_AZ.cpp

StopPow_Mehlhorn.o: $(DIR)StopPow.cpp $(DIR)StopPow_Mehlhorn.cpp
	$(compiler) $(opts) $(INCLUDE) $(DIR)StopPow.cpp $(DIR)StopPow_Mehlhorn.cpp

PlotGen.o: $(DIR)PlotGen.cpp
	$(compiler) $(opts) $(INCLUDE) $(DIR)PlotGen.cpp

AtomicData.o: $(DIR)AtomicData.cpp
	$(compiler) $(opts) $(INCLUDE) $(DIR)AtomicData.cpp

Spectrum.o: $(DIR)Spectrum.cpp
	$(compiler) $(opts) $(INCLUDE) $(DIR)Spectrum.cpp

clean:
	$(rm) *.o
	$(rm_dir) $(INCLUDES_DST)
	$(rm_dir) $(LIB_DST)