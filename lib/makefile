
INCLUDE = -I../src
DIR = ../src/
INCLUDES_SRC = ../src/*.h
#LIB = -L/usr/local/Cellar/boost/1.53.0/lib

UNAME = $(shell uname)
ifeq ($(UNAME), Darwin)
	compiler = clang++
	#opts = -c -Wall  -I/usr/local/Cellar/boost/1.53.0/include -lboost_math
	opts = -c -Wall
	LIB_FILE = StopPow.so
	cp = cp
	mv = mv
	mkdir = mkdir -p
	rm = rm
	rm_dir = rm -r
else ifeq ($(UNAME), Linux)
	compiler = g++
	opts = -c -Wall -fPIC
	LIB_FILE = StopPow.so
	cp = cp
	mv = mv
	mkdir = mkdir -p
	rm = rm
	rm_dir = rm -r
else # assume Windows
	compiler = g++
	opts = -c -Wall
	LIB_FILE = StopPow.dll
	mkdir = mkdir
	cp = copy
	mv = move
	rm = del
	rm_dir = rmdir /s /q
	INCLUDES_SRC = ..\src\*.h
endif

objects = StopPow.o StopPow_LP.o StopPow_BetheBloch.o StopPow_SRIM.o

INCLUDES_DST = include
LIB_DST = lib

$(LIB_FILE): $(objects)
	$(compiler) -o $(LIB_FILE) -shared $(objects)
	$(mkdir) $(INCLUDES_DST)
	$(cp) $(INCLUDES_SRC) $(INCLUDES_DST)
	$(mkdir) $(LIB_DST)
	$(rm) $(objects)
	$(mv) $(LIB_FILE) $(LIB_DST)
	
StopPow.o: $(DIR)StopPow.cpp 
	$(compiler) $(opts) $(INCLUDE) $(DIR)StopPow.cpp

StopPow_LP.o: $(DIR)StopPow.cpp $(DIR)StopPow_LP.cpp
	$(compiler) $(opts) $(INCLUDE) $(DIR)StopPow.cpp $(DIR)StopPow_LP.cpp

StopPow_BetheBloch.o: $(DIR)StopPow.cpp $(DIR)StopPow_BetheBloch.cpp
	$(compiler) $(opts) $(INCLUDE) $(DIR)StopPow.cpp $(DIR)StopPow_BetheBloch.cpp

StopPow_SRIM.o: $(DIR)StopPow.cpp $(DIR)StopPow_SRIM.cpp
	$(compiler) $(opts) $(INCLUDE) $(DIR)StopPow.cpp $(DIR)StopPow_SRIM.cpp

clean:
	$(rm_dir) $(INCLUDES_DST)
	$(rm_dir) $(LIB_DST)