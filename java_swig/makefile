# name of binaries for swig
SWIG = swig

SWIG_FILE = StopPow.i
SWIG_OPTS = -java -c++
SWIG_WRAP_SRC = StopPow_wrap.cxx

DEST_DIR = dist

# system dependent options:
UNAME = $(shell uname)
ifeq ($(UNAME), Darwin)
	rm = rm
	ext = so
	compiler = clang++
	JAVA_INCLUDE = -I/System/Library/Frameworks/JavaVM.framework/Versions/Current/Headers/
	CC_OPTS = 
	L_OPTS = -framework JavaVM -bundle
else ifeq ($(UNAME), Linux)
	rm = rm
	ext = so
	compiler = g++
	JAVA_INCLUDE = -I/usr/lib/jvm/java-7-openjdk-amd64/include/
	CC_OPTS = -fPIC
	L_OPTS = -shared
else # Windows
	rm = del
	ext = dll
	compiler = g++
	JAVA_INCLUDE = 
	CC_OPTS = 
endif

SRC_FILES = ../src/StopPow.cpp ../src/StopPow_SRIM.cpp ../src/StopPow_LP.cpp ../src/StopPow_BetheBloch.cpp
OBJ_FILES = StopPow_wrap.o StopPow.o StopPow_SRIM.o StopPow_LP.o StopPow_BetheBloch.o 

OUT_FILE = libcStopPow.jnilib

StopPow: 
	$(SWIG) $(SWIG_OPTS) $(SWIG_FILE)
	$(compiler) $(JAVA_INCLUDE) $(CC_OPTS) -c $(SRC_FILES) $(SWIG_WRAP_SRC)
	$(compiler) $(L_OPTS) $(OBJ_FILES) -o $(OUT_FILE)
	mkdir -p $(DEST_DIR)
	mv *.java $(DEST_DIR)
	mv $(OUT_FILE) $(DEST_DIR)
	$(rm) $(OBJ_FILES)
	$(rm) $(SWIG_WRAP_SRC)
	
clean:
	$(rm) -r $(DEST_DIR)