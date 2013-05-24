# name of binaries for swig
SWIG = swig

SWIG_FILE = StopPow.i
SWIG_OPTS = -java -package cStopPow -c++
SWIG_WRAP_SRC = StopPow_wrap.cxx

DEST_DIR = StopPow

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

SRC_FILES = ../src/StopPow.cpp ../src/StopPow_SRIM.cpp ../src/StopPow_LP.cpp ../src/StopPow_BetheBloch.cpp ../src/PlotGen.cpp
OBJ_FILES = StopPow_wrap.o StopPow.o StopPow_SRIM.o StopPow_LP.o StopPow_BetheBloch.o PlotGen.o

OUT_FILE = libcStopPow.jnilib

JAR_TEMP_DIR = cStopPow
JAR_FILE = libcStopPow.jar

StopPow: 
	$(SWIG) $(SWIG_OPTS) $(SWIG_FILE)
	$(compiler) $(JAVA_INCLUDE) $(CC_OPTS) -c $(SRC_FILES) $(SWIG_WRAP_SRC)
	$(compiler) $(L_OPTS) $(OBJ_FILES) -o $(OUT_FILE)
	mkdir -p $(DEST_DIR)
	mv *.java $(DEST_DIR)
	mv $(OUT_FILE) $(DEST_DIR)
	$(rm) $(OBJ_FILES)
	$(rm) $(SWIG_WRAP_SRC)
	mkdir -p $(JAR_TEMP_DIR)
	cp -r $(DEST_DIR)/* $(JAR_TEMP_DIR)
	javac $(JAR_TEMP_DIR)/*.java
	$(rm) $(JAR_TEMP_DIR)/*.java
	jar cf $(JAR_FILE) $(JAR_TEMP_DIR)
	$(rm) -r $(JAR_TEMP_DIR)
	
	
clean:
	$(rm) -r $(DEST_DIR)
	$(rm) $(JAR_FILE)