# name of binaries for swig
SWIG = swig

SWIG_FILE = StopPow.i
SWIG_OPTS = -java -package cStopPow -c++
SWIG_WRAP_SRC = StopPow_wrap.cxx

DEST_DIR = StopPow

SRC_FILES = ../src/StopPow.cpp ../src/StopPow_SRIM.cpp ../src/StopPow_LP.cpp ../src/StopPow_BetheBloch.cpp ../src/StopPow_AZ.cpp ../src/PlotGen.cpp ../src/AtomicData.cpp
OBJ_FILES = StopPow_wrap.o StopPow.o StopPow_SRIM.o StopPow_LP.o StopPow_BetheBloch.o StopPow_AZ.o PlotGen.o AtomicData.o

OUT_FILE = libcStopPow.jnilib

JAR_TEMP_DIR = cStopPow
JAR_FILE = libcStopPow.jar

# system dependent options:
UNAME = $(shell uname)
ifeq ($(UNAME), Darwin)
	ext = so
	compiler = clang++
	JAVA_INCLUDE = -I/System/Library/Frameworks/JavaVM.framework/Versions/Current/Headers/
	CC_OPTS = -c -Wall -std=c++11 -stdlib=libc++ -O3
	L_OPTS = -std=c++11 -stdlib=libc++ -framework JavaVM -bundle
	cp = cp
	mv = mv
	mkdir = mkdir -p
	rm = rm -f
	rm_dir = rm -rf
else ifeq ($(UNAME), Linux)
	ext = so
	compiler = g++
	linker = g++
	JAVA_INCLUDE = -I/usr/lib/jvm/java-7-oracle/include/ -I/usr/lib/jvm/java-7-oracle/include/linux/
	CC_OPTS = -fPIC -std=c++11 -O3
	L_OPTS = -shared
	cp = cp
	mv = mv
	mkdir = mkdir -p
	rm = rm -f
	rm_dir = rm -rf
	wincomp = i686-w64-mingw32-g++
	WIN_CC_OPTS = -std=c++0x -O3 -D_JNI_IMPLEMENTATION_ -m32 -mwindows -I/usr/lib/jvm/java-7-oracle/include/ -I/usr/lib/jvm/java-7-oracle/include/linux/
	WIN_L_OPTS = -shared -lmingw32 -static-libstdc++ -static-libgcc
else ifeq ($(UNAME), CYGWIN_NT-6.1-WOW64)
	ext = so
	compiler = g++
	linker = g++
	JAVA_INCLUDE = -I/cygdrive/c/Program\ Files/Java/jdk1.7.0_25/include/ -I/cygdrive/c/Program\ Files/Java/jdk1.7.0_25/include/win32/
	CC_OPTS = -fPIC -std=c++11 -O3 
	L_OPTS = -static
	cp = cp
	mv = mv
	mkdir = mkdir -p
	rm = rm -f
	rm_dir = rm -rf
else # Windows
	rm = del
	ext = dll
	compiler = g++
	linker = g++
	JAVA_INCLUDE = -I"C:/Progra~1/Java/jdk1.7.0_25/include/" -I"C:/Progra~1/Java/jdk1.7.0_25/include/win32/"
	CC_OPTS = -std=c++11 -O3 -I"C:/MinGW/include/boost/math/"
	L_OPTS = -shared -lmingw32 -mconsole
	cp = xcopy
	mv = xcopy
	rm = del
	rm_dir = rmdir /s /q
	mkdir = mkdir
endif

StopPow: 
	$(SWIG) $(SWIG_OPTS) $(SWIG_FILE)
	$(compiler) $(JAVA_INCLUDE) $(CC_OPTS) -c $(SRC_FILES) $(SWIG_WRAP_SRC)
	$(compiler) $(L_OPTS) $(OBJ_FILES) -o $(OUT_FILE)
	$(mkdir) $(DEST_DIR)
	$(mv) *.java $(DEST_DIR)
	$(mv) $(OUT_FILE) $(DEST_DIR)
	$(rm) $(OBJ_FILES)
	$(rm) $(SWIG_WRAP_SRC)
	$(mkdir) $(JAR_TEMP_DIR)
	$(cp) -r $(DEST_DIR)/* $(JAR_TEMP_DIR)
	javac $(JAR_TEMP_DIR)/*.java
	$(rm) $(JAR_TEMP_DIR)/*.java
	jar cf $(JAR_FILE) $(JAR_TEMP_DIR)
	$(rm_dir) $(JAR_TEMP_DIR)

windows:
	$(SWIG) $(SWIG_OPTS) $(SWIG_FILE)
	$(wincomp) $(JAVA_INCLUDE) $(WIN_CC_OPTS) -c $(SRC_FILES) $(SWIG_WRAP_SRC)
	$(wincomp) $(WIN_L_OPTS) $(OBJ_FILES) -o $(OUT_FILE)
	$(mkdir) $(DEST_DIR)
	$(mv) *.java $(DEST_DIR)
	$(mv) $(OUT_FILE) $(DEST_DIR)
	$(rm) $(OBJ_FILES)
	$(rm) $(SWIG_WRAP_SRC)
	$(mkdir) $(JAR_TEMP_DIR)
	$(cp) -r $(DEST_DIR)/* $(JAR_TEMP_DIR)
	javac $(JAR_TEMP_DIR)/*.java
	$(rm) $(JAR_TEMP_DIR)/*.java
	jar cf $(JAR_FILE) $(JAR_TEMP_DIR)
	$(rm_dir) $(JAR_TEMP_DIR)

	
clean:
	$(rm_dir) $(DEST_DIR)
	$(rm) $(JAR_FILE)
	$(rm) *.o *.cxx *.java *.jnilib