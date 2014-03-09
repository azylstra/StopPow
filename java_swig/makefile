# name of binaries for swig
SWIG = swig

SWIG_FILE = StopPow.i
SWIG_OPTS = -java -package cStopPow -c++
SWIG_WRAP_SRC = StopPow_wrap.cxx

DEST_DIR = StopPow

SRC_FILES = ../src/StopPow.cpp ../src/StopPow_SRIM.cpp ../src/StopPow_LP.cpp ../src/StopPow_BetheBloch.cpp ../src/StopPow_AZ.cpp ../src/StopPow_Mehlhorn.cpp ../src/PlotGen.cpp ../src/AtomicData.cpp ../src/Spectrum.cpp
OBJ_FILES = StopPow_wrap.o StopPow.o StopPow_SRIM.o StopPow_LP.o StopPow_BetheBloch.o StopPow_AZ.o StopPow_Mehlhorn.o PlotGen.o AtomicData.o Spectrum.o

OUT_FILE = libcStopPow.jnilib

JAR_TEMP_DIR = cStopPow
JAR_FILE = libcStopPow.jar

# system dependent options:
UNAME = $(shell uname)
ifeq ($(UNAME), Darwin)
	ext = so
	compiler = clang++
	linker = clang++
	JAVA_INCLUDE = -I/System/Library/Frameworks/JavaVM.framework/Versions/Current/Headers/
	CC_OPTS = -c -Wall -std=c++11 -stdlib=libc++ -O3
	L_OPTS = -std=c++11 -stdlib=libc++ -framework JavaVM -bundle -o 
	cp = cp
	mv = mv
	mkdir = mkdir -p
	rm = rm -f
	rm_dir = rm -rf
	slashstar = /*
else ifeq ($(UNAME), Linux)
	ext = so
	compiler = g++
	linker = g++
	JAVA_INCLUDE = -I/usr/lib/jvm/java-7-oracle/include/ -I/usr/lib/jvm/java-7-oracle/include/linux/
	CC_OPTS = -fPIC -std=c++11 -O3
	L_OPTS = -shared -o 
	cp = cp
	mv = mv
	mkdir = mkdir -p
	rm = rm -f
	rm_dir = rm -rf
	wincomp = i686-w64-mingw32-g++
	WIN_CC_OPTS = -std=c++0x -O3 -D_JNI_IMPLEMENTATION_ -m32 -mwindows -I/usr/lib/jvm/java-7-oracle/include/ -I/usr/lib/jvm/java-7-oracle/include/linux/
	WIN_L_OPTS = -shared -lmingw32 -static-libstdc++ -static-libgcc
	slashstar = /*
else # Windows
	rm = del
	ext = dll
	compiler = cl
	linker = link
	JAVA_INCLUDE = -I"C:\Program Files (x86)\Java\jdk1.7.0_45\include" -I"C:\Program Files (x86)\Java\jdk1.7.0_45\include\win32" -I"C:\Program Files (x86)\Microsoft Visual Studio 12.0\VC\include" -I"C:\boost\include\boost-1_55"
	CC_OPTS = /O2 /EHsc
	OBJ_FILES = StopPow_wrap.obj StopPow.obj StopPow_SRIM.obj StopPow_LP.obj StopPow_BetheBloch.obj StopPow_AZ.obj StopPow_Mehlhorn.obj PlotGen.obj AtomicData.obj
	L_OPTS = /DLL /LIBPATH:C:\boost\lib /LIBPATH:"C:\Program Files (x86)\Microsoft Visual Studio 12.0\VC\lib" /LIBPATH:"C:\Program Files (x86)\Windows Kits\8.1\Lib\winv6.3\um\x86" /OUT:
	cp = copy
	mv = move
	rm = del
	rm_dir = rmdir /s /q
	mkdir = mkdir
	slashstar = \*
endif

StopPow: 
	$(SWIG) $(SWIG_OPTS) $(SWIG_FILE)
	$(compiler) $(JAVA_INCLUDE) $(CC_OPTS) -c $(SRC_FILES) $(SWIG_WRAP_SRC)
	$(linker) $(L_OPTS)$(OUT_FILE) $(OBJ_FILES)
	$(mkdir) $(DEST_DIR)
	$(cp) *.java $(DEST_DIR)
	$(cp) $(OUT_FILE) $(DEST_DIR)
	$(rm) $(OBJ_FILES)
	$(rm) $(SWIG_WRAP_SRC)
	$(mkdir) $(JAR_TEMP_DIR)
	$(mv) $(OUT_FILE) $(JAR_TEMP_DIR)
	javac *.java
	$(rm) *.java
	$(mv) *.class $(JAR_TEMP_DIR)
	jar cf $(JAR_FILE) $(JAR_TEMP_DIR)
	$(rm_dir) $(JAR_TEMP_DIR)
	
clean:
	$(rm_dir) $(DEST_DIR)
	$(rm) $(JAR_FILE)
	$(rm) *.o *.cxx *.java *.jnilib