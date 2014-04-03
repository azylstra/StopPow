
INCLUDE = -I../src
DIR = ../src/
INCLUDES_SRC = ../src/*.h

UNAME = $(shell uname)
ifeq ($(UNAME), Darwin)
	compiler = clang++
	opts = -c -Wall -std=c++11 -stdlib=libc++ -O3
	linker = libtool
	link_opts = -static -lgsl -o 
	LIB_FILE = StopPow.a
	cp = cp
	mv = mv
	mkdir = mkdir -p
	rm = rm -f
	rm_dir = rm -rf
	obj_ext = .o
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
	obj_ext = .o
else ifeq ($(UNAME), CYGWIN_NT-6.1-WOW64)
	compiler = g++
	opts = -c -Wall -fPIC -std=c++11 -O3
	linker = ar
	link_opts = rvs -o 
	LIB_FILE = StopPow.a
	cp = cp
	mv = mv
	mkdir = mkdir -p
	rm = rm -f
	rm_dir = rm -rf
	obj_ext = .o
else # assume Windows
	compiler = cl
	opts = /c /nologo /O2 /MD /W3 /GS- -IC:\gsl\x86\include -I"C:\Program Files (x86)\Microsoft Visual Studio 12.0\VC\include" /EHsc
	linker = link
	link_opts = /DLL /LIBPATH:"C:\Program Files (x86)\Windows Kits\8.1\Lib\winv6.3\um\x86" /LIBPATH:C:\gsl\x86\lib /LIBPATH:"C:\Program Files (x86)\Microsoft Visual Studio 12.0\VC\lib" /DEFAULTLIB:gsl.lib /DEFAULTLIB:cblas.lib /OUT:
	LIB_FILE = StopPow.dll
	mkdir = mkdir
	cp = copy
	mv = move
	rm = del
	rm_dir = rmdir /s /q
	INCLUDES_SRC = ..\src\*.h
	obj_ext = .obj
endif

objects = StopPow$(obj_ext) StopPow_Plasma$(obj_ext) StopPow_PartialIoniz$(obj_ext) StopPow_LP$(obj_ext) StopPow_BetheBloch$(obj_ext) StopPow_SRIM$(obj_ext) StopPow_AZ$(obj_ext) StopPow_Mehlhorn$(obj_ext) StopPow_Grabowski$(obj_ext) PlotGen$(obj_ext) AtomicData$(obj_ext) Spectrum$(obj_ext)

INCLUDES_DST = include
LIB_DST = lib

$(LIB_FILE): $(objects)
	$(linker) $(link_opts)$(LIB_FILE) $(objects)
	$(mkdir) $(INCLUDES_DST)
	$(cp) $(INCLUDES_SRC) $(INCLUDES_DST)
	$(mkdir) $(LIB_DST)
	$(rm) $(objects)
	$(mv) $(LIB_FILE) $(LIB_DST)$(slash)
	
StopPow$(obj_ext): $(DIR)StopPow.cpp 
	$(compiler) $(opts) $(INCLUDE) $(DIR)StopPow.cpp
	
StopPow_Plasma$(obj_ext): $(DIR)StopPow_Plasma.cpp 
	$(compiler) $(opts) $(INCLUDE) $(DIR)StopPow_Plasma.cpp
	
StopPow_PartialIoniz$(obj_ext): $(DIR)StopPow_PartialIoniz.cpp 
	$(compiler) $(opts) $(INCLUDE) $(DIR)StopPow_PartialIoniz.cpp

StopPow_LP$(obj_ext): $(DIR)StopPow.cpp $(DIR)StopPow_LP.cpp
	$(compiler) $(opts) $(INCLUDE) $(DIR)StopPow.cpp $(DIR)StopPow_LP.cpp

StopPow_BetheBloch$(obj_ext): $(DIR)StopPow.cpp $(DIR)StopPow_BetheBloch.cpp
	$(compiler) $(opts) $(INCLUDE) $(DIR)StopPow.cpp $(DIR)StopPow_BetheBloch.cpp

StopPow_SRIM$(obj_ext): $(DIR)StopPow.cpp $(DIR)StopPow_SRIM.cpp
	$(compiler) $(opts) $(INCLUDE) $(DIR)StopPow.cpp $(DIR)StopPow_SRIM.cpp

StopPow_AZ$(obj_ext): $(DIR)StopPow.cpp $(DIR)StopPow_AZ.cpp
	$(compiler) $(opts) $(INCLUDE) $(DIR)StopPow.cpp $(DIR)StopPow_AZ.cpp

StopPow_Mehlhorn$(obj_ext): $(DIR)StopPow.cpp $(DIR)StopPow_Mehlhorn.cpp
	$(compiler) $(opts) $(INCLUDE) $(DIR)StopPow.cpp $(DIR)StopPow_Mehlhorn.cpp

StopPow_Grabowski$(obj_ext): $(DIR)StopPow.cpp $(DIR)StopPow_Grabowski.cpp
	$(compiler) $(opts) $(INCLUDE) $(DIR)StopPow.cpp $(DIR)StopPow_Grabowski.cpp

PlotGen$(obj_ext): $(DIR)PlotGen.cpp
	$(compiler) $(opts) $(INCLUDE) $(DIR)PlotGen.cpp

AtomicData$(obj_ext): $(DIR)AtomicData.cpp
	$(compiler) $(opts) $(INCLUDE) $(DIR)AtomicData.cpp

Spectrum$(obj_ext): $(DIR)Spectrum.cpp
	$(compiler) $(opts) $(INCLUDE) $(DIR)Spectrum.cpp

clean:
	$(rm) *.o
	$(rm_dir) $(INCLUDES_DST)
	$(rm_dir) $(LIB_DST)