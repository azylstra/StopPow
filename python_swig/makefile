# name of binaries for python and swig
PY = python3
SWIG = swig

PY_SCRIPT = setup.py
PY_OPTS = build_ext --inplace
PY_DEST_FILE = StopPow.py

SWIG_FILE = StopPow.i
SWIG_OPTS = -python -c++

DEST_DIR = dist

UNAME = $(shell uname)
ifeq ($(UNAME), Darwin)
	rm = rm
	rmdir = rm -r
	mkdir = mkdir -p
	move = mv
	ext = so
else ifeq ($(UNAME), Linux)
	rm = rm
	rmdir = rm -r
	mkdir = mkdir -p
	move = mv
	ext = so
else # Windows
	rm = del
	rmdir = rmdir /s /q
	mkdir = mkdir
	move = move
	ext = pyd
	PY = python
	PY_OPTS = build_ext --inplace
endif

_StopPow: 
	$(SWIG) $(SWIG_OPTS) $(SWIG_FILE)
	$(PY) $(PY_SCRIPT) $(PY_OPTS)
	$(rmdir) build
	$(rm) *.cxx
	$(mkdir) $(DEST_DIR)
	$(move) $(PY_DEST_FILE) $(DEST_DIR)
	$(move) *.$(ext) $(DEST_DIR)

clean:
	$(rmdir) $(DEST_DIR)