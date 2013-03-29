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
	ext = so
else ifeq ($(UNAME), Linux)
	rm = rm
	ext = so
else # Windows
	rm = del
	ext = dll
endif

_StopPow: 
	$(SWIG) $(SWIG_OPTS) $(SWIG_FILE)
	$(PY) $(PY_SCRIPT) $(PY_OPTS)
	$(rm) -r build
	$(rm) *.cxx
	mkdir -p $(DEST_DIR)
	mv $(PY_DEST_FILE) $(DEST_DIR)
	mv *.$(ext) $(DEST_DIR)

clean:
	$(rm) -r $(DEST_DIR)