#!/usr/bin/env python

"""
setup.py file for StopPow library
"""

from distutils.core import setup, Extension


StopPow_module = Extension('_StopPow',
                           	sources=['StopPow_wrap.cxx', '../src/StopPow.cpp', '../src/StopPow_SRIM.cpp', '../src/StopPow_BetheBloch.cpp', '../src/StopPow_LP.cpp', '../src/PlotGen.cpp'],
							extra_compile_args = ['-std=c++11'],
							language="c++" )

setup (name = 'StopPow',
       version = '0.2',
       author      = "Alex Zylstra",
       description = """Stopping power library""",
       ext_modules = [StopPow_module],
       py_modules = ["StopPow","StopPow_SRIM","StopPow_BetheBloch","StopPow_LP","PlotGen"],
       )