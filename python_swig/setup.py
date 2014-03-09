#!/usr/bin/env python

"""
setup.py file for StopPow library
"""

# detect the OS type:
import platform

# options for *nix:
cargs = ['-O3']
if platform.system() == 'Linux':
    cargs = ['-O3','-fPIC','-std=c++11']
    largs = []
elif platform.system() == 'Darwin':
    cargs = ['-stdlib=libc++','-std=c++11','-O3']
    largs = []
elif platform.system() == 'Windows':
    cargs = ['/O2',r'-IC:\boost\include\boost-1_55', '/EHsc']
    largs = [r'/LIBPATH:C:\boost\lib','libboost_exception-vc120-mt-gd-1_55.lib']

from distutils.core import setup, Extension


StopPow_module = Extension('_StopPow',
                            sources=['StopPow_wrap.cxx', '../src/StopPow.cpp', '../src/StopPow_SRIM.cpp', '../src/StopPow_BetheBloch.cpp', '../src/StopPow_LP.cpp', '../src/StopPow_AZ.cpp','../src/StopPow_Mehlhorn.cpp','../src/AtomicData.cpp','../src/PlotGen.cpp'],
                            extra_compile_args = cargs,
                            extra_link_args = largs,
                            language="c++" )

setup (name = 'StopPow',
       version = '0.2',
       author      = "Alex Zylstra",
       description = """Stopping power library""",
       ext_modules = [StopPow_module],
       py_modules = ["StopPow","StopPow_SRIM","StopPow_BetheBloch","StopPow_LP","StopPow_AZ","StopPow_Mehlhorn","AtomicData","PlotGen"],
       )