from distutils.core import setup, Extension
from distutils.sysconfig import get_python_lib
import os

numpy_inc = os.path.join(get_python_lib(plat_specific=1), 'numpy/core/include')

module1 =  Extension('FilamentPaint',
                     sources = ['code/FilamentPaint.c','code/FilamentPaint_mod.c'],
                     include_dirs = [numpy_inc,'code'],
                     libraries=['gsl'],
                     library_dirs = ["lib"],
                     extra_compile_args=['-fPIC','-Wall','-g'])

setup (name = 'FilamentPaint',
       version = '0.1',
       url='https://github.com/chervias/3dfilament-project-healpix',
       description = '',
       ext_modules = [module1],
       )
