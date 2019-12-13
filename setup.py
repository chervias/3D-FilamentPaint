from distutils.core import setup, Extension
from distutils.sysconfig import get_python_lib
import os

numpy_inc = os.path.join(get_python_lib(plat_specific=1), 'numpy/core/include')

module1 =  Extension('FilamentPaint',
                     sources = ['code/FilamentPaint.c','code/FilamentPaint_mod.c','code/query_polygon_wrapper.cpp'],
             		include_dirs = [numpy_inc,'code','/home/chervias/Software/miniconda3/envs/healpy/include/healpix_cxx/'],
                     libraries=['gsl','gslcblas','fftw3','healpix_cxx','cxxsupport','sharp','fftpack','c_utils','cfitsio'],
                     library_dirs = ["lib"],
                     extra_compile_args=['-fPIC','-Wall','-g'])

setup (name = 'FilamentPaint',
	   version = '0.1',
	   url='https://github.com/chervias/3dfilament-project-healpix',
	   description = '',
	   ext_modules = [module1],
#	   packages=['FilamentPaint' ],
#	   package_dir={'FilamentPaint':'code'},
#	   install_requires=['healpy', 'numpy','scipy','Bpowspec'],
#	   zip_safe=False
	  )
