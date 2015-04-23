from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize
from numpy import get_include

'''
compile by:
    python adaptive_moments_setup.py build_ext --inplace
'''
print('python adaptive_moments_setup.py build_ext --inplace')
## setup(
##     cmdclass = {'build_ext': build_ext},
##     ext_modules = [Extension("adaptive_moments", ["adaptive_moments.pyx"],
##         include_dirs=[get_include()])]
## )

extensions = [
    Extension('adaptive_moments', ['adaptive_moments.pyx'],
        include_dirs=[get_include()]),
]
setup(
        ext_modules = cythonize(extensions, annotate=True),
)

