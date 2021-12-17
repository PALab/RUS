from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [
    Extension("rus_forward", ["rus_forward.pyx"]),
    Extension("rus_inverse", ["rus_inverse.pyx"]),
    Extension("rus_parser", ["rus_parser.pyx"]),
    Extension("rus_tools", ["rus_tools.pyx"]),
    Extension("rus", ["rus.pyx"]),
]



setup(
	  name = 'rus_proj',
	  cmdclass = {'build_ext': build_ext},
	  ext_modules = ext_modules
)