from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize
from Cython.Compiler import Options
import numpy

Options.annotate = True
lang = 'C'


############################################################
# fortran python
############################################################

ext_modules = [Extension('fortran_python',
                         sources=['fortran_python.pyx'],include_dirs=[numpy.get_include()],
                         language=lang,annotate=True
                        )]

# use these to get annotation -> need to work out how to use it with numpy
#ext_modules=cythonize(
#        ['fortran_python.pyx'],
#        annotate=True)
setup(
name = 'fortran_python',
cmdclass = {'build_ext': build_ext},
ext_modules = ext_modules
)

############################################################
# four
############################################################

ext_modules = [Extension('four',
                         sources=['four.pyx'],include_dirs=[numpy.get_include()],
                         language=lang,annotate=True
                        )]

setup(
name = 'four',
cmdclass = {'build_ext': build_ext},
ext_modules = ext_modules
)

############################################################
# bayes_C
############################################################

ext_modules = [Extension('bayes_C',
                         sources=['bayes_C.pyx'],include_dirs=[numpy.get_include()],
                         language=lang,annotate=True
                        )]

setup(
name = 'bayes_C',
cmdclass = {'build_ext': build_ext},
ext_modules = ext_modules
)
############################################################
# qldata_subs
############################################################

#ext_modules = [Extension('qldata_subs',
#                         sources=['qldata_subs.pyx'],include_dirs=[numpy.get_include()],
#                         language=lang,annotate=True
#                        )]

#setup(
#name = 'qldata_subs',
#cmdclass = {'build_ext': build_ext},
#ext_modules = ext_modules
#)
# python setup.py build_ext -i clean