"""
Setuptools support for building the Fortran extensions with
numpy.f2py
"""
# from os.path import join
# Importing setuptools modifies the behaviour of setup from distutils
# to support building wheels. It will be marked as unused by IDEs/static
# analysis.
# import setuptools
# import sys
# from typing import Sequence, Tuple

from numpy.distutils.core import setup

from fortran_setup import get_fortran_extensions, FortranExtensionBuilder
from cython_setup import get_cython_extensions
from numpy.distutils.command.build_src import build_src as _build_src


PACKAGE_NAME = 'quasielasticbayes'

extensions = (get_fortran_extensions(PACKAGE_NAME)
              + get_cython_extensions(PACKAGE_NAME))


class build_source(_build_src):

    def filter_files(self, sources, exts=[]):
        # cython causes an empty 'list' to be passed
        # so we check that the 'list' is not empty
        if sources == [[]]:
            return [], []
        return super().filter_files(sources, exts)


# from Cython.Distutils import build_ext as cython_build_ext
# import copy

class extension_builder(FortranExtensionBuilder):
    def initialize_options(self):
        super().initialize_options()
        # self.cython_dist = copy.deepcopy(self.distribution)
        # self.cython_build = cython_build_ext(self.cython_dist)
        # self.cython_build.initialize_options()

    def finialize_options(self):
        super().finalize_options()
        # self.cython_build.finialize_options()

    def run(self):
        super().run()

        # tmp = self.extensions
        # c_ext = []
        # for ext in tmp:
        #    name = ext.name
        #    if "stuff" in name:
        #        c_ext.append(ext)

        # self.cython_build.extensions = c_ext
        # self.cython_build.run()


setup(
    name=PACKAGE_NAME,
    install_requires=['numpy>=1.12'],
    packages=[PACKAGE_NAME],
    description='A Bayesian fitting package used for '
                'fitting quasi-elastic neutron scattering data.',
    long_description='This package wraps fortran Bayesian '
                     'fitting libraries using f2py. '
                     'An application of this package is '
                     'to fit quasi-elastic '
                     'neutron scattering data in Mantid '
                     '(https://www.mantidproject.org)',
    author='Mantid Team',
    ext_modules=extensions,
    author_email="mantid-help@mantidproject.org",
    url='https://www.mantidproject.org',
    version="0.1.1",
    license='BSD',
    package_dir={'': 'src'},  # allows setup to find py and f90 files
    cmdclass={'build_ext': extension_builder, 'build_src': build_source}
    # cmdclass={'build_ext': FortranExtensionBuilder}
)
