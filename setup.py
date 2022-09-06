"""
Setuptools support for building the Fortran extensions with
numpy.f2py and cython extensions
"""
import os
from setuptools.command.build_py import build_py as _build_py
import sysconfig
from numpy.distutils.core import setup
from numpy.distutils.command.build_src import build_src as _build_src
from cython_setup import get_cython_extensions
from fortran_setup import get_fortran_extensions, FortranExtensionBuilder


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


# compile the fortran code and the build py and src classes
# will ensure that the cython is handled correctly
class extension_builder(FortranExtensionBuilder):
    def initialize_options(self):
        super().initialize_options()


# noinspection PyPep8Naming
class build_py(_build_py):
    # ignore py files if there are compiled extensions with the same name
    def find_package_modules(self, package, package_dir):
        ext_suffix = sysconfig.get_config_var('EXT_SUFFIX')
        modules = super().find_package_modules(package, package_dir)
        filtered_modules = []
        for (pkg, mod, filepath) in modules:
            if os.path.exists(filepath.replace('.py', ext_suffix)):
                continue
            filtered_modules.append((pkg, mod, filepath, ))
        return filtered_modules


setup(
    name=PACKAGE_NAME,
    install_requires=['cython', 'numpy>=1.12'],
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
    cmdclass={'build_ext': extension_builder,
              'build_src': build_source,
              'build_py': build_py}
)
