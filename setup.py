"""
Setuptools support for building the Fortran extensions with
numpy.f2py
"""
from pathlib import PurePosixPath
# Importing setuptools modifies the behaviour of setup from distutils
# to support building wheels. It will be marked as unused by
# IDEs/static analysis.
import sys
from typing import Sequence, Tuple
import setuptools # noqa F401
from numpy.distutils.core import (Extension as FortranExtension, setup)
from numpy.distutils.command.build_ext import build_ext as _build_ext

PACKAGE_NAME = 'quasielasticbayes'


def create_fortran_extension(fq_name: str,
                             sources: Sequence[str]) -> FortranExtension:
    """
    Create an extension module to be built by f2py
    :param fq_name: The final fully-qualified name of the module
    :param sources: List of relative paths from this file to the sources
    :return: An Extension class to be built
    """
    (extra_compile_args, extra_f90_compile_args,
     extra_link_args) = compiler_flags()
    return FortranExtension(name=fq_name,
                            sources=sources,
                            extra_f90_compile_args=extra_f90_compile_args,
                            extra_link_args=extra_link_args,
                            extra_compile_args=extra_compile_args)


def compiler_flags() -> Tuple[Sequence[str], Sequence[str], Sequence[str]]:
    """
    :return: The compiler flags appropriate for this platform
    """
    if sys.platform == 'win32':
        # Compile static libs on windows
        extra_compile_args = []
        extra_f90_compile_args = ["-O0"]
        extra_link_args = ["-static-libgfortran", "-static-libgcc"]
    elif sys.platform == 'darwin':
        extra_compile_args = ['-Wno-argument-mismatch']
        extra_f90_compile_args = ['-ff2c']
        extra_link_args = ['-dynamiclib', '-lgfortran',
                           '-Bstatic', '-static-libgfortran',
                           "-static-libgcc"]
        # extra_link_args = ["-static", "-static-libgcc"]
    else:
        # On Linux we build a manylinux2010
        # (https://www.python.org/dev/peps/pep-0571/#the-manylinux2010-policy)
        # wheel that assumes compatible versions of bases
        # libraries are installed.
        extra_compile_args = []
        extra_f90_compile_args = ["-fallow-argument-mismatch"]
        extra_link_args = []

    return extra_compile_args, extra_f90_compile_args, extra_link_args


def source_paths(dirname: PurePosixPath,
                 filenames: Sequence[str]) -> Sequence[str]:
    """
    :param dirname: A relative path to the list of source files
    :return: A list of relative paths to the given sources in the directory
    """
    return [str(dirname / filename) for filename in filenames]


class FortranExtensionBuilder(_build_ext):
    """Custom extension builder to use specific compilers"""

    def finalize_options(self):
        _build_ext.finalize_options(self)
        # If we don't do this on windows, when we do bdist_wheel we wont
        # get a static link
        # this is because it misses the compiler flags to f2py which means
        # it ignores the static flags we try to pass
        if sys.platform == 'win32':
            self.fcompiler = 'gnu95'
            self.compiler = 'mingw32'


# Start setup
# Create extension builders
module_source_map = {
    f'{PACKAGE_NAME}.ResNorm': ['ResNorm_main.f90',
                                'ResNorm_subs.f90',
                                'BlrRes.f90',
                                'Bayes.f90',
                                'Four.f90',
                                'Util.f90'],
    f'{PACKAGE_NAME}.Quest': ['Quest_main.f90',
                              'Quest_subs.f90',
                              'BlrRes.f90',
                              'Bayes.f90',
                              'Four.f90',
                              'Util.f90',
                              'Simopt.f90'],
    f'{PACKAGE_NAME}.QLse':
        ['QLse_main.f90',
         'QLse_subs.f90',
         'BlrRes.f90',
         'Bayes.f90',
         'Four.f90',
         'Util.f90',
         'Simopt.f90'],
    f'{PACKAGE_NAME}.QLres':
        ['QLres_main.f90',
         'QLres_subs.f90',
         'BlrRes.f90',
         'Bayes.f90',
         'Four.f90',
         'Util.f90'],
    f'{PACKAGE_NAME}.QLdata': ['QLdata_main.f90',
                               'QLdata_subs.f90',
                               'Bayes.f90',
                               'Four.f90',
                               'Util.f90'],
    f'{PACKAGE_NAME}.Four': ['Four_main.f90',
                             'Four.f90']
}
extensions = [create_fortran_extension(name,
              source_paths(PurePosixPath('quasielasticbayes'), sources)) for
              name, sources in module_source_map.items()]

setup(
    name=PACKAGE_NAME,
    install_requires=['numpy>=1.12'],
    packages=[PACKAGE_NAME],
    description='A Bayesian fitting package used for fitting '
                'quasi-elastic neutron scattering data.',
    long_description='This package wraps fortran Bayesian '
                     'fitting libraries using f2py. '
                     'An application of this package is to '
                     'fit quasi-elastic'
                     'neutron scattering data in Mantid '
                     '(https://www.mantidproject.org)',
    author='Mantid Team',
    ext_modules=extensions,
    author_email="mantid-help@mantidproject.org",
    url='https://www.mantidproject.org',
    version="0.1.1",
    license='BSD',
    cmdclass={'build_ext': FortranExtensionBuilder}
)
