# from numpy.distutils.command.build_ext import build_ext as _build_ext
# from numpy.distutils.command.build_src import build_src as _build_src
# from distutils.command.build import build
from distutils.core import Extension
from Cython.Compiler import Options
from os.path import join
# Importing setuptools modifies the behaviour of setup from distutils
# to support building wheels. It will be marked as unused by IDEs/static
# analysis.
# import setuptools
# import sys
from typing import Sequence
from Cython.Build import cythonize
# import copy
import numpy


def create_extension(fq_name: str,
                             sources: Sequence[str]) -> Extension:
    """
    Create an extension module to be built by f2py
    :param fq_name: The final fully-qualified name of the module
    :param sources: List of relative paths from this file to the sources
    :return: An Extension class to be built
    """
    lang = 'C'
    Options.annotate = True
    return Extension(name=fq_name,
                     sources=sources,
                     include_dirs=[numpy.get_include()],
                     language=lang)


def source_paths(dirname: str, filenames: Sequence[str]) -> Sequence[str]:
    """
    :param dirname: A relative path to the list of source files
    :return: A list of relative paths to the given sources in the directory
    """
    return [join(dirname, filename) for filename in filenames]


def get_cython_extensions(PACKAGE_NAME):

    module_source_map = {
        f'{PACKAGE_NAME}.stuff':
            ['stuff.pyx']}
    return cythonize([create_extension(name,
                      source_paths(str(join('src', 'c_python')), sources)) for
                      name, sources in module_source_map.items()])
