import numpy
from distutils.core import Extension
from os.path import join
from typing import Sequence
from Cython.Build import cythonize
from Cython.Compiler import Options


def create_extension(fq_name: str,
                     sources: Sequence[str]) -> Extension:
    """
    Create an extension module to be built by f2py
    :param fq_name: The final fully-qualified name of the module
    :param sources: List of relative paths from this file to the sources
    :return: An Extension class to be built
    """
    lang = 'C'
    # set to True for profiling
    Options.annotate = False
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
