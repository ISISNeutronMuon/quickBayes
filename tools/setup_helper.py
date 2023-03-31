from distutils.core import Extension
import numpy
from os.path import join
from typing import Sequence, Dict
from Cython.Build import cythonize
from pathlib import Path


def read_makefiles(package_name: str) -> Dict[str, str]:
    """
    Automatically locates all of the makefiles in src
    and reads the contents.
    :param package_name: the name of the package
    :return a dict of the python imports and file locations
    """
    makefile_data = []
    module_source_map = {}
    for path in Path('src').glob('**/*makefile.txt'):
        with open(path, 'r') as file:
            lines = file.readlines()
            for line in lines:
                makefile_data.append(line)

    for line in makefile_data:
        py_path, file_path = line.split(':')
        file_path = str(eval(file_path))
        module_source_map[f'{package_name}.{py_path}'] = [file_path]
    return module_source_map


def create_extension(fq_name: str,
                     sources: Sequence[str]) -> Extension:
    """
    Create an extension module
    :param fq_name: The final fully-qualified name of the module
    :param sources: List of relative paths from this file to the sources
    :return: An Extension class to be built
    """
    return Extension(name=fq_name,
                     sources=sources,
                     include_dirs=[numpy.get_include()])


def source_paths(dirname: str, filenames: Sequence[str]) -> Sequence[str]:
    """
    :param dirname: A relative path to the list of source files
    :param filenames: the names of the source files
    :return: A list of relative paths to the given sources in the directory
    """
    return [join(dirname, filename) for filename in filenames]


def get_extensions(package_name):
    module_source_map = read_makefiles(package_name)
    path = join('src', package_name)
    return cythonize([create_extension(name,
                      source_paths(str(path), sources)) for
                      name, sources in module_source_map.items()])
