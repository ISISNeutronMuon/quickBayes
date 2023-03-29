import numpy
from distutils.core import Extension
from os.path import join
from typing import Sequence
from Cython.Build import cythonize
from src.quickBayes.fit_functions.make import get_fit_functions
from src.quickBayes.fit_engines.make import get_fit_engines
from src.quickBayes.test_helpers.make import get_test_helpers
from src.quickBayes.workflows.make import get_workflows
from src.quickBayes.utils.make import get_utils
from src.quickBayes.make import get_misc


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
    :return: A list of relative paths to the given sources in the directory
    """
    return [join(dirname, filename) for filename in filenames]


def get_v2_extensions(PACKAGE_NAME):
    module_source_map = {
                         **get_fit_functions(PACKAGE_NAME),
                         **get_fit_engines(PACKAGE_NAME),
                         **get_misc(PACKAGE_NAME),
                         **get_test_helpers(PACKAGE_NAME),
                         **get_utils(PACKAGE_NAME),
                         **get_workflows(PACKAGE_NAME),
                         }
    path = join('src', PACKAGE_NAME)
    return cythonize([create_extension(name,
                      source_paths(str(path), sources)) for
                      name, sources in module_source_map.items()])
