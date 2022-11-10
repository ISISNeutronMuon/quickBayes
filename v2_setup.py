import numpy
from distutils.core import Extension
from os.path import join
from typing import Sequence
from Cython.Build import cythonize


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
        f'{PACKAGE_NAME}.v2.functions.base':
            [join('fit_functions', 'base.py')],
        f'{PACKAGE_NAME}.v2.functions.BG':
            [join('fit_functions', 'BG.py')],
        f'{PACKAGE_NAME}.v2.functions.delta':
            [join('fit_functions', 'delta_function.py')],
        f'{PACKAGE_NAME}.v2.functions.lorentz':
            [join('fit_functions', 'lorentz.py')],
        f'{PACKAGE_NAME}.v2.functions.gaussian':
            [join('fit_functions', 'gaussian.py')],
        f'{PACKAGE_NAME}.v2.functions.composite':
            [join('fit_functions', 'composite_fun.py')],
        f'{PACKAGE_NAME}.v2.functions.convolution':
            [join('fit_functions', 'conv_with_res.py')],
        f'{PACKAGE_NAME}.v2.functions.qldata_function':
            [join('fit_functions', 'qldata_function.py')],
        f'{PACKAGE_NAME}.v2.functions.SE':
            [join('fit_functions', 'stretch_exp.py')],
        f'{PACKAGE_NAME}.v2.functions.qse_function':
            [join('fit_functions', 'qse.py')],

        f'{PACKAGE_NAME}.v2.fitting.scipy_fit':
            [join('fit_engines', 'scipy_fit.py')],

        f'{PACKAGE_NAME}.v2.log_likelihood':
            ['log_likelihood.py'],

        f'{PACKAGE_NAME}.v2.QlData':
            ['qldata_main.py'],

        f'{PACKAGE_NAME}.v2.utils.general':
            [join('utils', 'general.py')],
        f'{PACKAGE_NAME}.v2.utils.spline':
            [join('utils', 'spline.py')],
        f'{PACKAGE_NAME}.v2.utils.crop_data':
            [join('utils', 'crop_data.py')]

            }
    path = join('src', PACKAGE_NAME)
    return cythonize([create_extension(name,
                     source_paths(str(join(path, 'v2_python')), sources)) for
                     name, sources in module_source_map.items()])
