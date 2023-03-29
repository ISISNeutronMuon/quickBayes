from os.path import join


"""
File to generate the list of
python files to be compiled
(using Cython) for
the package
"""


def get_fit_functions(PACKAGE_NAME):
    module_source_map = {
        f'{PACKAGE_NAME}.functions.base':
            [join('fit_functions', 'base.py')],
        f'{PACKAGE_NAME}.functions.BG':
            [join('fit_functions', 'BG.py')],
        f'{PACKAGE_NAME}.functions.delta':
            [join('fit_functions', 'delta_function.py')],
        f'{PACKAGE_NAME}.functions.lorentz':
            [join('fit_functions', 'lorentz.py')],
        f'{PACKAGE_NAME}.functions.gaussian':
            [join('fit_functions', 'gaussian.py')],
        f'{PACKAGE_NAME}.functions.composite':
            [join('fit_functions', 'composite_fun.py')],
        f'{PACKAGE_NAME}.functions.convolution':
            [join('fit_functions', 'conv_with_res.py')],
        f'{PACKAGE_NAME}.functions.qldata_function':
            [join('fit_functions', 'qldata_function.py')],
        f'{PACKAGE_NAME}.functions.qe_function':
            [join('fit_functions', 'quasielastic_function.py')],
        f'{PACKAGE_NAME}.functions.SE_fix':
            [join('fit_functions', 'stretch_exp_fixed.py')],
        f'{PACKAGE_NAME}.functions.SE':
            [join('fit_functions', 'stretch_exp.py')],
        f'{PACKAGE_NAME}.functions.qse_function':
            [join('fit_functions', 'qse.py')],
        f'{PACKAGE_NAME}.functions.exp_decay':
            [join('fit_functions', 'exp_decay.py')],
        f'{PACKAGE_NAME}.functions.qse_fixed':
            [join('fit_functions', 'qse_fixed.py')]
            }
    return module_source_map
