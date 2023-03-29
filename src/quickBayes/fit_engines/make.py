from os.path import join


"""
File to generate the list of
python files to be compiled
(using Cython) for
the package
"""


def get_fit_engines(PACKAGE_NAME):
    module_source_map = {
        f'{PACKAGE_NAME}.fitting.fit_utils':
            [join('fit_engines', 'fit_utils.py')],
        f'{PACKAGE_NAME}.fitting.fit_engine':
            [join('fit_engines', 'fit_engine.py')],
        f'{PACKAGE_NAME}.fitting.scipy_engine':
            [join('fit_engines', 'scipy_fit_engine.py')],
        f'{PACKAGE_NAME}.fitting.gofit_engine':
            [join('fit_engines', 'gofit_engine.py')]
            }
    return module_source_map
