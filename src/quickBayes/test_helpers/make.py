from os.path import join


"""
File to generate the list of
python files to be compiled
(using Cython) for
the package
"""


def get_test_helpers(PACKAGE_NAME):
    module_source_map = {
        f'{PACKAGE_NAME}.test_helpers.template_fit_test':
            [join('test_helpers', 'template_test_fit.py')],
        f'{PACKAGE_NAME}.test_helpers.fitting_data':
            [join('test_helpers', 'fitting_data.py')],
        f'{PACKAGE_NAME}.test_helpers.template_scipy_fit':
            [join('test_helpers', 'template_scipy_fit_test.py')],
        f'{PACKAGE_NAME}.test_helpers.workflows':
            [join('test_helpers', 'workflow_helper.py')],
           }
    return module_source_map
