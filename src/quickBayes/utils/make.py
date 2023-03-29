from os.path import join


"""
File to generate the list of
python files to be compiled
(using Cython) for
the package
"""


def get_utils(PACKAGE_NAME):
    module_source_map = {
        f'{PACKAGE_NAME}.utils.general':
            [join('utils', 'general.py')],
        f'{PACKAGE_NAME}.utils.spline':
            [join('utils', 'spline.py')],
        f'{PACKAGE_NAME}.utils.crop_data':
            [join('utils', 'crop_data.py')]
            }
    return module_source_map
