"""
File to generate the list of
python files to be compiled
(using Cython) for
the package
"""


def get_misc(PACKAGE_NAME):
    module_source_map = {
        f'{PACKAGE_NAME}.log_likelihood':
            ['log_likelihood.py']
            }
    return module_source_map
