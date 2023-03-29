from os.path import join


"""
File to generate the list of
python files to be compiled
(using Cython) for
the package
"""


def get_grid_search(PACKAGE_NAME):
    module_source_map = {
        f'{PACKAGE_NAME}.workflow.grid_template':
            [join('workflows', 'grid_search',
                  'grid_search_template.py')],
        f'{PACKAGE_NAME}.workflow.qse_search':
            [join('workflows', 'grid_search',
                  'qse_grid_search.py')],
            }
    return module_source_map
