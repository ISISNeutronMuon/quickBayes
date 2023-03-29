from os.path import join
from src.quickBayes.workflows.grid_search.make import get_grid_search


"""
File to generate the list of
python files to be compiled
(using Cython) for
the package
"""


def get_workflows(PACKAGE_NAME):
    module_source_map = {
        f'{PACKAGE_NAME}.workflow.template':
            [join('workflows', 'workflow_template.py')],
            }
    return {**module_source_map,
            **get_grid_search(PACKAGE_NAME)
            }
