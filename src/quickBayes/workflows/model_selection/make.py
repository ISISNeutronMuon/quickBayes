from os.path import join


"""
File to generate the list of
python files to be compiled
(using Cython) for
the package
"""


def get_fit_functions(PACKAGE_NAME):
    module_source_map = {
        f'{PACKAGE_NAME}.workflow.model_template':
            [join('workflows', 'model_selection',
                  'model_template.py')],
        f'{PACKAGE_NAME}.workflow.QlData':
            [join('workflows', 'model_selection',
                  'qldata_main.py')],
        f'{PACKAGE_NAME}.workflow.QSE':
            [join('workflows', 'model_selection',
                  'qse_main.py')],
        f'{PACKAGE_NAME}.workflow.MuonExpDecay':
            [join('workflows', 'model_selection',
                  'muon_exp_decay_main.py')]
            }
    return module_source_map
