from typing import Dict, List

"""
There are no direct tests for this class.
This is because all of the functionality
is tested by the classes that inherit
from it
"""


class BaseFitFunction(object):
    """
    A basic class outline for fit functions
    """
    def __init__(self, N_params: int, prefix: str):
        self._N_params = N_params
        self._prefix = prefix
        return

    def add_to_prefix(self, to_add: str):
        tmp = self._prefix.split('.')
        name = ''
        for j in range(len(tmp)-1):
            name += tmp[j] + '.'
        name += to_add + tmp[-1] + '.'
        self._prefix = name

    @property
    def N_params(self):
        return self._N_params

    def _add_to_report(self, name: str, value: float,
                       report_dict:
                       Dict[str, List[float]]) -> Dict[str, List[float]]:
        if name not in report_dict.keys():
            report_dict[name] = [value]
        else:
            report_dict[name].append(value)
        return report_dict

    def report(self):
        raise NotImplementedError()

    def __call__(self):
        raise NotImplementedError()
