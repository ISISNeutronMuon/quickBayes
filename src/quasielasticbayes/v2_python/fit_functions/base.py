from typing import Dict, List
from abc import ABC, abstractmethod
from numpy import ndarray


"""
There are no direct tests for this class.
This is because all of the functionality
is tested by the classes that inherit
from it
"""


class BaseFitFunction(ABC):
    """
    A basic class outline for fit functions
    """
    def __init__(self, N_params: int, prefix: str):
        """
        Base class for fit function
        :param N_params: number of parameters in function
        :param prefix: prefix for parameters when reporting
        """
        self._N_params = N_params
        self._prefix = prefix
        return

    def update_prefix(self, new: str) -> None:
        """
        Updates the begining of the prefix (before ":")
        Assume there are only 1 or 0 instances of ":"
        :param new: the new text to go before the ":"
        """
        if ":" not in self._prefix:
            self._prefix = new + self._prefix
        else:
            tmp = self._prefix.split(":")
            self._prefix = new + tmp[1]

    def add_to_prefix(self, to_add: str) -> None:
        """
        Used to update the prefix
        :param to_add: the name to insert
        """
        tmp = self._prefix.split('.')
        name = ''
        for j in range(len(tmp)-1):
            name += tmp[j] + '.'
        name += to_add + tmp[-1] + '.'
        self._prefix = name

    @property
    def N_params(self) -> int:
        """
        :return the number of parameters in function
        """
        return self._N_params

    @staticmethod
    def _add_to_report(name: str, value: float,
                       report_dict:
                       Dict[str, List[float]]) -> Dict[str, List[float]]:
        """
        Method for adding parameters to dict of results.
        If the param is present it will append the list
        :param name: name of the parameter
        :param value: the value for the parameter
        :param report_dict: the results dict
        :return the modified results dict
        """
        if name not in report_dict.keys():
            report_dict[name] = [value]
        else:
            report_dict[name].append(value)
        return report_dict

    def _read_report(self, report_dict: Dict[str, List[float]],
                     name: str, index: int) -> float:
        if name not in report_dict.keys():
            raise ValueError(f"parameter {name} not in results")
        tmp = report_dict[name]
        if index >= len(tmp):
            raise ValueError("Not enough parameters for this index")
        return tmp[index]

    @abstractmethod
    def read_from_report(self, report_dict: Dict[str, List[float]],
                         index: int = 0) -> List[float]:
        """
        Read the parameters from the results dict
        :param report_dict: the dict of results
        :param index: the index to get results from
        :return the parameters
        """
        raise NotImplementedError()

    @abstractmethod
    def report(self, results: Dict[str, List[float]],
               *kwargs: float) -> Dict[str, List[float]]:
        """
        This method is for accumalating the parameters
        into a single dict. This is useful for multiple
        calls where the aim is to get the results for
        different conditions (e.g. Q value).
        :param results: the dict of results
        :param kwargs: the parameters for the function
        :return the updated results dict
        """
        raise NotImplementedError()

    def report_errors(self, results: Dict[str, List[float]],
                      errors: List[float],
                      params: List[float]) -> Dict[str, List[float]]:
        """
        This method is for accumalating the parameter errors
        into a single dict. This is useful for multiple
        calls where the aim is to get the results for
        different conditions (e.g. Q value).
        :param results: the dict of parameter errors
        :param errors: the parameter errors
        :param params: the parameter values
        :return the updated results dict
        """
        return self.report(results, *errors)

    @abstractmethod
    def __call__(self, x: ndarray,
                 *kwargs: float) -> ndarray:
        """
        Implement the function call.
        Need to follow the expected
        form for scipy
        :param x: x values for function evaluation
        :param kwargs: parameters for the function
        :return y values for function evaluation
        """
        raise NotImplementedError()

    @abstractmethod
    def get_guess(self) -> List[float]:
        """
        Generates a guess for the fit values
        :return a list of guesses
        """
        raise NotImplementedError()

    @abstractmethod
    def get_bounds(self) -> (List[float], List[float]):
        """
        Gets the bounds for the fit
        :retun lists of the lower and upper bounds
        """
        raise NotImplementedError()
