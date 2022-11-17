from quasielasticbayes.v2.functions.base import BaseFitFunction
from numpy import ndarray
import numpy as np
from typing import Dict, List


class CompositeFunction(BaseFitFunction):
    def __init__(self, prefix: str = ''):
        """
        Defines a function to wrap a sum of functions
        :param prefix: the prefix for parameters
        """
        super().__init__(0, prefix)
        self._funcs = []

    def add_function(self, func: BaseFitFunction) -> None:
        """
        Adds a function to the sum
        :param func: the function to add
        """
        func.add_to_prefix(f'f{len(self._funcs)+1}')
        self._funcs.append(func)
        self._N_params += func.N_params

    def split_args(self, args: List[float]) -> List[List[float]]:
        """
        Split the single args list into a list of lists
        for use with individual functions
        :param args: list of all the parameters
        :return list of the parameters for each individual function
        """
        j = 0
        split = []
        for func in self._funcs:
            N = func.N_params
            split.append(list(args[j:j+N]))
            j += N
        return split

    def __call__(self, x: ndarray, *args: float) -> ndarray:
        """
        Implement a sum of functions.
        Need to follow the expected
        form for scipy
        :param x: x values for function evaluation
        :param args: parameters for functions
        :return y values for evaluated function
        """
        if len(self._funcs) == 0:
            return np.zeros(len(x))
        elif len(args) != self.N_params:
            raise ValueError(f"Expected {self.N_params} args, got {len(args)}")

        fun_args = self.split_args(args)
        result = np.zeros(len(x))
        for j, func in enumerate(self._funcs):
            result += func(x, *fun_args[j])
        return result

    def read_from_report(self, report_dict: Dict[str, List[float]],
                         index: int = 0) -> List[float]:
        """
        Read the parameters from the results dict
        :param report_dict: the dict of results
        :param index: the index to get results from
        :return the parameters
        """
        params = []
        for fun in self._funcs:
            params += fun.read_from_report(report_dict, index)
        return params

    def report(self, report_dict: Dict[str, List[float]],
               *args: float) -> Dict[str, List[float]]:
        """
        report the results
        :param report_dic: dict of results
        :param args: parameters for functions
        :return updated dict of results
        """
        if len(args) != self.N_params:
            raise ValueError(f"Expected {self.N_params} args, got {len(args)}")

        fun_args = self.split_args(args)
        for j, func in enumerate(self._funcs):
            report_dict = func.report(report_dict, *fun_args[j])
        return report_dict

    def get_guess(self) -> List[float]:
        """
        Get the starting guess for a fit function
        :return the initial guess
        """
        guess = []
        for func in self._funcs:
            guess += func.get_guess()
        return guess

    def get_bounds(self) -> (List[float], List[float]):
        """
        Get the fitting bounds
        :return lists for lower and upper bounds
        """
        upper, lower = [], []
        for func in self._funcs:
            bounds = func.get_bounds()
            lower += bounds[0]
            upper += bounds[1]
        return lower, upper
