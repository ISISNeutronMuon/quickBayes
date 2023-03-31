from quickBayes.functions.base import BaseFitFunction
from quickBayes.functions.BG import (LinearBG,
                                     FlatBG,
                                     NoBG)
from typing import List


def get_background_function(BG_type: str) -> (BaseFitFunction):
    """
    A simple function for getting the BG type
    :param BG_type: string of the BG type
    :return the BG function
    """
    if BG_type.lower() == 'linear':
        return LinearBG()
    elif BG_type.lower() == "flat":
        return FlatBG()
    elif BG_type.lower() == "none":
        return NoBG()
    else:
        raise ValueError("invalid BG function")


def update_guess(params: List[float], func: BaseFitFunction) -> List[float]:
    """
    Get an updated list of guesses, using the known params
    :param params: the known parameters
    :param func: the fitting function
    :return a list of updated guesses (keeps values from fit)
    """
    if len(params) > len(func.get_guess()):
        raise ValueError("Too many parameters")
    return params + func.get_guess()[len(params):]
