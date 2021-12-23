"""Utility code to help with testing
"""
import base64
import json
import numpy as np


def _json_numpy_obj_hook(dct):
    """
    Decodes a previously encoded numpy ndarray
    with proper shape and dtype
    :param dct: (dict) json encoded ndarray
    :return: (ndarray) if input was an encoded ndarray
    """
    if isinstance(dct, dict) and '__ndarray__' in dct:
        data = base64.b64decode(dct['__ndarray__'])
        return np.frombuffer(data, dct['dtype']).reshape(dct['shape'])
    return dct


def load_json(*args, **kwargs):
    """Loads a json-encoded file-like object to a dictionary.
    Adds supports for decoding numpy arrays
    See json.load.
    """
    kwargs.setdefault('object_hook', _json_numpy_obj_hook)
    return json.load(*args, **kwargs)