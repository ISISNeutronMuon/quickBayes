import argparse
import sys
from conda_dict_to_yml import write_conda_yml_from_dict


"""
f strings have not been used as the mac github action
does not recognise them and complains about
syntax errors.
"""


supported = ['windows', 'ubuntu', 'windows-latest', 'ubuntu-latest',
             'mac', 'macOS-latest']
exp = []

# Cannot move to 3.12 until gofit does
versions = ['3.8', '3.9', '3.10', '3.11']


def get_input():
    """
    get the OS and version from the command line arg
    :returns the OS and Python version
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('OS',
                        help='the OS for the yml file'
                        ' (windows, windows-latest, ubuntu, ubuntu-latest,'
                        ' mac, macOS-latest', type=str)
    parser.add_argument('version',
                        help='the Python version'
                        ' (3.8, 3.9, 3.10, 3.11)', type=str)
    args = parser.parse_args()

    if args.OS not in supported and args.OS not in exp:
        raise ValueError(str(args.OS) + ' is not a supported OS.')

    if args.version not in versions:
        raise ValueError(str(args.version) + ' is not a supported'
                                             ' Python version.')
    print(args.OS, args.version)
    return args.OS, args.version


def get_OS_info(OS, version):
    """
    Gets the yml_dict and file name from the OS
    :param OS: the OS the file is for
    :param version: the Python version
    :Return a dict of contents for the yml file and file name
    """
    default_yml = create_default(version)
    if OS == 'mac' or OS == 'macOS-latest':
        yml_dict = for_mac(default_yml)
    elif OS == 'windows' or OS == 'windows-latest':
        yml_dict = for_windows(default_yml)
    elif OS == 'ubuntu' or OS == 'ubuntu-latest':
        yml_dict = for_linux(default_yml)
    return yml_dict, str(yml_dict["name"])+'.yml'


"""
Here we list the conda recipies.
We start with a default, then write
methods to modify it for specific OS's.
"""


def create_default(version):
    """
    The default yml file for all OS's.
    Should only change these values if there is a
    good reason.
    :param version: Python version to build (exclude patch)
    """
    default_yml = {}

    pip_dict = {'readthedocs-sphinx-ext': ''}

    default_yml['name'] = 'quickBayes-dev'
    default_yml['channels'] = 'conda-forge'
    default_yml['dependencies'] = {'python': '=' + version + '.*',
                                   'numpy': '<2.0.0',
                                   'scipy': '',
                                   'pytest': '',
                                   'pre-commit': '>=2.15',
                                   'joblib': '',
                                   'Cython': '',
                                   'sphinx': '',
                                   'jupyter-book': '',
                                   'nbsphinx': '',
                                   '"pybind11[global]"': '',
                                   'eigen': '',
                                   'pip': pip_dict}
    return default_yml


def for_windows(yml_dict):
    """
    Edits the yml_dict to have Windows options
    :param yml_dict: the input yml_dict to edit
    :return the updated yml_dict
    """
    return yml_dict


def for_linux(yml_dict):
    """
    Edits the yml_dict to have ubuntu options
    :param yml_dict: the input yml_dict to edit
    :return the updated yml_dict
    """
    return yml_dict


def for_mac(yml_dict):
    """
    Edits the yml_dict to have mac options
    :param yml_dict: the input yml_dict to edit
    :return the updated yml_dict
    """
    return yml_dict


if __name__ == "__main__":
    try:
        OS, version = get_input()
        yml_dict, file_name = get_OS_info(OS, version)
        with open(file_name, 'w') as outfile:
            write_conda_yml_from_dict(yml_dict, outfile)

    except ValueError:
        error = sys.exc_info()[1]
        print(error)
