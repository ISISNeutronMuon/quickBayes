import argparse
import sys
from conda_dict_to_yml import write_conda_yml_from_dict


supported = ['windows', 'ubuntu', 'windows-latest', 'ubuntu-latest',
             'mac', 'macOS-latest']
exp = ['mac', 'macOS-latest']


def get_OS():
    """
    get the OS from the command line arg
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('OS',
                        help='the OS for the yml file'
                        ' (windows, windows-latest, ubuntu, ubuntu-latest,'
                        ' mac, macOS-latest',
                        type=str)
    args = parser.parse_args()

    if args.OS not in supported and args.OS not in exp:
        raise ValueError(str(args.OS) + ' is not a supported OS.')
    print(args.OS)
    return args.OS


def get_OS_info(OS):
    """
    Gets the yml_dict and file name from the OS
    :param OS: the OS the file is for
    :Return a dict of contents for the yml file and file name
    """
    default_yml = create_default()
    if OS == 'mac' or OS == 'macOS-latest':
        yml_dict = for_mac(default_yml)
        file_name = f'{yml_dict["name"]}-mac.yml'


    elif OS == 'windows' or OS == 'windows-latest':
        yml_dict = for_windows(default_yml)
        file_name = f'{yml_dict["name"]}-win.yml'
    elif OS == 'ubuntu' or OS == 'ubuntu-latest':
        yml_dict = for_linux(default_yml)
        file_name = f'{yml_dict["name"]}-linux.yml'
    return yml_dict, file_name


"""
Here we list the conda recipies.
We start with a default, then write
methods to modify it for specific OS's.
"""


def create_default():
    """
    The default yml file for all OS's.
    Should only change these values if there is a
    good reason.
    """
    default_yml = {}

    # pip_dict = {"cython":
    # ">=0.29.32 # stops conda getting the wrong version",
    pip_dict = {'cython': '',
                'gofit': ''""}

    default_yml['name'] = 'quickBayes-dev'
    default_yml['channels'] = 'conda-forge'
    default_yml['dependencies'] = {'python': '=3.8.*',
                                   'numpy': '=1.16.*',
                                   'scipy': '',
                                   'pytest': '',
                                   'pre-commit': '>=2.15',
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
        OS = get_OS()
        yml_dict, file_name = get_OS_info(OS)
        with open(file_name, 'w') as outfile:
            write_conda_yml_from_dict(yml_dict, outfile)

    except ValueError:
        error = sys.exc_info()[1]
        print(error)
