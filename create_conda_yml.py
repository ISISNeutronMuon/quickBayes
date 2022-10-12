import argparse
import sys
from conda_dict_to_yml import write_conda_yml_from_dict


supported = ["w", "u", "windows-latest", "ubuntu-latest"]


def create_default():
    default_yml = {}

    pip_dict = {"cython": ">=0.29.32 # stops conda getting the wrong version"}

    default_yml['name'] = "quasielasticbayes-dev"
    default_yml['channels'] = 'conda-forge'
    default_yml['dependencies'] = {'python': '=3.8.*',
                                   'numpy': '=1.16.*',
                                   'scipy': '',
                                   'pytest': '',
                                   'pre-commit': '>=2.15',
                                   'pip': pip_dict}
    return default_yml


def for_windows(yml_dict):
    yml_dict['dependencies']['python'] = '=3.16.*'
    return yml_dict


def for_linux(yml_dict):
    yml_dict['dependencies']['gfortran'] = '>=2.2'
    return yml_dict


try:

    parser = argparse.ArgumentParser()
    parser.add_argument("OS",
                        help="the OS for the yml file"
                        " (w=windows-latest, u=ubuntu-latest",
                        type=str)
    args = parser.parse_args()

    if args.OS not in supported:
        raise ValueError("This OS is not supported. Please use a supported OS")
    print(args.OS)
    default_yml = create_default()
    if args.OS == 'w' or args.OS == "windows-latest":
        yml_dict = for_windows(default_yml)
        file_name = f'{yml_dict["name"]}-win.yml'
    elif args.OS == 'u' or args.OS == 'ubuntu-latest':
        yml_dict = for_linux(default_yml)
        file_name = f'{yml_dict["name"]}-linux.yml'

    with open(file_name, "w") as outfile:
        write_conda_yml_from_dict(yml_dict, outfile)

except ValueError:
    error = sys.exc_info()[1]

    print(error)
