"""
Written a custom passer to write yml file
the default did not have nice formatting
and didn't allow for nice managment of
dependencies
"""


INDENT = " "


def record_name(yml, outfile):
    """
    Adds the name tag from the yml to the file
    :param yml: a dict containing the desired contents of the yml file
    :param outfile: the file the data will be written to
    """
    outfile.write(f"name: {yml['name']} \n \n")


def record_channels(yml, outfile):
    """
    Adds the channels tag from the yml to the file
    :param yml: a dict containing the desired contents of the yml file
    :param outfile: the file the data will be written to
    """
    outfile.write("channels: \n")
    outfile.write(f"{INDENT} - {yml['channels']} \n \n")


def record_pip(pip_dict, outfile):
    """
    Adds the packagesto be installed by pip to the file
    :param pip_dict: a dict containing the packages for pip
    :param outfile: the file the data will be written to
    """

    for package in pip_dict.keys():
        big_indent = INDENT + INDENT + INDENT
        outfile.write(f'{big_indent} - {package} {pip_dict[package]} \n')


def record_dependencies(yml, outfile):
    """
    Adds the dependencies to be listed in the file
    :param yml: a dict containing the desired contents of the yml file
    :param outfile: the file the data will be written to
    """

    outfile.write("dependencies: \n")

    deps = yml['dependencies']
    for package in deps.keys():
        if package == 'pip':
            outfile.write(f'{INDENT} - pip: \n')
            record_pip(deps['pip'], outfile)
        else:
            outfile.write(f'{INDENT} - {package} {deps[package]} \n')


def write_conda_yml_from_dict(yml, outfile):
    """
    Writes the file based on the input dict yml
    :param yml: a dict containing the desired contents of the yml file
    :param outfile: the file the data will be written to
    """

    record_name(yml, outfile)
    record_channels(yml, outfile)

    record_dependencies(yml, outfile)
