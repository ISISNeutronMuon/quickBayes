"""
Written a custom passer to write yml file
the default did not have nice formatting
and didn't allow for nice managment of
dependencies

f strings have not been used as the mac github action
does not recognise them and complains about
syntax errors.
"""


INDENT = ' '


def record_name(yml, outfile):
    """
    Adds the name tag from the yml to the file
    :param yml: a dict containing the desired contents of the yml file
    :param outfile: the file the data will be written to
    """
    outfile.write('name: ' + str(yml['name']) + ' \n \n')


def record_channels(yml, outfile):
    """
    Adds the channels tag from the yml to the file
    :param yml: a dict containing the desired contents of the yml file
    :param outfile: the file the data will be written to
    """
    outfile.write('channels: \n')
    outfile.write(str(INDENT) + ' - ' + str(yml['channels']) + ' \n \n')


def record_pip(pip_dict, outfile):
    """
    Adds the packagesto be installed by pip to the file
    :param pip_dict: a dict containing the packages for pip
    :param outfile: the file the data will be written to
    """

    for package in pip_dict.keys():
        big_indent = INDENT + INDENT + INDENT
        outfile.write(str(big_indent) + ' - ' + str(package) +
                      ' ' + str(pip_dict[package]) + ' \n')


def record_dependencies(yml, outfile):
    """
    Adds the dependencies to be listed in the file
    :param yml: a dict containing the desired contents of the yml file
    :param outfile: the file the data will be written to
    """

    outfile.write('dependencies: \n')

    deps = yml['dependencies']
    for package in deps.keys():
        if package == 'pip' and deps['pip'] != {}:
            outfile.write(str(INDENT) + ' - pip: \n')
            record_pip(deps['pip'], outfile)
        elif package != 'pip':
            outfile.write(str(INDENT) + ' - ' + str(package) +
                          ' ' + str(deps[package]) + ' \n')


def write_conda_yml_from_dict(yml, outfile):
    """
    Writes the file based on the input dict yml
    :param yml: a dict containing the desired contents of the yml file
    :param outfile: the file the data will be written to
    """

    record_name(yml, outfile)
    record_channels(yml, outfile)

    record_dependencies(yml, outfile)
