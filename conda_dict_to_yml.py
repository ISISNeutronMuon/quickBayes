"""
Written a custom passer to write yml file
the default did not have nice formatting
and didn't allow for nice managment of
dependencies
"""


INDENT = " "


def record_name(yml, outfile):
    outfile.write(f"name: {yml['name']} \n \n")


def record_channels(yml, outfile):
    outfile.write("channels: \n")
    outfile.write(f"{INDENT} - {yml['channels']} \n \n")


def record_pip(pip_dict, outfile):
    for package in pip_dict.keys():
        big_indent = INDENT + INDENT + INDENT
        outfile.write(f'{big_indent} - {package} {pip_dict[package]} \n')


def record_dependencies(yml, outfile):
    outfile.write("dependencies: \n")

    deps = yml['dependencies']
    for package in deps.keys():
        if package == 'pip':
            outfile.write(f'{INDENT} - pip: \n')
            record_pip(deps['pip'], outfile)
        else:
            outfile.write(f'{INDENT} - {package} {deps[package]} \n')


def write_conda_yml_from_dict(yml, outfile):
    record_name(yml, outfile)
    record_channels(yml, outfile)

    record_dependencies(yml, outfile)
