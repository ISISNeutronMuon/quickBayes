import argparse
import sys


VERSION_STRING = "version = "
VERSION_MARKS = '"'
ORDER = ['major', 'minor', 'patch', 'beta']
VERSION_SEP = ['.', '.', 'b', '']


def make_str(version):
    """
    Simple function to convert version dict into a string
    """
    txt = VERSION_MARKS
    for j in range(len(ORDER)):
        txt += version[ORDER[j]] + VERSION_SEP[j]
    txt += VERSION_MARKS + '\n'
    return txt


def update_version_str(txt: str, bump: str) -> str:
    """
    Extracts the version from text and then applies
    an appropriate bump
    """
    tmp = txt.split(VERSION_MARKS)
    tmp = tmp[1]
    version = {}

    version[ORDER[0]] = tmp.split(VERSION_SEP[0])[0]
    version[ORDER[1]] = tmp.split(VERSION_SEP[1])[1]
    version[ORDER[2]] = tmp.split(VERSION_SEP[2])[0].split(VERSION_SEP[1])[2]
    version[ORDER[3]] = tmp.split(VERSION_SEP[2])[1]

    version[bump] = str(int(version[bump]) + 1)
    index = ORDER.index(bump)
    # reset below current bump
    for j in range(index+1, len(ORDER)):
        version[ORDER[j]] = str(0)
    return VERSION_STRING + make_str(version)


def update_version(file_name: str, bump: str):
    """
    Reads the file and then updates the version
    according to the specified bump
    """
    if bump not in ORDER:
        raise ValueError("Specified bump is not allowed")

    lines = []
    with open(file_name, 'r') as file:
        lines = file.readlines()

    with open(file_name, 'w') as file:
        for txt in lines:
            if VERSION_STRING in txt:
                txt = update_version_str(txt, bump)
            file.write(txt)


def get_input():
    """
    get the bump for the version from the command line arg
    :returns the bump for the version
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('bump',
                        help='the bump to the setup file', type=str)
    args = parser.parse_args()
    print(args)
    return args.bump


if __name__ == "__main__":
    try:
        bump = get_input()
        update_version('pyproject.toml', bump)
    except ValueError:
        error = sys.exc_info()[1]
        print(error)
