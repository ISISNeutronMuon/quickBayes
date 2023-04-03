# quickBayes

This package provides a convenient API
for performing Bayesian analysis
on muon and neutron-scattering data.

## Installing

The library is available on [PyPi](https://pypi.org/project/quickBayes/#description) and can be installed with pip:

```sh
> python -m pip install quickBayes
```

## Development

### Pre-commit (Windows)

To activate pre-commit, setup the Conda environment and then type ``mamba update pre-commit``.
Then once its updated type ``pre-commit install``, this generates a bash style file.
However, this will produce an error when you try to run it (a known issue with pre-commit).
In your terminal, but outside of your environment, type `` mamba env update --file quickBayes-dev-win.yml --prune``
Once it has completed, close the terminal and open a new one and activate the Conda environment.
In your terminal (Gitbash on Windows) you will now be able to run pre-commit checks (outside of your Conda environment).

### Setup


This project uses a minimal conda environment for development called `quickBayes-dev` (you can use another name if you prefer):

To create the environment run this command from the root directory of the project:

```sh
> mamba env create -f quickBayes-dev-<OS>.yml
```

where, ``<OS>`` is either `windows`, `linux` or `mac` depending on your operating system.
To activate the environment:

```sh
> conda activate quickBayes-dev
> which python  # should produce something in $HOME/mambaforge
```

To compile and install the library in development mode:

```sh
> cd <path-to-repository-clone>  # this directory containing setup.py
> python -m pip install -v --editable .
```

You should now be able to run the tests:

```sh
> python -m pytest .
```

## Building for PyPi

If this is your first time interacting with PyPi then please see [here](https://packaging.python.org/en/latest/tutorials/packaging-projects/#uploading-the-distribution-archives) for instructions of how to setup accounts and API tokens.
The package will automatically be uploaded when a release is created in Github, see [here](https://cibuildwheel.readthedocs.io/en/stable/options/) for the details.

### Linux Notes

Linux wheels require a docker image to produce a `manylinux2010` wheel. For more details see this blog <https://uwekorn.com/2019/09/15/how-we-build-apache-arrows-manylinux-wheels.html>
