# quickBayes

This package provides a convenient API
for performing Bayesian analysis
on muon and neutron-scattering data.

## Installing

The library is available on [PyPi](https://pypi.org/) and can be installed with pip:

```sh
> python -m pip install quickBayes
```

## Development

### Windows

To activate pre-commit, setup the Conda environment and then type ``mamba update pre-commit``.
Then once its updated type ``pre-commit install``, this generates a bash style file.
However, this will produce an error when you try to run it (a known issue with pre-commit).
In your terminal, but outside of your environment, type `` mamba env update --file quickBayes-dev-win.yml --prune``
Once it has completed, close the terminal and open a new one and activate the Conda environment.
In your terminal (Gitbash on Windows) you will now be able to run pre-commit checks (outside of your Conda environment).

### Linux

Create a minimal conda environment named `quickBayes-dev` (you can use another name if you prefer):

```sh
> mamba env create -f quickBayes-dev-linux.yml
```

### General

```sh
> conda activate quickBayes-dev
> which python  # should produce something in $HOME/mambaforge
```

Now the conda environment is activated, compile and install the library in development mode:

```sh
> cd <path-to-repository-clone>  # this directory containing setup.py
> python -m pip install -v --editable .
```

You should now be able to run the tests:

```sh
> python -m pytest .
```


## Building for PyPi

If this is your first time interacting with PyPi then please see [here](https://packaging.python.org/en/latest/tutorials/packaging-projects/#uploading-the-distribution-archives) for instructions of how to setup accounts and API tokens. The package will automatically be uploaded when a relase is created in github.

### Linux Notes

Linux wheels require a docker image to produce a `manylinux2010` wheel. For more details see this blog <https://uwekorn.com/2019/09/15/how-we-build-apache-arrows-manylinux-wheels.html>

### macOS Notes

Unfortunately we cannot avoid a dynamic link to libquadmath on OSX. This appears to be a licencing issue with GCC and the conda gfortran package doesn't include the static version of quadmath.
So in order to use this package you'll need to have a system install of gfortran. We're not the only ones with this issue, e.g. <https://github.com/aburrell/apexpy/issues/69>
