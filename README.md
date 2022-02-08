# quasielasticbayes

This package provides a convenient set of Python wrappers
for a set of routines used to perform Bayesian analysis
on quasi-elastic neutron-scattering data.

## Installing

The library is available on [PyPi](https://pypi.org/) and can be installed with pip:

```sh
> python -m pip install quasielasticbayes
```

## Development

The library code is primarily Fortran based.
[numpy.distutils.core.Extension](https://numpy.org/devdocs/reference/generated/numpy.distutils.core.Extension.html) is used to automate
the creation of Python-wrappers around each of the modules.

To build the package a fortran compiler such as `gfortran` is required
along with `numpy`.
A system package manager can be used but [mamba](https://github.com/conda-forge/miniforge#mambaforge) (a faster replacement for `conda`) is the preferred way to setup an enclosed development environment.

Download and install the latest version of `mamba` for your platform from the link above by running the script you download and following the instructions.

### Linux

Create a minimal conda environment named `quasielasticbayes-dev` (you can use another name if you prefer):

```sh
> mamba create --name quasielasticbayes-dev python=3.6 fortran-compiler numpy pytest
> conda activate quasielasticbayes-dev
> which python  # should produce something in $HOME/mambaforge
```

We currently use Python 3.6 to match the most common running environments for the
code: RHEL/Cent OS 7.

Now the conda environment is activated, compile and install the library in development mode:

```sh
> cd <path-to-repository-clone>  # this directory containing setup.py
> python -m pip install -v --editable .
```

You should now be able to run the tests:

```sh
> python -m pytest .
```

**Note:** There is currently a single test for `QLdata` and the output is based
on a version compiled with mingw 4.3 on Windows as this has been a version running
in production for years. This test fails for any other version of the compiler and
operating system.

To recompile the Fortran after editing it rerun:

```sh
python -m pip install -v --editable .
```

## Building for PyPi

If this is your first time interacting with PyPi then please see [here](https://packaging.python.org/en/latest/tutorials/packaging-projects/#uploading-the-distribution-archives) for instructions of how to setup accounts and API tokens.

Once built the wheel can be uploaded using twine:

```
twine upload ./dist/name_of_wheel
```

### Linux Notes

Linux wheels require a docker image to produce a `manylinux2010` wheel. For more details see this blog <https://uwekorn.com/2019/09/15/how-we-build-apache-arrows-manylinux-wheels.html>

### macOS Notes

Unfortunately we cannot avoid a dynamic link to libquadmath on OSX. This appears to be a licencing issue with GCC and the conda gfortran package doesn't include the static version of quadmath.
So in order to use this package you'll need to have a system install of gfortran. We're not the only ones with this issue, e.g. <https://github.com/aburrell/apexpy/issues/69>
