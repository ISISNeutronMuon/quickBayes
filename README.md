# quasielasticbayes

This package provides a convenient set of Python wrappers
for a set of routines used to perform Bayesian analysis
on quasi-elastic neutron-scattering data.

## Building the f2py python extensions

The simplest way to build the distribution is to use [conda](https://docs.conda.io/en/latest/miniconda.html) to create
a separate environment to build the distribution.

Create a minimal conda environment:

*Windows*

```

conda create --name fortran python=3.8
conda activate fortran
conda install numpy
conda install -c msys2 m2w64-gcc-fortran

```

*OSX and Linux*

```

conda create --name fortran python=3.8
conda activate fortran
conda install numpy
conda install -c conda-forge fortran-compiler
```

NOTE: If you're building on OSX with the conda compilers it likely you'll need to export the compiler flags export LDFLAGS="-undefined dynamic_lookup -bundle"

To build a wheel, run

```
python setup.py bdist_wheel
```

from the root directory of the repository. 



## Building for PyPi

If this is your first time interacting with PyPi then please see [here](https://packaging.python.org/en/latest/tutorials/packaging-projects/#uploading-the-distribution-archives) for instructions of how to setup accounts and API tokens. 

Once built the wheel can be uploaded using twine:

```
twine upload ./dist/name_of_wheel
```

### Linux Notes

Linux wheels require a docker image to produce a `manylinux2010` wheel. For more details see this blog https://uwekorn.com/2019/09/15/how-we-build-apache-arrows-manylinux-wheels.html

### macOS Notes

Unfortunately we cannot avoid a dynamic link to libquadmath on OSX. This appears to be a licencing issue with GCC and the conda gfortran package doesn't include the static version of quadmath.
So in order to use this package you'll need to have a system install of gfortran. We're not the only ones with this issue, e.g. https://github.com/aburrell/apexpy/issues/69 
