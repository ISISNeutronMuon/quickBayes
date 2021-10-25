# Building the f2py python extensions

You can build this package using the setup.py file. To use this, create a minimal conda environment as shown below:
```

conda create --name fortran python=3.8
conda activate fortran
conda install numpy
conda install -c msys2 m2w64-gcc-fortran

```
And on OSX and Linux

```

conda create --name fortran python=3.8
conda activate fortran
conda install numpy
conda install -c conda-forge fortran-compiler
```

NOTE: If you're building on OSX with the conda compilers it likely you'll need to export the compiler flags export LDFLAGS="-undefined dynamic_lookup -bundle"

Next do,

```
python setup.py bdist_wheel
```

to build a wheel. This wheel can be uploaded using twine:

```
twine upload ./dist/name_of_wheel
```

Linux wheels require a docker image. For more details see this blog https://uwekorn.com/2019/09/15/how-we-build-apache-arrows-manylinux-wheels.html


Unfortunately we cannot avoid a dynamic link to libquadmath on OSX. This appears to be a licencing issue with GCC and the conda gfortran package doesn't include the static version of quadmath.
So in order to use this package you'll need to have a system install of gfortran. We're not the only ones with this issue, e.g. https://github.com/aburrell/apexpy/issues/69 
