# Building the f2py python extensions

You can build this package using the setup.py file. To use this, create a minimal conda environment as shown below:

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

Next do,

```
python setup.py bdist_wheel
```

to build a wheel. This wheel can be uploaded using twine:

```
twine upload ./dist/name_of_wheel
```

Linux wheels require a docker image. For more details see this blog https://uwekorn.com/2019/09/15/how-we-build-apache-arrows-manylinux-wheels.html
