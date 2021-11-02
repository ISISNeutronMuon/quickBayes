import setuptools
import numpy.distutils.core
import sys

if sys.platform == 'win32':
    # Compile static libs on windows
    extra_link_args = ["-static", "-static-libgfortran", "-static-libgcc"]
elif sys.platform == 'darwin':
    extra_compiler_args = ['-Wno-argument-mismatch']
    extra_link_args = ["-static", "-static-libgfortran", "-static-libgcc"]
else:
    # In the docker container there are no issues building the libraries, so we can assume we don't need any flags
    extra_compiler_args = []
    extra_link_args = []   

ResNorm = numpy.distutils.core.Extension(name="bayesfitting.ResNorm",
                                         sources=['bayesfitting/ResNorm_main.f90',
                                                  'bayesfitting/ResNorm_subs.f90',
                                                  'bayesfitting/BlrRes.f90',
                                                  'bayesfitting/Bayes.f90',
                                                  'bayesfitting/Four.f90',
                                                  'bayesfitting/Util.f90'],
                                         extra_link_args=extra_link_args,
                                         extra_compile_args=extra_compiler_args)

Quest = numpy.distutils.core.Extension(name="bayesfitting.Quest",
                                       sources=['bayesfitting/Quest_main.f90',
                                                'bayesfitting/Quest_subs.f90',
                                                'bayesfitting/BlrRes.f90',
                                                'bayesfitting/Bayes.f90',
                                                'bayesfitting/Four.f90',
                                                'bayesfitting/Util.f90',
                                                'bayesfitting/Simopt.f90'],
                                       extra_link_args=extra_link_args,
                                       extra_compile_args=extra_compiler_args)

QLse = numpy.distutils.core.Extension(name="bayesfitting.QLse",
                                      sources=['bayesfitting/QLse_main.f90',
                                               'bayesfitting/QLse_subs.f90',
                                               'bayesfitting/BlrRes.f90',
                                               'bayesfitting/Bayes.f90',
                                               'bayesfitting/Four.f90',
                                               'bayesfitting/Util.f90',
                                               'bayesfitting/Simopt.f90'],
                                        extra_link_args=extra_link_args,
                                        extra_compile_args=extra_compiler_args)

QLres = numpy.distutils.core.Extension(name="bayesfitting.QLres",
                                       sources=['bayesfitting/QLres_main.f90',
                                                'bayesfitting/QLres_subs.f90',
                                                'bayesfitting/BlrRes.f90',
                                                'bayesfitting/Bayes.f90',
                                                'bayesfitting/Four.f90',
                                                'bayesfitting/Util.f90'],
                                        extra_link_args=extra_link_args,
                                        extra_compile_args=extra_compiler_args)

QLdata = numpy.distutils.core.Extension(name="bayesfitting.QLdata",
                                        sources=['bayesfitting/QLdata_main.f90',
                                                 'bayesfitting/QLdata_subs.f90',
                                                 'bayesfitting/Bayes.f90',
                                                 'bayesfitting/Four.f90',
                                                 'bayesfitting/Util.f90'],
                                        extra_link_args=extra_link_args,
                                        extra_compile_args=extra_compiler_args)

QLwat = numpy.distutils.core.Extension(name="bayesfitting.QLwat",
                                        sources=['bayesfitting/QLwat_main.f90',
                                                 'bayesfitting/QLwat_subs.f90',
                                                 'bayesfitting/BlrRes.f90',
                                                 'bayesfitting/Bayes.f90',
                                                 'bayesfitting/Four.f90',
                                                 'bayesfitting/Util.f90'],
                                        extra_link_args=extra_link_args,
                                        extra_compile_args=extra_compiler_args)

from numpy.distutils.command.build_ext import build_ext as _build_ext
class build_ext(_build_ext):
    def finalize_options(self):
        _build_ext.finalize_options(self)
        # If we don't do this on windows, when we do bdist_wheel we wont get a static link
        # this is because it misses the compiler flags to f2py which means it ignores the static flags we try to pass
        if sys.platform == 'win32':
            self.fcompiler = 'gnu95'
            self.compiler = 'mingw32'

numpy.distutils.core.setup(
    name='bayesfitting',
    install_requires=['numpy>=1.17.5'],
    packages=['bayesfitting'],
    description='A Bayesian fitting package used for fitting quasi-elastic neutron scattering data.',
    long_description='This package wraps fortran Bayesian fitting libraries using f2py. '\
                     'An application of this package is to fit quasi elastic neutron scattering data in Mantid (https://www.mantidproject.org)',
    author='Mantid Team',
    ext_modules=[ResNorm, Quest, QLse, QLres, QLdata, QLwat],
    author_email="mantid-help@mantidproject.org",
    url='https://www.mantidproject.org',
    version="0.2.0",
    license='BSD',
    cmdclass={'build_ext': build_ext}
)
