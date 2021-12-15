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

ResNorm = numpy.distutils.core.Extension(name="quasielasticbayes.ResNorm",
                                         sources=['quasielasticbayes/ResNorm_main.f90',
                                                  'quasielasticbayes/ResNorm_subs.f90',
                                                  'quasielasticbayes/BlrRes.f90',
                                                  'quasielasticbayes/Bayes.f90',
                                                  'quasielasticbayes/Four.f90',
                                                  'quasielasticbayes/Util.f90'],
                                         extra_link_args=extra_link_args,
                                         extra_compile_args=extra_compiler_args)

Quest = numpy.distutils.core.Extension(name="quasielasticbayes.Quest",
                                       sources=['quasielasticbayes/Quest_main.f90',
                                                'quasielasticbayes/Quest_subs.f90',
                                                'quasielasticbayes/BlrRes.f90',
                                                'quasielasticbayes/Bayes.f90',
                                                'quasielasticbayes/Four.f90',
                                                'quasielasticbayes/Util.f90',
                                                'quasielasticbayes/Simopt.f90'],
                                       extra_link_args=extra_link_args,
                                       extra_compile_args=extra_compiler_args)

QLse = numpy.distutils.core.Extension(name="quasielasticbayes.QLse",
                                      sources=['quasielasticbayes/QLse_main.f90',
                                               'quasielasticbayes/QLse_subs.f90',
                                               'quasielasticbayes/BlrRes.f90',
                                               'quasielasticbayes/Bayes.f90',
                                               'quasielasticbayes/Four.f90',
                                               'quasielasticbayes/Util.f90',
                                               'quasielasticbayes/Simopt.f90'],
                                        extra_link_args=extra_link_args,
                                        extra_compile_args=extra_compiler_args)

QLres = numpy.distutils.core.Extension(name="quasielasticbayes.QLres",
                                       sources=['quasielasticbayes/QLres_main.f90',
                                                'quasielasticbayes/QLres_subs.f90',
                                                'quasielasticbayes/BlrRes.f90',
                                                'quasielasticbayes/Bayes.f90',
                                                'quasielasticbayes/Four.f90',
                                                'quasielasticbayes/Util.f90'],
                                        extra_link_args=extra_link_args,
                                        extra_compile_args=extra_compiler_args)

QLdata = numpy.distutils.core.Extension(name="quasielasticbayes.QLdata",
                                        sources=['quasielasticbayes/QLdata_main.f90',
                                                 'quasielasticbayes/QLdata_subs.f90',
                                                 'quasielasticbayes/Bayes.f90',
                                                 'quasielasticbayes/Four.f90',
                                                 'quasielasticbayes/Util.f90'],
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
    name='quasielasticbayes',
    install_requires=['numpy>=1.17.5'],
    packages=['quasielasticbayes'],
    description='A Bayesian fitting package used for fitting quasi-elastic neutron scattering data.',
    long_description='This package wraps fortran Bayesian fitting libraries using f2py. '\
                     'An application of this package is to fit quasi elastic neutron scattering data in Mantid (https://www.mantidproject.org)',
    author='Mantid Team',
    ext_modules=[ResNorm, Quest, QLse, QLres, QLdata],
    author_email="mantid-help@mantidproject.org",
    url='https://www.mantidproject.org',
    version="0.1.0",
    license='BSD',
    cmdclass={'build_ext': build_ext}
)
