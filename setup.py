from numpy.distutils.core import setup
from v2_setup import get_v2_extensions


VERSION = "1.0.0a21"
PACKAGE_NAME = 'quasielasticbayes'


extensions = (get_v2_extensions(PACKAGE_NAME))


setup(
    name=PACKAGE_NAME,
    install_requires=['numpy>=1.12', 'scipy', 'gofit'],
    packages=[PACKAGE_NAME],
    description='A Bayesian fitting package used for '
                'fitting quasi-elastic neutron scattering data.',
    long_description='This package wraps fortran Bayesian '
                     'fitting libraries using f2py. '
                     'An application of this package is '
                     'to fit quasi-elastic '
                     'neutron scattering data in Mantid '
                     '(https://www.mantidproject.org)',
    author='Mantid Team',
    ext_modules=extensions,
    author_email="mantid-help@mantidproject.org",
    url='https://www.mantidproject.org',
    version=VERSION,
    license='BSD',
    package_dir={'': 'src'},  # allows setup to find py and f90 files
)
