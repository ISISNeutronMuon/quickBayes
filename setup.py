from setuptools import find_packages, setup
from tools.setup_helper import get_extensions


VERSION = "1.0.0b19"


PACKAGE_NAME = 'quickBayes'


extensions = get_extensions(PACKAGE_NAME)


setup(
    name=PACKAGE_NAME,
    requires=['numpy'],
    setup_requires=['numpy>=1.12'],
    install_requires=['numpy>=1.12', 'scipy', 'gofit'],
    packages=find_packages(where='src'),
    description='A Bayesian fitting package used for '
                'model selection and grid searches of fits '
                'for neutron and muon data.',
    long_description='This package provides code for a Bayesian '
                     'workflow. The two options are '
                     'model selection and grid search. '
                     'This package replaces quasielasticbayes. '
                     'An application of this package is '
                     'to fit quasi-elastic '
                     'neutron scattering data in Mantid '
                     '(https://www.mantidproject.org)',
    author='Anthony Lim',
    ext_modules=extensions,
    author_email="anthony.lim@stfc.ac.uk",
    url='https://www.mantidproject.org',
    version=VERSION,
    license='BSD',
    package_dir={'': 'src'}
)
