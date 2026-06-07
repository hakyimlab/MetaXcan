import ez_setup
ez_setup.use_setuptools()

import setuptools
import os

# Use app's version
import metax

# Use the readme as the long description baked into the application itself
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setuptools.setup(name="MetaXcan",
                 version=metax.__version__,
                 author="Alvaro Barbeira, Eric Torstenson",
                 author_email='alvarobarbeira@gmail.com, eric.s.torstenson@vanderbilt.edu',
                 url="TBD",
                 packages=['metax', 'tests','metax.misc', 'metax.gwas','metax.metaxcan', 'metax.deprecated'],
                 license="TBD",
                 scripts=[  'M00_prerequisites.py',
                            'M01_covariances_correlations.py',
                            'M02_variances.py',
                            'M03_betas.py',
                            'M04_zscores.py',
                            'PrediXcan.py',
                            'SPrediXcan.py',
                            'MetaMany.py',
                            'MulTiXcan.py',
                            'SMulTiXcan.py'],
                 description=["TBD"],
                 install_requires=[
                     'numpy>=1.23,<2',
                     'scipy>=1.10,<2',
                     'pandas>=2.0,<3',
                     'patsy>=0.5.6',
                     'statsmodels>=0.14,<0.15',
                     'h5py>=3.10,<4',
                     'pyliftover>=0.4',
                     'bgen==1.9.10',
                     'cyvcf2>=0.30',
                 ],
                 extras_require={"test": ["sqlalchemy"]},
                 long_description=read('Readme.md'),
                 keywords=['TBD'],
                 test_suite='tests',
                 package_data={'tests/files':['*']},
                 python_requires='>=3.9',
                 classifiers=[
                    "Development Status :: 4 - Beta",
                    "Topic :: Utilities",
                    "Topic :: Scientific/Engineering :: Bio-Informatics",
                    "Topic :: Software Development :: Libraries",
                    "Programming Language :: Python :: 2.7"
                 ],
)
