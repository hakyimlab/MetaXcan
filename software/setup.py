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
                 packages=['metax','tests'],
                 license="TBD",
                 scripts=[  'GeneStats.py',
                            'M00_prerequisites.py',
                            'M01_covariances_correlations.py',
                            'M02_variances.py',
                            'M03_betas.py',
                            'M04_zscores.py',
                            'MetaXcan.py',
                            'MetaXcanUI.py',
                            'PPlots.py',
                            'PostProcessData.py',
                            'ProcessAlleleStatistics.py',
                            'ProcessZScoreFolder.py',
                            'VarianceComparison.py',
                            'scripts/BetaPlot.R',
                            'scripts/PlotGrid.R',
                            'scripts/PlotValuesWithSamples.R',
                            'scripts/Plots.R',
                            'scripts/QQUnifAllGrid.R',
                            'scripts/allele_stats.R'],
                 description=["TBD"],
                 install_requires=['scipy', 'numpy'],
                 long_description=read('Readme.md'),
                 keywords=['TBD'],
                 test_suite='tests',
                 package_data={'tests/files':['*']},
                 classifiers=[
                    "Development Status :: 4 - Beta",
                    "Topic :: Utilities",
                    "Topic :: Scientific/Engineering :: Bio-Informatics",
                    "Topic :: Software Development :: Libraries",
                    "Programming Language :: Python :: 2.7"
                 ],
)
