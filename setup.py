from setuptools import setup

def readme():
      with open('README.rst') as f:
	       return f.read()

setup(name='TEToolkit',
      version='1.4.11',
      description='Tools for estimating differential enrichment of Transposable Elements and other highly repetitive regions',
      long_description=readme(),
      classifiers=[
      ],
      keywords='TE transposable element differential enrichment',
      url='http://hammelllab.labsites.cshl.edu/software#TEToolkit',
      author='Ying Jin, Eric Paniagua, Oliver Tam and Molly Hammell',
      author_email='yjin@cshl.edu',
      license='GPLv3',
      packages=[
          'TEToolkit',
          'TEToolkit.IO',
          'TEToolkit.ShortRead'
      ],
      platforms=[
          'Linux',
          'MacOS'
      ],
      install_requires=[
          'argparse',
          'pysam>=0.8'
      ],

      include_package_data=True,
      zip_safe=False,
      scripts=[
          'bin/TEtranscripts',
          'bin/TEpeaks'
      ]
)
