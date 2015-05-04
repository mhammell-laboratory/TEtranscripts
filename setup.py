from setuptools import setup

try:
      from pypandoc import convert
      read_md = lambda f: convert(f, 'rst')
except ImportError:
      print("warning: pypandoc module not found, could not convert Markdown to RST")
      read_md = lambda f: open(f, 'r').read()

setup(name='TEToolkit',
      version='X.X.X',
      description='Tools for estimating differential enrichment of Transposable Elements and other highly repetitive regions',
      long_description=readme('README.md'),
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
