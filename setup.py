from setuptools import setup

def readme():
      with open('README.rst') as f:
               return f.read()

setup(name='TEtranscripts',
      version='2.2.3',
      description='Tools for estimating differential enrichment of Transposable Elements and other highly repetitive regions',
      long_description=readme(),
      classifiers=[
          'Development Status :: 5 - Production/Stable',
          'Intended Audience :: Science/Research',
          'Environment :: Console',
          'Natural Language :: English',
          'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
          'Programming Language :: Python :: 2.7',
          'Programming Language :: Python :: 3',
          'Operating System :: MacOS :: MacOS X',
          'Operating System :: Unix'
      ],
      keywords='TE transposable element differential enrichment',
      url='http://hammelllab.labsites.cshl.edu/software#TEtranscripts',
      author='Ying Jin, Eric Paniagua, Talitha Forcier, Oliver Tam, Molly Hammell',
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
          'pysam>=0.9'
      ],
      include_package_data=True,
      zip_safe=False,
      scripts=[
          'bin/TEtranscripts',
          'bin/TEcount'
      ]
)
