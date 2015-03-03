from setuptools import setup

def readme():
    with open('README') as f:
        return f.read()

setup(name='TEToolkit',
      version='1.3',
      description='Tools for estimating differential enrichment of Transposable Elements or other highly repetetive regions',
      long_description=readme(),
      classifiers=[
      ],
      keywords='TE transposable element differential enrichment',
      url='http://www.example.org',
      author='Ying Jin, Oliver Tam',
      author_email='yjin@cshl.edu',
      license='MIT',
      packages=[
          'TEToolkit',
          'TEToolkit.IO',
          'TEToolkit.ShortRead'
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
