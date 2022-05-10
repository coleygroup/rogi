"""Measures of roughness for molecular property landscapes
"""

import versioneer
from setuptools import setup


# readme file
def readme():
    with open('README.md') as f:
        return f.read()


# -----
# Setup
# -----
setup(name='rogi',
      version=versioneer.get_version(),
      cmdclass=versioneer.get_cmdclass(),
      description='Measures of roughness for molecular property landscapes',
      long_description=readme(),
      long_description_content_type='text/markdown',
      classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',
      ],
      url='https://github.com/coleygroup/rogi',
      author='Matteo Aldeghi',
      author_email='maldeghi@mit.edu',
      license='MIT',
      packages=['rogi'],
      package_dir={'': 'src'},
      zip_safe=False,
      tests_require=['pytest'],
      install_requires=['numpy', 'scipy>=1.4', 'fastcluster', 'pandas', 'scikit-learn>=1'],
      python_requires=">=3.7",
      ext_modules=[]
      )
