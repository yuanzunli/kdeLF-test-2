import os
import sys

version = {}
with open('kdeLF/__init__.py') as kde:
    exec(kde.read(), version)
__version__ = version['__version__']

try:
    import setuptools
except ImportError:
    pass

from numpy.distutils.core import setup, Extension
from numpy.distutils.system_info import get_info

# Fortran extension
fortran_t = Extension(name = 'kdeLF.kde_fortran_t', 
		sources = ['kdeLF/kdeLF.f90'])

description = 'A flexible method for estimating luminosity functions via Kernel Density Estimation'

this_directory = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

if __name__ == "__main__":
	setup(name = 'kdeLF',     
	      packages=['kdeLF'],
	      ext_modules=[fortran_t],
	      version=__version__,
	      description=description,
	      #long_description = '', # 
	      long_description=long_description,    #
	      long_description_content_type="text/markdown",
	      author='Zunli Yuan',
	      author_email='yzl@hunnu.edu.cn',
              url='https://github.com/yuanzunli/kdeLF',
	      install_requires=[
	      			'numpy >= 1.22',
	      			'scipy',
	      			'astropy',
	      			'emcee',
	      			'matplotlib'
	      			],
	      python_requires=">=3.6")
#!/usr/bin/env python
