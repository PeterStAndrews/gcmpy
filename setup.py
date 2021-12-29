# Setup for gcmpy
#
# Copyright (C) 2021 Peter Mann
#
# This file is part of gcmpy, generalised configuration model networks in Python.
#
# gcmpy is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# gcmpy is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with gcmpy. If not, see <http://www.gnu.org/licenses/gpl.html>.


from setuptools import setup

with open('README.rst') as f:
    longDescription = f.read()

setup(name = 'gcmpy',
      version = '0.0.2',
      description = 'Generalised Configuration Model random Graphs in Python',
      long_description = longDescription,
      url = 'https://github.com/PeterStAndrews/gcmpy',
      author = 'Peter Mann',
      author_email = 'pm78@st-andrews.ac.uk',
      license = 'License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)',
      classifiers = [ 'Intended Audience :: Science/Research',
                      'Intended Audience :: Developers',
                      'Programming Language :: Python :: 3.8',
                      'Programming Language :: Python :: 3.9',
                      'Topic :: Scientific/Engineering' ],
      python_requires = '>=3.8',
      packages = [ 'gcmpy' ],
      package_data = { 'gcmpy': [ 'py.typed' ] },
      zip_safe = False,
      install_requires = ["numpy >= 1.21.4"],
      extra_requires = { ':python_version < 3.8': [ 'typing_extensions' ] },
)
