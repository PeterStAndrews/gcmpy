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


from setuptools import setup, find_packages

# the current version of the package
__version__ = "0.2.21"

with open('README.rst') as f:
    longDescription = f.read()

setup(
    name='gcmpy',
    version=__version__,
    description='Generalised Configuration Model random graphs in Python',
    long_description=longDescription,
    url='https://github.com/PeterStAndrews/gcmpy',
    author='Peter Mann',
    author_email='pm78@st-andrews.ac.uk',
    license='License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)',
    classifiers=[
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Topic :: Scientific/Engineering'
    ],
    python_requires='>=3.8',
    packages=['gcmpy', 'gcmpy.covers','gcmpy.distributions','gcmpy.gcm_algorithm',
    'gcmpy.joint_degree','gcmpy.joint_degree.joint_degree_loaders', 'gcmpy.motif_generators','gcmpy.names','gcmpy.network', 'gcmpy.tools'],
    package_data={'gcmpy': ['py.typed']},
    zip_safe=False,
    install_requires=[
        "iteration_utilities==0.11.0",
        "networkx==2.6.3",
        "numpy==1.22.3"
    ]
)
