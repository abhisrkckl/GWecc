GWecc computes pulsar TOA delays due to gravitational waves from eccentric supermassive binary sources. 
enterprise_GWecc is a python wrapper for GWecc to be used with enterprise.

This code is based on Susobhanan et al. 2020. If you use this code in your work please cite the original article 

- Abhimanyu Susobhanan, Achamveedu Gopakumar, George Hobbs, and Stephen Taylor, 2020, *"Pulsar timing array signals induced by black hole binaries in relativistic eccentric orbits"*, Physical Review D, 101, 4,  043022, DOI: 10.1103/PhysRevD.101.043022, 	arXiv: 2002.03285



============
Dependencies
============

C/C++ Libraries
***************
* libgsl-dev
* libeigen3-dev

Python Packages
***************
* numpy
* enterprise
* setuptools

Other
*****
* swig

============
Installation
============

To install, type

> python3 setup.py install --user

Note that the C++ extension uses C++17 features and requires a new C++ compiler. I have tested this only with g++.

=====
Usage
=====

Examples of usage are given in the `examples/` directory.
