GWecc computes pulsar TOA delays due to gravitational waves from eccentric supermassive binary sources. 
enterprise_GWecc is a python wrapper for GWecc to be used with enterprise.

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
> python3 setup.py install

Note that the C++ extension uses C++17 features and requires a new C++ compiler.
I have tested this only with g++.

=====
Usage
=====
Exposes the function "enterprise_GWecc.eccentric_cw_delay(...)".
See examples/Example.py for usage.