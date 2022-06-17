GWecc computes pulsar TOA delays due to gravitational waves from eccentric supermassive binary sources. 
enterprise_GWecc is a python wrapper for GWecc to be used with enterprise.

This code is based on Susobhanan et al. 2020. If you use this code in your work please cite the original article 

- Abhimanyu Susobhanan, Achamveedu Gopakumar, George Hobbs, and Stephen Taylor, 2020, *"Pulsar timing array signals induced by black hole binaries in relativistic eccentric orbits"*, Physical Review D, 101, 4,  043022, DOI: 10.1103/PhysRevD.101.043022, 	arXiv: 2002.03285



============
Dependencies
============

C/C++ Libraries
***************
* GSL
* Eigen

Python Packages
***************
* numpy
* enterprise-pulsar
* setuptools

Other
*****
* swig

============
Installation
============

I suggest installing this in a conda environment.
To install, type

$ conda install enterprise-pulsar gsl eigen

$ pip install .

Note that the C++ extension uses C++17 features and requires a new C++ compiler. I have tested this only with g++-8 and later.

=====
Usage
=====

Examples of usage are given in the `examples/` directory.

================
Acknowledgements
================
Lankeswar Dey, Yannick Boetzel, Belinda Cheeseboro, Nidhi Pant, Amit Jit Singh
