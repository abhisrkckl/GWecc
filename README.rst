GWecc computes pulsar TOA delays due to gravitational waves from eccentric supermassive binary sources. 
enterprise_GWecc is a python wrapper for GWecc to be used with enterprise.

This code is based on Susobhanan et al. 2020 and Susobhanan 2022. If you use this code in your work please cite the original articles 

- Abhimanyu Susobhanan, Achamveedu Gopakumar, George Hobbs, and Stephen Taylor, 2020, *"Pulsar timing array signals induced by black hole binaries in relativistic eccentric orbits"*, Physical Review D, 101, 4,  043022, DOI: 10.1103/PhysRevD.101.043022, 	arXiv:2002.03285
- Abhimanyu Susobhanan, 2022, *"Post-Newtonian-accurate pulsar timing array signals induced by inspiralling eccentric binaries: accuracy and computational cost"*, arXiv e-prints, arXiv:2210.11454

The code itself can be cited using the ASCL record

- Abhimanyu Susobhanan, *"GWecc: Calculator for pulsar timing array signals due to eccentric supermassive binaries"*, Astrophysics Source Code Library, ascl:2002.013, https://github.com/abhisrkckl/GWecc


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
* astropy
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

``$ conda install -c conda-forge enterprise-pulsar gsl eigen``

``$ pip install .``

Note that the C++ extension uses C++17 features and requires a new C++ compiler. I have tested this only with g++-8 and later.

=====
Usage
=====

Examples of usage are given in the `examples` directory.

================
Acknowledgements
================
Lankeswar Dey, Yannick Boetzel, Belinda Cheeseboro, Nidhi Pant, Amit Jit Singh
