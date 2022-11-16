%module GWecc
%{
#include "GWecc/binary.hpp"
#include "GWecc/orbital_evolution.hpp"
#include "GWecc/eccentric_residuals.hpp"
#include "GWecc/GWecc.hpp"
%}

%include "std_vector.i"

namespace std {
   %template(DoubleVector) vector<double>;
   %template(DoubleVectorVector) vector<vector<double> >;
}

%include "GWecc/binary.hpp"
%include "GWecc/orbital_evolution.hpp"
%include "GWecc/eccentric_residuals.hpp"
%include "GWecc/GWecc.hpp"
