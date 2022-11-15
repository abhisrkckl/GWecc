%module GWecc
%{
#include "GWecc/residual_options.hpp"
#include "GWecc/binary.hpp"
#include "GWecc/OrbitalEvolution.hpp"
#include "GWecc/GWecc.hpp"
%}

%include "std_vector.i"

namespace std {
   %template(DoubleVector) vector<double>;
   %template(DoubleVectorVector) vector<vector<double> >;
}

%include "GWecc/residual_options.hpp"
%include "GWecc/binary.hpp"
%include "GWecc/OrbitalEvolution.hpp"
%include "GWecc/GWecc.hpp"
