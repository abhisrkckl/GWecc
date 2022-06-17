%module GWecc
%{
#include "GWecc/GWecc.hpp"
#include "GWecc/ResidualsOptions.hpp"
#include "GWecc/Binary.hpp"
#include "GWecc/OrbitalEvolution.hpp"
%}

%include "std_vector.i"

namespace std {
   %template(DoubleVector) vector<double>;
   %template(DoubleVectorVector) vector<vector<double> >;
}

%include "GWecc/GWecc.hpp"
%include "GWecc/ResidualsOptions.hpp"
%include "GWecc/Binary.hpp"
%include "GWecc/OrbitalEvolution.hpp"
