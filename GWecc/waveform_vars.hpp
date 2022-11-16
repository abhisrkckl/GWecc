#ifndef _waveform_vars_hpp_
#define _waveform_vars_hpp_

#include <cmath>
#include <array>
#include <tuple>

#include "binary.hpp"
#include "mikkola.h"
#include "post_newtonian.hpp"
#include "ipow.hpp"

struct SinCos{
    const double s, c;
    
    SinCos(double x) : s(sin(x)), c(cos(x)) {}
};

std::tuple<SinCos, double, double> get_sincosu_phi_omega(const BinaryMass& bin_mass, const BinaryState& bin_now);

std::array<double, 3> get_residual_PQR(double e, SinCos scu);

std::array<double, 3> get_residual_A012(double e, SinCos scu, double w);

std::array<double, 3> get_residual_a012(const double ci);

std::array<double, 2> get_residual_sAB(const std::array<double,3>& A, const std::array<double,3>& a, const double S);

std::array<double, 2> get_residual_spx(const std::array<double,2>& sAB, const SinCos cs2psi);

std::array<double, 3> get_waveform_A012(double e, SinCos scu, double phi);

#endif