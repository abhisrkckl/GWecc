#include <cmath>
#include "binary.hpp"
#include "mikkola.h"
#include "post_newtonian.hpp"
#include "ipow.hpp"
#include "waveform_vars.hpp"

std::tuple<SinCos, double, double> get_sincosu_phi_omega(const BinaryMass& bin_mass, const BinaryState& bin_now){
    const auto  e = bin_now.e,
                g = bin_now.gamma,
                l = bin_now.l;

    const double u = mikkola(l, e); 
    const SinCos scu{u};
    const auto [su, cu] = scu;

    const auto  L = l+g,
                
                k = advance_of_periastron(bin_mass, bin_now),
                ephi = angular_eccentricity(bin_mass, bin_now),
                
                betaphi = (ephi>1e-15)  ? (1-sqrt(1-ephi*ephi))/ephi 
                                        : (ephi/2. + ipow(ephi,3)/8. + ipow(ephi,5)/16.),
                v_u = 2*atan2(betaphi*su, 1-betaphi*cu),
                v_l = v_u + e*su,
                v = v_l + l,
                        
                W = (1+k)*v_l,
                phi = L+W,
                w = phi-v;
    
    return {{u}, phi, w};
}

std::array<double, 3> get_residual_PQR(double e, SinCos scu){
    const auto [su, cu] = scu;
    
    const double c2u = cu*cu - su*su,
                 OTS = sqrt(1-e*e);
    
    const double P = (OTS*(c2u - e*cu))/(1 - e*cu),
                 Q = (((e*e - 2)*cu + e)*su)/(1 - e*cu),
                 R = e*su;
    
    return {P, Q, R};
}

std::array<double, 3> get_residual_A012(double e, SinCos scu, double w){

    const auto [P, Q, R] = get_residual_PQR(e, scu);

    const auto [s2w, c2w] = SinCos(2*w);
    
    const auto A0 = R,
               A1 = -P*s2w + Q*c2w,
               A2 = P*c2w + Q*s2w;
    
    return {A0, A1, A2};
}

std::array<double, 3> get_residual_a012(const double ci){
    return {1-ci*ci, 1+ci*ci, 2*ci};
}

std::array<double, 2> get_residual_sAB(const std::array<double,3>& A, const std::array<double,3>& a, const double S){
    return {S*(a[1]*A[1] + a[0]*A[0]), S*a[2]*A[2]};
}

std::array<double, 2> get_residual_spx(const std::array<double,2>& sAB, const SinCos cs2psi){
    const auto [sA, sB] = sAB;
    const auto [s2psi, c2psi] = cs2psi;
    const auto sp = c2psi*sA - s2psi*sB,
               sx = s2psi*sA + c2psi*sB;
    return {sp, sx};
}

std::array<double, 3> get_waveform_PQR(double e, SinCos scu){
    const auto OTS = sqrt(1-e*e);
    const auto [su, cu] = scu;
    const auto denom = ipow(1-e*cu, 2);
    const auto P = 2*OTS*e*su / denom,
               Q = (2*e*e - e*cu*e*cu + e*cu - 2) / denom,
               R = (1-e*cu)*e*cu / denom;
    return {P, Q, R};
}

std::array<double, 3> get_waveform_A012(double e, SinCos scu, double phi){
    const auto [P, Q, R] = get_waveform_PQR(e, scu);
    const auto [s2phi, c2phi] = SinCos(2*phi);
    const auto A1 = -P*s2phi + Q*c2phi,
               A2 =  P*c2phi + Q*s2phi,
               A0 =  R;  
    return {A0, A1, A2};
}

