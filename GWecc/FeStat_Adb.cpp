#include "EccentricResiduals.hpp"
#include "OrbitalEvolution.hpp"
#include "antenna_pattern.hpp"
#include "NumericalWaveform.hpp"
#include "PN.hpp"
#include "mikkola.h"
#include "ipow.hpp"
#include "FeStat.hpp"
#include <sstream>

std::array<double, 3> fe_stat_A_fn_pt_Adb(double t, void *_params){

    const auto &wf_params = *static_cast<WaveformParams*>(_params);
   
    const auto &bin_mass = wf_params.bin_mass;
    const auto &bin_init = wf_params.bin_init;
    const auto &ev_coeffs = wf_params.ev_coeffs;
    
    const auto bin_now = solve_orbit_equations(bin_init, ev_coeffs, t-bin_init.t);

    if(bin_now.merged){
        return {0,0,0};
    }
    
    const auto  &n = bin_now.n,
                &e = bin_now.e,
                &g = bin_now.gamma,
                &l = bin_now.l,
                        
                u = mikkola(l,e),
                L = l+g,
                
                su = sin(u),
                cu = cos(u),
                c2u = cu*cu - su*su,
                        
                esu = e*su,
                ecu = e*cu,
                
                k = advance_of_periastron(bin_mass, bin_now),
                ephi = angular_eccentricity(bin_mass, bin_now),
                
                betaphi = (ephi>1e-15)  ? (1-sqrt(1-ephi*ephi))/ephi 
                                        : (e/2. + ipow(e,3)/8. + ipow(e,5)/16.),
                v_u = 2*atan2(betaphi*su, 1-betaphi*cu),
                v_l = v_u + esu,
                v = v_l + l,
                        
                W = (1+k)*v_l,
                phi = L+W,
                w = phi-v,
                                
                s2w = sin(2*w),
                c2w = cos(2*w),
                        
                OTS = sqrt(1-e*e);
    
    const double P = OTS*(c2u - e*cu)/(1-e*cu),
                 Q = ((e*e-2)*cu + e)*su/(1-e*cu),
                 R = e*su;

    return {P*c2w + Q*s2w, P*s2w - Q*c2w, R};
}

std::array<Signal1D, 6> fe_stat_funcs_Adb(const BinaryMass &bin_mass,
                                        const BinaryState &bin_init,
                                        const SkyPosition &bin_pos,
                                        const SkyPosition &psr_pos,
                                        const Signal1D &ts){

    const auto [cosmu, Fp, Fx] = antenna_pattern(bin_pos, psr_pos);

    const EvolveCoeffs_t ev_coeffs = compute_evolve_coeffs(bin_mass, bin_init);

    WaveformParams params { Fp, Fx, 1,
                            bin_mass,
                            bin_init,
                            ev_coeffs };
    
    const size_t nts = ts.size();
    Signal1D f1s(nts), f2s(nts), f3s(nts);
    for(unsigned idx=0; idx<nts; idx++){
        const auto [f1, f2, f3] = fe_stat_A_fn_pt_Adb(ts[idx], &params);
        f1s[idx] = f1;
        f2s[idx] = f2;
        f3s[idx] = f3;
    }    

    const double H0 = gw_amplitude(bin_mass, bin_init, bin_pos.DL);

    return {H0*Fp*f1s, H0*Fp*f2s, H0*Fp*f3s,
            H0*Fx*f1s, H0*Fx*f2s, H0*Fx*f3s};
}
