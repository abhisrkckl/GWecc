#include "NumericalWaveform.hpp"
#include "AntennaPattern.hpp"
#include "mikkola.h"
#include "PN.hpp"
#include "ipow.hpp"
#include <iostream>

auto EccentricResiduals_px_fn_pt_Adb(const BinaryMass &bin_mass,
                                     const BinaryState &bin_init,
                                     const double DGW,
                                     const EvolveCoeffs_t& ev_coeffs,
                                     const double t){
    
    //const auto [cosmu, Fp, Fx] = AntennaPattern(bin_pos, psr_pos);
    //const auto ev_coeffs = compute_evolve_coeffs(bin_mass, bin_init);

    const auto bin_now = solve_orbit_equations(bin_init, ev_coeffs, t-bin_init.t);

    if(bin_now.merged){
        throw std::runtime_error("The binary has already merged.");
    }

    const auto  &n = bin_now.n,
                &e = bin_now.e,
                &g = bin_now.gamma,
                &l = bin_now.l,
                        
                u = MIKKOLA(l,e),
                L = l+g,
                
                su = sin(u),
                cu = cos(u),
                        
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

                si = ev_coeffs.sini,
                ci = ev_coeffs.cosi,
                                
                s2w = sin(2*w),
                c2w = cos(2*w),
                        
                OTS = sqrt(1-e*e),
                
                c2Om = ev_coeffs.cos2Omega,
                s2Om = ev_coeffs.sin2Omega,
                
                H0 = GWAmplitude(bin_mass, bin_now, DGW);
    
    const auto P = (esu*esu + ecu - 1)*2*OTS/(e*e)/(1-ecu),
               Q = (2*ecu/(e*e) - ecu - 1)*esu/(1-ecu),
               R = esu;
    
    const auto sA = (H0/n) * (  (1+ci*ci)*(P*s2w - Q*c2w)  +  si*si*R  ),
               sB = (H0/n) * (      -2*ci*(P*c2w + Q*s2w));
    
    const auto sp = c2Om*sA - s2Om*sB,
               sx = c2Om*sB + s2Om*sA;

    return std::make_tuple(sp, sx);
}

auto EccentricResiduals_fn_pt_Adb(const BinaryMass &bin_mass,
                                  const BinaryState &bin_init,
                                  const double Fp, const double Fx, const double DGW,
                                  const EvolveCoeffs_t& ev_coeffs,
                                  const double t) {
    const auto [sp, sx] = EccentricResiduals_px_fn_pt_Adb(bin_mass, bin_init, DGW, ev_coeffs, t);
    return Fp*sp + Fx*sx;
}

Signal1D EccentricResiduals_Adb(const BinaryMass &bin_mass,
                                const BinaryState &bin_init,
                                const SkyPosition &bin_pos,
                                const SkyPosition &psr_pos,
                                const ResidualsTerms residuals_terms,
                                const Signal1D &ts){

    const auto [cosmu, Fp, Fx] = AntennaPattern(bin_pos, psr_pos);
    const auto DGW = bin_pos.DL;
    const auto ev_coeffs = compute_evolve_coeffs(bin_mass, bin_init);

    const auto nts = ts.size();
    Signal1D Rs(nts);
    for(unsigned idx=0; idx<nts; idx++){
        const auto& t = ts[idx];

        if(residuals_terms == ResidualsTerms::Earth){
            Rs[idx] = -EccentricResiduals_fn_pt_Adb(bin_mass, bin_init, Fp, Fx, DGW, ev_coeffs, t);
        }
        else{
            const auto delay = -psr_pos.DL*(1-cosmu) / (1+bin_pos.z);
            if(residuals_terms == ResidualsTerms::Pulsar){
                Rs[idx] = EccentricResiduals_fn_pt_Adb(bin_mass, bin_init, Fp, Fx, DGW, ev_coeffs, t + delay);
            }
            else{
                Rs[idx] =   EccentricResiduals_fn_pt_Adb(bin_mass, bin_init, Fp, Fx, DGW, ev_coeffs, t + delay) 
                          - EccentricResiduals_fn_pt_Adb(bin_mass, bin_init, Fp, Fx, DGW, ev_coeffs, t);
            }
        }
    }

    return Rs;
}

std::tuple<Signal1D, Signal1D> EccentricResiduals_px_Adb(const BinaryMass &bin_mass,
                                                         const BinaryState &bin_init,
                                                         const double DGW,
                                                         const Signal1D &ts){

    const auto ev_coeffs = compute_evolve_coeffs(bin_mass, bin_init);

    const auto nts = ts.size();
    Signal1D Rps(nts), Rxs(nts);
    for(unsigned idx=0; idx<nts; idx++){
        const auto [spE, sxE] = EccentricResiduals_px_fn_pt_Adb(bin_mass, bin_init, DGW, ev_coeffs, ts[idx]);
        Rps[idx] = -spE;
        Rxs[idx] = -sxE;
    }

    return std::make_tuple(Rps, Rxs);
}