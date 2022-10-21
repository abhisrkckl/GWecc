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

                si = ev_coeffs.sini,
                ci = ev_coeffs.cosi,
                                
                s2w = sin(2*w),
                c2w = cos(2*w),
                        
                OTS = sqrt(1-e*e),
                
                c2Om = ev_coeffs.cos2Omega,
                s2Om = ev_coeffs.sin2Omega,
                
                H0 = GWAmplitude(bin_mass, bin_now, DGW);
    
    /*const auto P = ((e + (-2 + e*e)*cu)*su)/(1 - ecu),
               Q = (OTS*(ecu - c2u))/(1 - ecu),
               R = esu;
    
    const auto A0 = R,
               A1 = P*c2w + Q*s2w,
               A2 = P*s2w - Q*c2w;*/
    
    const auto P = (OTS*(c2u - ecu))/(1 - ecu),
               Q = (((e*e - 2)*cu + e)*su)/(1 - ecu),
               R = esu;
    
    const auto A0 = R,
               A1 = -P*s2w + Q*c2w,
               A2 = P*c2w + Q*s2w;

    const auto a0 = si*si,
               a1 = 1+ci*ci,
               a2 = 2*ci;
    
    const auto sA = (H0/n) * ( a1*A1 + a0*A0 ),
               sB = (H0/n) * ( a2*A2 );
    
    const auto sp = c2Om*sA - s2Om*sB,
               sx = s2Om*sA + c2Om*sB;

    return std::make_tuple(sp, sx);
}

auto EccentricResidualsAndWaveform_px_fn_pt_Adb(const BinaryMass &bin_mass,
                                                const BinaryState &bin_init,
                                                const double DGW,
                                                const EvolveCoeffs_t& ev_coeffs,
                                                const double t){

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

                si = ev_coeffs.sini,
                ci = ev_coeffs.cosi,
                                
                s2w = sin(2*w),
                c2w = cos(2*w),

                s2phi = sin(2*phi),
                c2phi = cos(2*phi),
                        
                OTS = sqrt(1-e*e),
                
                c2Om = ev_coeffs.cos2Omega,
                s2Om = ev_coeffs.sin2Omega,
                
                H0 = GWAmplitude(bin_mass, bin_now, DGW);
    
    const auto P = (OTS*(c2u - ecu))/(1 - ecu),
               Q = (((e*e - 2)*cu + e)*su)/(1 - ecu),
               R = esu;
    
    const auto A0 = R,
               A1 = -P*s2w + Q*c2w,
               A2 = P*c2w + Q*s2w;

    const auto a0 = si*si,
               a1 = 1+ci*ci,
               a2 = 2*ci;
    
    const auto sA = (H0/n) * ( a1*A1 + a0*A0 ),
               sB = (H0/n) * ( a2*A2 );
    
    const auto sp = c2Om*sA - s2Om*sB,
               sx = s2Om*sA + c2Om*sB;

    const auto h_mq_A   = 1./ipow(1-ecu,2) * (    - (ci*ci+1)*(2*OTS*esu)         *s2phi
                            + (ci*ci+1)*(2*e*e - ecu*ecu + ecu - 2)    *c2phi
                            + (si*si)  *(1-ecu)*ecu          ),
               h_mq_B   = 1./ipow(1-ecu,2) * (             (2*OTS*esu)        *c2phi
                            +           (2*e*e - ecu*ecu + ecu - 2) *s2phi
                                    )*2*ci;
    
    const auto hA       = H0 * h_mq_A, 
               hB       = H0 * h_mq_B;
     
    const auto hp       = c2Om*hA - s2Om*hB,
               hx       = c2Om*hB + s2Om*hA;

    return std::make_tuple(sp, sx, hp, hx);
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

std::tuple<Signal1D,Signal1D> EccentricResidualsAndWaveform(const BinaryMass &bin_mass,
                                                            const BinaryState &bin_init,
                                                            const SkyPosition &bin_pos,
                                                            const SkyPosition &psr_pos,
                                                            const ResidualsTerms residuals_terms,
                                                            const Signal1D &ts){
    
    const auto [cosmu, Fp, Fx] = AntennaPattern(bin_pos, psr_pos);
    const auto DGW = bin_pos.DL;
    const auto ev_coeffs = compute_evolve_coeffs(bin_mass, bin_init);

    const auto nts = ts.size();
    Signal1D Rs(nts), hs(nts);
    for(unsigned idx=0; idx<nts; idx++){
        const auto& t = ts[idx];

        if(residuals_terms == ResidualsTerms::Earth){
            const auto [sp, sx, hp, hx] = EccentricResidualsAndWaveform_px_fn_pt_Adb(bin_mass, bin_init, DGW, ev_coeffs, t);
            Rs[idx] = -(Fp*sp + Fx*sx);
            hs[idx] = -(Fp*hp + Fx*hx);
        }
        else{
            const auto delay = -psr_pos.DL*(1-cosmu) / (1+bin_pos.z);
            if(residuals_terms == ResidualsTerms::Pulsar){
                const auto [sp, sx, hp, hx] = EccentricResidualsAndWaveform_px_fn_pt_Adb(bin_mass, bin_init, DGW, ev_coeffs, t+delay);
                Rs[idx] = (Fp*sp + Fx*sx);
                hs[idx] = (Fp*hp + Fx*hx);
            }
            else{
                const auto [spP, sxP, hpP, hxP] = EccentricResidualsAndWaveform_px_fn_pt_Adb(bin_mass, bin_init, DGW, ev_coeffs, t+delay);
                const auto [spE, sxE, hpE, hxE] = EccentricResidualsAndWaveform_px_fn_pt_Adb(bin_mass, bin_init, DGW, ev_coeffs, t);
                Rs[idx] = Fp*(spP-spE) + Fx*(sxP-sxE);
                hs[idx] = Fp*(hpP-hpE) + Fx*(hxP-hxE);
            }
        }
    }

    return {Rs, hs};
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
        Rps[idx] = spE;
        Rxs[idx] = sxE;
    }

    return std::make_tuple(Rps, Rxs);
}