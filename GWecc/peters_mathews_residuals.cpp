#include <iostream>
#include <gsl/gsl_sf_bessel.h>
#include "eccentric_residuals.hpp"
#include "orbital_evolution.hpp"
#include "antenna_pattern.hpp"
#include "mikkola.h"
#include "post_newtonian.hpp"
#include "ipow.hpp"

int get_nharm(const double e){
    if(e < 0.009){
        return 2;
    }  
    else if(e < 0.063){
        return 3;
    }
    else{
        double OTS = sqrt(1-e*e),         
               OTS3 = OTS*OTS*OTS;
        return int(round(18.64801851/OTS3 - 14.04695398));
    }
}

auto eccentric_residuals_px_fn_pt_PM(const BinaryMass &bin_mass,
                                    const BinaryState &bin_init,
                                    const double DGW,
                                    const EvolveCoeffs& ev_coeffs,
                                    const double t){

    const auto bin_now = solve_orbit_equations(bin_init, ev_coeffs, t-bin_init.t);

    if(bin_now.merged){
        throw std::runtime_error("The binary has already merged.");
    }

    const auto  &n = bin_now.n,
                &e = bin_now.e,
                &g = bin_now.gamma,
                &l = bin_now.l,
                        
                u = mikkola(l,e),
                L = l+g,
                
                su = sin(u),
                cu = cos(u),
                        
                esu = e*su,
                
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
                
                c2psi = ev_coeffs.cos2psi,
                s2psi = ev_coeffs.sin2psi,
                
                H0 = gw_amplitude(bin_mass, bin_now, DGW);

    const int pmax = get_nharm(e);
    double P=0, Q=0, R=0;
    for(int p=1; p<=pmax; p++){
        const double Jpm2 = gsl_sf_bessel_Jn(p-2, p*e),
                     Jpm1 = gsl_sf_bessel_Jn(p-1, p*e),
                     Jp   = gsl_sf_bessel_Jn(p, p*e),
                     Jpp1 = gsl_sf_bessel_Jn(p+1, p*e),
                     Jpp2 = gsl_sf_bessel_Jn(p+2, p*e);

        const double spl = sin(p*l), cpl = cos(p*l);

        P += OTS*(Jpm2 + Jpp2 - 2*Jp)*cpl;
        Q += -(Jpm2 - Jpp2 - 2*e*(Jpm1 - Jpp1) +(2/p)*Jp)*spl;
        R += (2/p)*Jp*spl;
    }

    /*const auto P = ((e + (-2 + e*e)*cu)*su)/(1 - ecu),
               Q = (OTS*(ecu - c2u))/(1 - ecu),
               R = esu;*/
    
    const auto A0 = R,
               A1 = -P*s2w + Q*c2w,
               A2 = P*c2w + Q*s2w;
    
    const auto a0 = si*si,
               a1 = 1+ci*ci,
               a2 = 2*ci;
    
    const auto sA = (H0/n) * ( a1*A1 + a0*A0 ),
               sB = (H0/n) * ( a2*A2 );
    
    const auto sp = c2psi*sA - s2psi*sB,
               sx = s2psi*sA + c2psi*sB;

    return std::make_tuple(sp, sx);
}

auto eccentric_residuals_fn_pt_PM(const BinaryMass &bin_mass,
                                 const BinaryState &bin_init,
                                 const double Fp, const double Fx, const double DGW,
                                 const EvolveCoeffs& ev_coeffs,
                                 const double t) {
    const auto [sp, sx] = eccentric_residuals_px_fn_pt_PM(bin_mass, bin_init, DGW, ev_coeffs, t);
    return Fp*sp + Fx*sx;
}

Signal1D eccentric_residuals_PM(const BinaryMass &bin_mass,
                               const BinaryState &bin_init,
                               const SkyPosition &bin_pos,
                               const SkyPosition &psr_pos,
                               const ResidualsTerms residuals_terms,
                               const Signal1D &ts){
    
    const auto [cosmu, Fp, Fx] = antenna_pattern(bin_pos, psr_pos);
    const auto DGW = bin_pos.DL;
    const auto ev_coeffs = compute_evolve_coeffs(bin_mass, bin_init);

    const auto nts = ts.size();
    Signal1D Rs(nts);
    for(unsigned idx=0; idx<nts; idx++){
        const auto& t = ts[idx];

        if(residuals_terms == ResidualsTerms::Earth){
            Rs[idx] = -eccentric_residuals_fn_pt_PM(bin_mass, bin_init, Fp, Fx, DGW, ev_coeffs, t);
        }
        else{
            const auto delay = -psr_pos.DL*(1-cosmu) / (1+bin_pos.z);
            if(residuals_terms == ResidualsTerms::Pulsar){
                Rs[idx] = eccentric_residuals_fn_pt_PM(bin_mass, bin_init, Fp, Fx, DGW, ev_coeffs, t + delay);
            }
            else{
                Rs[idx] =   eccentric_residuals_fn_pt_PM(bin_mass, bin_init, Fp, Fx, DGW, ev_coeffs, t + delay) 
                          - eccentric_residuals_fn_pt_PM(bin_mass, bin_init, Fp, Fx, DGW, ev_coeffs, t);
            }
        }
    }

    return Rs;
}

std::tuple<Signal1D, Signal1D> eccentric_residuals_px_PM(const BinaryMass &bin_mass,
                                                        const BinaryState &bin_init,
                                                        const double DGW,
                                                        const Signal1D &ts){

    const auto ev_coeffs = compute_evolve_coeffs(bin_mass, bin_init);

    const auto nts = ts.size();
    Signal1D Rps(nts), Rxs(nts);
    for(unsigned idx=0; idx<nts; idx++){
        const auto [spE, sxE] = eccentric_residuals_px_fn_pt_PM(bin_mass, bin_init, DGW, ev_coeffs, ts[idx]);
        Rps[idx] = spE;
        Rxs[idx] = sxE;
    }

    return std::make_tuple(Rps, Rxs);
}