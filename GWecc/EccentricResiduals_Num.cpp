#include "EccentricResiduals.hpp"
#include "NumericalWaveform.hpp"
#include "AntennaPattern.hpp"
#include "OrbitalEvolution.hpp"
#include "PN.hpp"
#include "mikkola.h"
#include "ipow.hpp"
#include <iostream>

double EccentricResiduals_fn_pt(double t, void *_params);

Signal1D EccentricResiduals_Num(const BinaryMass &bin_mass,
                                const BinaryState &bin_init,
                                const SkyPosition &bin_pos,
                                const SkyPosition &psr_pos,
                                const ResidualsTerms residuals_terms,
                                const Signal1D &ts){
    
    const auto [cosmu, Fp, Fx] = AntennaPattern(bin_pos, psr_pos);
    
    if(residuals_terms==ResidualsTerms::Earth){
        return -EccentricResiduals_fn_Num(bin_mass, bin_init, Fp, Fx, bin_pos.DL, ts);
    }
    else if(residuals_terms==ResidualsTerms::Pulsar){
        
        const auto delay = -psr_pos.DL*(1-cosmu) / (1+bin_pos.z);
        //const auto bin_psrterm = solve_orbit_equations(bin_mass, bin_init, delay);
        
        return EccentricResiduals_fn_Num(bin_mass, /*bin_psrterm*/bin_init, Fp, Fx, bin_pos.DL, ts + delay);
    }
    else{
        const auto delay = -psr_pos.DL*(1-cosmu);
        //const auto bin_psrterm = solve_orbit_equations(bin_mass, bin_init, delay);
        
        return    EccentricResiduals_fn_Num(bin_mass, /*bin_psrterm*/bin_init, Fp, Fx, bin_pos.DL, ts + delay)
                - EccentricResiduals_fn_Num(bin_mass, bin_init,    Fp, Fx, bin_pos.DL, ts);
    }
}

Signal1D EccentricResiduals_fn_Num(const BinaryMass &bin_mass,
                                   const BinaryState &bin_init,
                                   const double Fp, const double Fx, const double DGW,
                                   const Signal1D &ts){
    
    const auto integrator_size = 100;
    const auto integ_eps_abs = 1e-12,
               integ_eps_rel = 1e-12;
    const GSL_QAG_Integrator qag_integrator(integrator_size, integ_eps_abs, integ_eps_rel, GSL_INTEG_GAUSS15); 
    
    const EvolveCoeffs_t ev_coeffs = compute_evolve_coeffs(bin_mass, bin_init);
    
    WaveformParams params { Fp, Fx, DGW,
                            bin_mass,
                            bin_init,
                            ev_coeffs };
    
    
    gsl_function EccentricResiduals_gsl_func {&EccentricResiduals_fn_pt, &params};
    
    return qag_integrator.eval_noerr(EccentricResiduals_gsl_func, ts[0], ts);
    
    /*
    const int length = ts.size();
    Signal1D result(length);
    for(int i=0;i<length;i++){
        result[i] = EccentricResiduals_fn_pt(ts[i], &params);
    }
    
    return result;
    */
    
}

double EccentricResiduals_fn_pt(double t, void *_params){

    const auto [Fp, Fx, DGW,
                bin_mass, 
                bin_init,
                ev_coeffs ] = *reinterpret_cast<WaveformParams*>(_params);
            
    const auto bin_now = solve_orbit_equations(bin_init, ev_coeffs, t-bin_init.t);
    
    const auto //&n     = bin_now.n,
               &e       = bin_now.e,
               &g       = bin_now.gamma,
               &l       = bin_now.l,
                    
               u        = MIKKOLA(l,e),
               L        = l+g,
               
               su       = sin(u),
               cu       = cos(u),
                       
               esu      = e*su,
               ecu      = e*cu,
               
               x        = PN_param_x(bin_mass, bin_now),
               sqrtx    = sqrt(x),
               k        = advance_of_periastron(bin_mass, bin_now),
               ephi     = angular_eccentricity(bin_mass, bin_now),
               
               betaphi  = (ephi>1e-15)  ? (1-sqrt(1-ephi*ephi))/ephi 
                                        : (e/2. + ipow(e,3)/8. + ipow(e,5)/16.),
               v_u      = 2*atan2(betaphi*su, 1-betaphi*cu),
               v_l      = v_u + esu,
                    
               W        = (1+k)*v_l,
               phi      = L+W,
               
               si       = ev_coeffs.sini,    //sin(bin_now.i),
               ci       = ev_coeffs.cosi,    //cos(bin_now.i),
               
               delta    = bin_mass.differential_mass_ratio(),
               
               s1phi    = sin(phi),
               c1phi    = cos(phi),
               s2phi    = 2*s1phi*c1phi,
               c2phi    = c1phi*c1phi-s1phi*s1phi,
               s3phi    = s1phi*(-s1phi*s1phi+3*c1phi*c1phi),
               c3phi    = c1phi*( c1phi*c1phi-3*s1phi*s1phi),
                    
               OTS      = sqrt(1-e*e),
               
               c2Om     = ev_coeffs.cos2Omega,    //cos(2*bin_now.Omega),
               s2Om     = ev_coeffs.sin2Omega,    //sin(2*bin_now.Omega),
               
               H0       = (DGW==-1) ?1 :GWAmplitude(bin_mass, bin_now, DGW);    
    
    const auto h_mq_A   = 1./ipow(1-ecu,2) * (    - (ci*ci+1)*(2*OTS*esu)         *s2phi
                            + (ci*ci+1)*(2*e*e - ecu*ecu + ecu - 2)    *c2phi
                            + (si*si)  *(1-ecu)*ecu          ),
               h_mq_B   = 1./ipow(1-ecu,2) * (             (2*OTS*esu)        *c2phi
                                   +         (2*e*e - ecu*ecu + ecu - 2) *s2phi
                                    )*2*ci;
    
    const auto h_cq_A   = 1./ipow(1-ecu,3) * (   (1-ecu)*(-5-ipow(ci,2) + (-2+6*ipow(ci,2))*ecu)*OTS        *c1phi
                          +  2*(1-3*ipow(ci,2))*ipow(1-ecu,2)*esu            *s1phi 
                          + (1+ipow(ci,2))*OTS*(9-8*ipow(e,2)-7*ecu+6*ipow(ecu,2))    *c3phi 
                          + (1+ipow(ci,2))*2*esu*(5-4*ipow(e,2)-2*ecu+ipow(ecu,2))    *s3phi
                          )*delta*si/4,
               h_cq_B   = 1./ipow(1-ecu,3) * (   2*ipow(1-ecu,2)*esu                        *c1phi 
                                 + (1-ecu)*(-3+2*ecu)*OTS                    *s1phi           
                                 + 2*esu*(-5+4*ipow(e,2)+2*ecu-ipow(ecu,2))             *c3phi
                                 + (9-8*ipow(e,2)-7*ecu+6*ipow(ecu,2))*OTS            *s3phi
                               )*delta*ci*si/2;
    
    const auto hA       = h_mq_A + sqrtx*h_cq_A,
               hB       = h_mq_B + sqrtx*h_cq_B;
    
    const auto hp       = c2Om*hA - s2Om*hB,
               hx       = c2Om*hB + s2Om*hA;
    
    return H0*(Fp*hp + Fx*hx);
}

std::tuple<Signal1D, Signal1D> EccentricResiduals_px_Num(const BinaryMass &bin_mass,
                                                         const BinaryState &bin_init,
                                                         const double DGW,
                                                         const Signal1D &ts){

        const auto integrator_size = 100;
        const auto integ_eps_abs = 1e-6,
                   integ_eps_rel = 1e-6;
        const GSL_QAG_Integrator qag_integrator(integrator_size, integ_eps_abs, integ_eps_rel, GSL_INTEG_GAUSS15);

        const EvolveCoeffs_t ev_coeffs = compute_evolve_coeffs(bin_mass, bin_init);

        WaveformParams params_p { 1, 0, DGW,
                                  bin_mass,
                                  bin_init,
                                  ev_coeffs },
                       params_x { 0, 1, DGW,
                                  bin_mass,
                                  bin_init,
                                  ev_coeffs };


        gsl_function EccentricResiduals_gsl_func_p {&EccentricResiduals_fn_pt, &params_p},
                     EccentricResiduals_gsl_func_x {&EccentricResiduals_fn_pt, &params_x};

        Signal1D Rp = qag_integrator.eval_noerr(EccentricResiduals_gsl_func_p, ts[0], ts),
                 Rx = qag_integrator.eval_noerr(EccentricResiduals_gsl_func_x, ts[0], ts);
        
        return std::make_tuple(Rp,Rx);

}

std::tuple<Signal1D, Signal1D> EccentricWaveform_px(const BinaryMass &bin_mass,
                                                    const BinaryState &bin_init,
                                                    const double DGW,
                                                    const Signal1D &ts){
    const size_t length = ts.size();                                                
    Signal1D hp(length), hx(length);
    
    const EvolveCoeffs_t ev_coeffs = compute_evolve_coeffs(bin_mass, bin_init);
    
    WaveformParams params_p { 1, 0, DGW,
                              bin_mass,
                              bin_init,
                              ev_coeffs },
                   params_x { 0, 1, DGW,
                              bin_mass,
                              bin_init,
                              ev_coeffs };
    
    for(unsigned i=0; i<length; i++){
        hp[i] = EccentricResiduals_fn_pt(ts[i], &params_p);
        hx[i] = EccentricResiduals_fn_pt(ts[i], &params_x);    
    }
    
    return std::make_tuple(hp, hx);
}

