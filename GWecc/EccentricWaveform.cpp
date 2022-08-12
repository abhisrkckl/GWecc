#include "GWecc.hpp"
#include "mikkola.h"
#include "PN.hpp"
#include "ipow.hpp"
#include "EccentricResiduals.hpp"
#include "NumericalWaveform.hpp"
#include "AntennaPattern.hpp"

std::vector<double> EccentricWaveform_fn(const double M, const double q,
                                         const double Omega, const double i,
                                         const double n, const double e, const double l, const double gamma,
                                         const double DGW){

    const BinaryMass bin_mass{M, q};
    const BinaryState bin_now { 0,
                                Omega, i,
                                n, e, l, gamma};
    
    const double  u        = MIKKOLA(l,e),
               	  L        = l+gamma,
               
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
                  
                  si       = sin(i),    //sin(bin_now.i),
                  ci       = cos(i),    //cos(bin_now.i),
                  
                  delta    = fabs(1-q)/(1+q),
                  
                  s1phi    = sin(phi),
                  c1phi    = cos(phi),
                  s2phi    = 2*s1phi*c1phi,
                  c2phi    = c1phi*c1phi-s1phi*s1phi,
                  s3phi    = s1phi*(-s1phi*s1phi+3*c1phi*c1phi),
                  c3phi    = c1phi*( c1phi*c1phi-3*s1phi*s1phi),
                        
                  OTS      = sqrt(1-e*e),
                  
                  c2Om     = cos(2*Omega),    //cos(2*bin_now.Omega),
                  s2Om     = sin(2*Omega),    //sin(2*bin_now.Omega),
                  
                  H0       = (DGW==-1) ?1 :GWAmplitude(bin_mass, bin_now, DGW);
    
    const auto h_mq_A   = 1./ipow(1-ecu,2) * (    - (ci*ci+1)*(2*OTS*esu)         *s2phi
                            + (ci*ci+1)*(2*e*e - ecu*ecu + ecu - 2)    *c2phi
                            + (si*si)  *(1-ecu)*ecu          ),
               h_mq_B   = 1./ipow(1-ecu,2) * (             (2*OTS*esu)        *c2phi
                            +           (2*e*e - ecu*ecu + ecu - 2) *s2phi
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
    
    return {H0*hp, H0*hx};
 
}

Signal1D EccentricWaveform( const BinaryMass &bin_mass,
                            const BinaryState &bin_init,
                            const SkyPosition &bin_pos,
                            const SkyPosition &psr_pos,
                            const ResidualsTerms residuals_terms,
                            const Signal1D &ts){
    
    const auto [cosmu, Fp, Fx] = AntennaPattern(bin_pos, psr_pos);

    const size_t length = ts.size();                                                
    Signal1D h(length);
    
    const EvolveCoeffs_t ev_coeffs = compute_evolve_coeffs(bin_mass, bin_init);
    
    WaveformParams params { Fp, Fx, bin_pos.DL,
                            bin_mass,
                            bin_init,
                            ev_coeffs };
    
    if(residuals_terms==ResidualsTerms::Earth){
        for(unsigned i=0; i<length; i++){
            h[i] = -EccentricResiduals_fn_pt(ts[i], &params);
        }
    }
    else if(residuals_terms==ResidualsTerms::Pulsar){
        const auto delay = -psr_pos.DL*(1-cosmu) / (1+bin_pos.z);
        for(unsigned i=0; i<length; i++){
            h[i] = EccentricResiduals_fn_pt(ts[i] + delay, &params);
        }
    }
    else{
        const auto delay = -psr_pos.DL*(1-cosmu) / (1+bin_pos.z);
        for(unsigned i=0; i<length; i++){
            h[i] = EccentricResiduals_fn_pt(ts[i] + delay, &params) - EccentricResiduals_fn_pt(ts[i], &params);
        }
    }

    return h;
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