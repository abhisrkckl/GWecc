#include "EccentricResiduals.hpp"
#include "OrbitalEvolution.hpp"
#include "AntennaPattern.hpp"
#include "NumericalWaveform.hpp"
#include "PN.hpp"
#include "mikkola.h"
#include "ipow.hpp"

template<int i>
double FeStatFunc_fn_pt(double t, void *_params){

    const auto &wf_params = *static_cast<WaveformParams*>(_params);
   
    const auto &bin_mass = wf_params.bin_mass;
    const auto &bin_init = wf_params.bin_init;
    const auto &ev_coeffs = wf_params.ev_coeffs;
    
    const auto bin_now = solve_orbit_equations(bin_init, ev_coeffs, t-bin_init.t);

    if(bin_now.merged){
        return 0;
    }
    
    const auto //&n     = bin_now.n,
               &e       = bin_now.e,
               &g       = bin_now.gamma,
               &l       = bin_now.l,
                    
               u        = MIKKOLA(l,e),
               L        = l+g,
               
               su       = sin(u),
               cu       = cos(u),
                       
               xi       = e*su,
               chi      = e*cu,
               
               k        = advance_of_periastron(bin_mass, bin_now),
               ephi     = angular_eccentricity(bin_mass, bin_now),
               
               betaphi  = (ephi>1e-15)  ? (1-sqrt(1-ephi*ephi))/ephi 
                                        : (e/2. + ipow(e,3)/8. + ipow(e,5)/16.),
               v_u      = 2*atan2(betaphi*su, 1-betaphi*cu),
               v_l      = v_u + xi,
                    
               W        = (1+k)*v_l,
               phi      = L+W,
               
               s2phi    = sin(2*phi),
               c2phi    = cos(2*phi),
                    
               OTS      = sqrt(1-e*e);
    
    double fi;
    switch(i){
        case 1:
            fi = 2*OTS*xi/(1-chi)/(1-chi)*c2phi;
            break;
        case 2:
            fi = 2*OTS*xi/(1-chi)/(1-chi)*s2phi;
            break;
        case 3:
            fi = (chi*(1-chi) - 2*OTS*OTS)/(1-chi)/(1-chi)*c2phi;
            break;
        case 4:
            fi = (chi*(1-chi) - 2*OTS*OTS)/(1-chi)/(1-chi)*s2phi;;
            break;
        case 5:
            fi = chi/(1-chi);
            break;
    }
    
    return fi;
}

/*
std::array<Signal1D, 5> FeStatFuncs_h(const BinaryMass &bin_mass,
                                      const BinaryState &bin_init,
                                      const SkyPosition &bin_pos,
                                      const SkyPosition &psr_pos,
                                      const Signal1D &ts){

    double cosmu, Fp, Fx;
    std::tie(cosmu, Fp, Fx) = AntennaPattern(bin_pos, psr_pos);

    //const auto integrator_size = 100;
    //const auto integ_eps_abs = 1e-6,
    //           integ_eps_rel = 1e-6;
    //const GSL_QAG_Integrator qag_integrator(integrator_size, integ_eps_abs, integ_eps_rel, GSL_INTEG_GAUSS15);

    const EvolveCoeffs_t ev_coeffs = compute_evolve_coeffs(bin_mass, bin_init);

    WaveformParams params { Fp, Fx, 1,
                            bin_mass,
                            bin_init,
                            ev_coeffs };

    const size_t length = ts.size();
    Signal1D B1(length), B2(length), B3(length), B4(length), B5(length);

    for(size_t i=0; i<length; i++){
        B1[i] = FeStatFunc_fn_pt<1>(ts[i], &params);
        B2[i] = FeStatFunc_fn_pt<2>(ts[i], &params);
        B3[i] = FeStatFunc_fn_pt<3>(ts[i], &params);
        B4[i] = FeStatFunc_fn_pt<4>(ts[i], &params);
        B5[i] = FeStatFunc_fn_pt<5>(ts[i], &params);
    }

    return {B1, B2, B3, B4, B5};
}
*/

std::array<Signal1D, 10> FeStatFuncs(const BinaryMass &bin_mass,
                                     const BinaryState &bin_init,
                                     const SkyPosition &bin_pos,
                                     const SkyPosition &psr_pos,
                                     const Signal1D &ts){

    double cosmu, Fp, Fx;
    std::tie(cosmu, Fp, Fx) = AntennaPattern(bin_pos, psr_pos);

    const auto integrator_size = 100;
    const auto integ_eps_abs = 1e-6,
               integ_eps_rel = 1e-6;
    const GSL_QAG_Integrator qag_integrator(integrator_size, integ_eps_abs, integ_eps_rel, GSL_INTEG_GAUSS15);

    const EvolveCoeffs_t ev_coeffs = compute_evolve_coeffs(bin_mass, bin_init);

    WaveformParams params { Fp, Fx, 1,
                            bin_mass,
                            bin_init,
                            ev_coeffs };

    gsl_function EccentricResiduals_gsl_func_1 {&FeStatFunc_fn_pt<1>, &params},
                 EccentricResiduals_gsl_func_2 {&FeStatFunc_fn_pt<2>, &params},
                 EccentricResiduals_gsl_func_3 {&FeStatFunc_fn_pt<3>, &params},
                 EccentricResiduals_gsl_func_4 {&FeStatFunc_fn_pt<4>, &params},
                 EccentricResiduals_gsl_func_5 {&FeStatFunc_fn_pt<5>, &params};

    Signal1D f1 = qag_integrator.eval_noerr(EccentricResiduals_gsl_func_1, ts[0], ts),
             f2 = qag_integrator.eval_noerr(EccentricResiduals_gsl_func_2, ts[0], ts),
             f3 = qag_integrator.eval_noerr(EccentricResiduals_gsl_func_3, ts[0], ts),
             f4 = qag_integrator.eval_noerr(EccentricResiduals_gsl_func_4, ts[0], ts),
             f5 = qag_integrator.eval_noerr(EccentricResiduals_gsl_func_5, ts[0], ts);
    
    return {Fp*f1, Fx*f1, Fp*f2, Fx*f2, Fp*f3, Fx*f3, Fp*f4, Fx*f4, Fp*f5, Fx*f5};

}
