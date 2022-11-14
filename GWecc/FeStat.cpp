#include "EccentricResiduals.hpp"
#include "OrbitalEvolution.hpp"
#include "AntennaPattern.hpp"
#include "NumericalWaveform.hpp"
#include "PN.hpp"
#include "mikkola.h"
#include "ipow.hpp"
#include "FeStat.hpp"
#include <sstream>

template<unsigned i>
double FeStat_A_fn_pt(double t, void *_params){

    static_assert(i>0 && i<4, "Invalid A function index.");

    const auto &wf_params = *static_cast<WaveformParams*>(_params);
   
    const auto &bin_mass = wf_params.bin_mass;
    const auto &bin_init = wf_params.bin_init;
    const auto &ev_coeffs = wf_params.ev_coeffs;
    
    const auto bin_now = solve_orbit_equations(bin_init, ev_coeffs, t-bin_init.t);

    if(bin_now.merged){
        return 0;
    }
    
    const auto &n       = bin_now.n,
               &e       = bin_now.e,
               &g       = bin_now.gamma,
               &l       = bin_now.l,
                    
               u        = mikkola(l,e),
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
    
    const double N = 2*OTS*xi/(1-chi)/(1-chi),
                 P = (chi*(1-chi) - 2*OTS*OTS)/(1-chi)/(1-chi),
                 Q = chi/(1-chi);

    double result;

    if constexpr(i==1){
        result = N*c2phi + P*s2phi;
    }
    else if constexpr(i==2){
        result = N*s2phi - P*c2phi;
    }
    else if constexpr(i==3){
        result = Q;
    }

    if(std::isnan(result)){
        std::ostringstream errstrm;
        errstrm<<"NaN encountered in FeStat_A_fn_pt<"<<i<<">. "
               <<"M = "<<bin_mass.mass()<<", "
               <<"eta = "<<bin_mass.symmetric_mass_ratio()<<", "
               <<"e = "<<e<<", "
               <<"n = "<<n<<", "
               <<"gamma = "<<g<<", "
               <<"l = "<<l<<", "
               <<"u = "<<u<<", "
               <<"L = "<<L<<", "
               <<"su = "<<su<<", "
               <<"cu = "<<cu<<", "
               <<"xi = "<<xi<<", "
               <<"chi = "<<chi<<", "
               <<"k = "<<k<<", "
               <<"ephi = "<<ephi<<", "
               <<"betaphi = "<<betaphi<<", "
               <<"v_u = "<<v_u<<", "
               <<"v_l = "<<v_l<<", "
               <<"W ="<<W<<", "
               <<"phi = "<<phi<<", "
               <<"s2phi = "<<s2phi<<", "
               <<"c2phi = "<<c2phi<<", "
               <<"OTS = "<<OTS<<", ";
        throw std::runtime_error(errstrm.str());
    }

    return result;
}

std::array<Signal1D, 6> FeStatFuncs(const BinaryMass &bin_mass,
                                    const BinaryState &bin_init,
                                    const SkyPosition &bin_pos,
                                    const SkyPosition &psr_pos,
                                    const Signal1D &ts){

    const auto [cosmu, Fp, Fx] = antenna_pattern(bin_pos, psr_pos);

    gsl_set_error_handler_off();

    const auto integrator_size = 100;
    const auto integ_eps_abs = 1e-6,
               integ_eps_rel = 1e-6;
    const GSL_CQUAD_Integrator cquad_integrator(integrator_size, integ_eps_abs, integ_eps_rel, GSL_INTEG_GAUSS15);

    const EvolveCoeffs_t ev_coeffs = compute_evolve_coeffs(bin_mass, bin_init);

    WaveformParams params { Fp, Fx, 1,
                            bin_mass,
                            bin_init,
                            ev_coeffs };

    gsl_function EccentricResiduals_gsl_func_1 {&FeStat_A_fn_pt<1>, &params},
                 EccentricResiduals_gsl_func_2 {&FeStat_A_fn_pt<2>, &params},
                 EccentricResiduals_gsl_func_3 {&FeStat_A_fn_pt<3>, &params};

    Signal1D f1 = cquad_integrator.eval_noerr(EccentricResiduals_gsl_func_1, ts[0], ts),
             f2 = cquad_integrator.eval_noerr(EccentricResiduals_gsl_func_2, ts[0], ts),
             f3 = cquad_integrator.eval_noerr(EccentricResiduals_gsl_func_3, ts[0], ts);
    
    const double H0 = GWAmplitude(bin_mass, bin_init, bin_pos.DL);

    return {H0*Fp*f1, H0*Fp*f2, H0*Fp*f3,
            H0*Fx*f1, H0*Fx*f2, H0*Fx*f3};

}
