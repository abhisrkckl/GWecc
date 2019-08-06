
#include "OrbitalEvolution.hpp"
#include <cmath>
#include <gsl/gsl_sf_hyperg.h>
#include "ipow.hpp"
#include "PN.hpp"
#include "mikkola.h"

BinaryState solve_orbit_equations(const BinaryMass &bin_mass, const BinaryState &bin_init, const double delay){
    return Evolve::instance().solve_orbit_equations(bin_mass, bin_init, delay);    
}

BinaryState solve_orbit_equations(const BinaryState &bin_init, const EvolveCoeffs_t &ev_coeffs, const double delay){
    return Evolve::instance().solve_orbit_equations(bin_init, ev_coeffs, delay);    
}

double phase_err(const BinaryMass &bin_mass, const BinaryState &bin_init, const double delay){
    return Evolve::instance().phase_err(bin_mass, bin_init, delay);
}

double compute_P_coeff(const double Mchirp, const double n0, const double e0){

    double A   = pow(Mchirp,5./3)/5;
    double e02 = e0*e0;    
    
    return A/3  * pow(n0,8./3)
                * pow(e0,48./19)
                * pow(304 + 121*e02, 3480./2299)
                / pow(1-e02,4);
}

double compute_alpha_coeff(const double Mchirp, const double n0, const double e0){
    
    double A = pow(Mchirp,5./3)/5;
    double e02 = e0*e0;
    
    return 3./A / pow(n0,5./3)
                / pow(e0,30./19)
                / pow(304 + 121*e02, 2175./2299)
                * pow(1-e02,2.5);
    
}

double compute_beta_coeff(const double Mchirp, const double M, const double n0, const double e0){
    
    double A = pow(Mchirp,5./3)/5;
    double e02 = e0*e0;
    
    return 9./A * pow(M, 2./3)
                / n0
                / pow(e0,18./19)
                / pow(304 + 121*e02, 1305./2299)
                * pow(1-e02,1.5);
    
}

double compute_beta2_coeff(const double Mchirp, const double M, const double n0, const double e0){
    
    double A = pow(Mchirp,5./3)/5;
    double e02 = e0*e0;
    
    return 3./(4*A) * pow(M, 4./3)
                    / cbrt(n0)
                    / pow(e0,6./19)
                    / pow(304 + 121*e02, 435./2299)
                    * sqrt(1-e02);
    
}

double n_from_e(const double n0, const double e0, const double e){
    double e02 = e0*e0,
           e2  = e*e; 
    return n0   * pow(e0/e,18./19)
                * sqrt(ipow((1-e2)/(1-e02),3))
                * pow((304 + 121*e02)/(304 + 121*e2),1305./2299);
}

double lbar_from_e(const double e){
    const/*expr*/ double coeff = pow(19,2175./2299)/30/pow(2, 496./2299);
    return  coeff * e * pow(e,11./19) * gsl_sf_hyperg_2F1(124./2299, 15./19, 34./19, -121.*e*e/304);
}

double gbar_from_e(const double e){
    const/*expr*/ double coeff = pow(19,1305./2299)/36/pow(2,1677./2299);
    return  coeff * pow(e,18./19) * gsl_sf_hyperg_2F1(994./2299,  9./19, 28./19, -121.*e*e/304);
}

double gbar2_from_e(const double e, const double eta){
    
    const/*expr*/ double coeff = 3 * pow(2,1740./2299) * pow(19,435./2299);
    
    double e2 = e*e;
    
    return pow(e,6./19)/336  * ( 4*(51-26*eta)*pow(304+121*e2,435./2299)
                                    + coeff*(23+ 2*eta)*gsl_sf_hyperg_2F1(3./19, 1864./2299, 22./19, -121*e2/304)
                               );
}

EvolveCoeffs_t compute_evolve_coeffs(const BinaryMass &bin_mass, const BinaryState &bin_init){

    static const Evolve &ev = Evolve::instance();
    
    const double  &n0    = bin_init.n,
                  &e0    = bin_init.e,
                  
                  Mchirp = bin_mass.chirp_mass(),
                  M      = bin_mass.mass(),
                  eta    = bin_mass.symmetric_mass_ratio();
    
    const double  A      = pow(Mchirp,5./3)/5,
                  AA     = 256.*A*pow(n0,8./3),
                  AG     = pow(M,2./3) / (32*A),
                  tau0   = ev.tau_from_e(e0),
                  P      = compute_P_coeff(Mchirp,n0,e0),
                  alpha  = compute_alpha_coeff(Mchirp, n0, e0),
                  lbar0  = lbar_from_e(e0),
                  beta   = compute_beta_coeff(Mchirp, M, n0, e0),
                  gbar0  = gbar_from_e(e0),
                  beta2  = compute_beta2_coeff(Mchirp, M, n0, e0),
                  gbar20 = gbar2_from_e(e0,eta),
                
                  sini   = sin(bin_init.i),
                  cosi   = cos(bin_init.i),
                
                  sin2Omega = sin(2*bin_init.Omega),
                  cos2Omega = cos(2*bin_init.Omega);

    return EvolveCoeffs_t {    A, AA, AG,
                                tau0, P,
                                lbar0, alpha,
                                gbar0, beta,
                                gbar20, beta2,
                                eta,
                                sini, cosi,
                                sin2Omega, cos2Omega  };    
    
}

double compute_phi(BinaryMass &bin_mass, BinaryState &bin_state){
    const auto e    = bin_state.e,
               g    = bin_state.gamma,
               l    = bin_state.l,
               
               u    = MIKKOLA(l,e),
               L    = l+g,
               
               su   = sin(u),
               cu   = cos(u),
                       
               //esu    = e*su,
               //ecu    = e*cu,
               
               k    = advance_of_periastron(bin_mass, bin_state),
               ephi = angular_eccentricity(bin_mass, bin_state),
               betaphi = (ephi>1e-15) ? (1-sqrt(1-ephi*ephi))/ephi 
                                      : (e/2. + ipow(e,3)/8. + ipow(e,5)/16.),
                               
               v_u = 2*atan2(betaphi*su, 1-betaphi*cu),
               v_l = v_u + e*su,
               
               W   = (1+k)*v_l,
               phi = L+W;
           
    return phi;
}
