#ifndef _OrbitalEvolution_hpp_
#define _OrbitalEvolution_hpp_ 1

#include "Binary.hpp"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <vector>
#include <tuple>

/*
 * [e, tau, lbar, gammabar]
 */
typedef std::tuple<std::vector<double>, std::vector<double> > OrbitIntegrals;

OrbitIntegrals precompute_orbit(const double emin, const double taumax);

/****************/

struct EvolveCoeffs_t{

    double  A, AT, AG,
            tau0, P,
            lbar0, alpha,
            gbar0, beta,
            gbar20, beta2,
            gbar30, beta3,
            eta,
            sini, cosi,
            sin2Omega, cos2Omega;
};

/****************/

class Evolve{

private:
    // These acceleration objects are only used for single point functions. 
    gsl_interp_accel *acc, *acci;                    
    gsl_spline       *spline, *splinei;

    double tau_max, e_max,
           tau_min, e_min;
                
    //Coefficients for evaluating tau->inf asymptotic expressions
    const double a  = 2*sqrt(2.)/5/pow(5,  63./2299)/pow(17,1181./2299);
    double b;

    void initialize(const std::vector<double> &taus_pre, const std::vector<double> &es_pre);
    
    // Don't allow instantiation
    Evolve(const char *datafilename);
    Evolve();
        
public:   
    // Don't allow copying
    Evolve(Evolve const&)          = delete;
    void operator=(Evolve const&)  = delete;
    
    // This class is a singleton.
    static const Evolve& instance(){
        static const Evolve ev("_evolve_.dat");
        return ev;
    }
    
    size_t rows;
    
    double e_from_tau(const double tau) const;
    double tau_from_e(const double e) const;
    
    BinaryState solve_orbit_equations(const BinaryMass &bin_mass, 
                                      const BinaryState &bin_init, 
                                      const double delay) const;
    
    BinaryState solve_orbit_equations(const BinaryState &bin_init,
                                      const EvolveCoeffs_t &ev_coeffs, 
                                      const double delay) const;
    
    double phase_err(const BinaryMass &bin_mass, 
                     const BinaryState &bin_init, 
                     const double delay) const;
    
    ~Evolve();
};

double compute_P_coeff(const double Mchirp, const double n0, const double e0);
double compute_alpha_coeff(const double Mchirp, const double n0, const double e0);
double compute_beta_coeff(const double Mchirp, const double M, const double n0, const double e0);
double compute_beta2_coeff(const double Mchirp, const double M, const double n0, const double e0);
double compute_beta3_coeff(const double Mchirp, const double M, const double n0, const double e0);

double n_from_e(const double n0, const double e0, const double e);
double lbar_from_e(const double e);
double gbar_from_e(const double e);
double gbar2_from_e(const double e, const double eta);
double gbar3_from_e(const double e, const double eta);

BinaryState solve_orbit_equations(const BinaryMass &bin_mass, const BinaryState &bin_init, const double delay);

EvolveCoeffs_t compute_evolve_coeffs(const BinaryMass &bin_mass, const BinaryState &bin_init);
BinaryState solve_orbit_equations(const BinaryState &bin_init, const EvolveCoeffs_t &ev_coeffs, const double delay);

double phase_err(const BinaryMass &bin_mass, const BinaryState &bin_init, const double delay);

double compute_phi(BinaryMass &bin_mass, BinaryState &bin_state);

#endif
