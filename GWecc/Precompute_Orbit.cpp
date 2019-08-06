#include <stdio.h>
#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <vector>
#include "OrbitalEvolution.hpp"

/*
 * Function containing expressions for the RHS of the system of ODEs.
 * 
 * tau      is the independent variable
 * y[]      are the dependent variables :    y=[e]
 * f[]      are the derivatives
 * params   are parameters appearing in the ODEs    :    params=[]
 */
int ODE_RHS(double tau, const double y[], double f[], void *params) {

    double e   = y[0],
           e2  = e*e;

    double de_dtau = pow(1-e2,1.5) / pow(e,29./19) / pow(121*e2+304,1181./2299);

    f[0] = de_dtau;

    return GSL_SUCCESS;
}

OrbitIntegrals precompute_orbit(const double emin, const double taumax) {

    std::vector<double> taus;
    std::vector<double> es;

    // Initial conditions
    // The integration cannot start from zero as the derivative is singular there.
    const double  taumin = 19*pow(19,1181./2299)*pow(emin,48./19)/(6*pow(2,2173./2299));

    // Start the integration here
    double tau   = taumin;
    double e; 
    double y[1]  = {emin};

    // Parametars appearing in the system of ODEs 
    // double params[] = {};

    // Variables for solving system of ODEs
    constexpr size_t neqs = 1;                            // Number of ODEs
    constexpr double eps_abs = 1e-20, eps_rel = 1e-22;    // desired tolerance
    double stepsize = 1e-20;                              // initial integration step size.
    int status;

    // output_length depends on the tolerances eps_abs and eps_rel.
    // If the tolerances are changed, output_length should also me changed accordingly.
    constexpr size_t output_length = 29788;
    taus.reserve(output_length);
    es.reserve(output_length);

    /*
     * Explicit embedded Runge-Kutta-Fehlberg (4,5) method. 
     * This method is a good general-purpose integrator.
     */
    gsl_odeiv2_step    *s = gsl_odeiv2_step_alloc(gsl_odeiv2_step_rkf45, neqs);
    gsl_odeiv2_control *c = gsl_odeiv2_control_y_new(eps_abs, eps_rel);
    gsl_odeiv2_evolve *ev = gsl_odeiv2_evolve_alloc(neqs);
    gsl_odeiv2_system sys = {ODE_RHS, NULL, neqs, NULL};

    taus.push_back(taumin);
    es.push_back(emin);
    
    // Evolution loop 
    while(tau<taumax){

        // y = [e]
        status = gsl_odeiv2_evolve_apply(ev, c, s, &sys, &tau, taumax, &stepsize, y);

        if(status != GSL_SUCCESS) {
            fprintf(stderr, "Stopping...\n");
            break;
        }

        e = y[0];

        taus.push_back(tau);
        es.push_back(e);
    }

    gsl_odeiv2_evolve_free(ev);
    gsl_odeiv2_control_free(c);
    gsl_odeiv2_step_free(s);
    
    return std::make_tuple(taus, es);
}

