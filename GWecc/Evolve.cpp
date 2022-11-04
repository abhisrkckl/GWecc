#include "OrbitalEvolution.hpp"
#include "PN.hpp"
#include <cstdio>
#include <stdexcept>


using std::runtime_error, std::nullopt;


void write_data_file(const char *datafilename, 
                     const vector<double> &taus, 
                     const vector<double> &es){
    const size_t rows = taus.size();

    FILE *outfile = fopen(datafilename, "wb");
    fwrite(&rows,    sizeof(rows),   1,    outfile);
    fwrite(&taus[0], sizeof(double), rows, outfile);
    fwrite(&es[0],   sizeof(double), rows, outfile);
    fclose(outfile);

    printf("Successfully wrote %s.\n",datafilename);
}


Evolve::Evolve(const char *datafilename){

    FILE *datafile = fopen(datafilename,"rb");
    if(datafile==NULL){
        
        fprintf(stderr, "ERROR : Unable to read pre-computed integral from %s. Re-computing integral...\n", datafilename);
        
        //Evolve();
        constexpr double emin   = 2.5e-9,  // Initial condition
                         taumax = 1000;    // Stopping condition

        const auto [taus, es] = precompute_orbit(emin,taumax);

        initialize(taus, es);

        write_data_file(datafilename, taus, es);
    }
    else{
        size_t _rows;
    
        fseek(datafile, 0, SEEK_END);
        size_t filesize = ftell(datafile);
        rewind(datafile);

        fread(&_rows,sizeof(_rows),1,datafile);

        size_t filesize_expected = sizeof(_rows) + 2*_rows*sizeof(double);
        if(filesize_expected!=filesize){
            fprintf(stderr, "ERROR : Unable to read pre-computed integral from %s. Re-computing integral...\n", datafilename);
            
            // Initial condition
            constexpr double emin   = 2.5e-9,  
                             taumax = 1000;
            
            const auto [taus, es] = precompute_orbit(emin,taumax);

            initialize(taus, es);

            write_data_file(datafilename, taus, es);
        }
        else{

            vector<double> taus_pre(_rows);
            vector<double> es_pre(_rows);

            fread(taus_pre.data() ,sizeof(double),_rows,datafile);
            fread(es_pre.data()   ,sizeof(double),_rows,datafile);

            initialize(taus_pre, es_pre);

            printf("Successfully initiated from %s.\n", datafilename);
        }

        fclose(datafile);
    }
}

Evolve::Evolve(){

    constexpr double emin   = 2.5e-9,  // Initial condition
                     taumax = 1000;    // Stopping condition
    
    const auto [taus, es] = precompute_orbit(emin,taumax);
    
    initialize(taus, es);
}

Evolve::~Evolve(){
    gsl_spline_free(spline);
    gsl_spline_free(splinei);

    gsl_interp_accel_free(acc);
    gsl_interp_accel_free(acci);
}

void Evolve::initialize(const vector<double> &taus_pre, const vector<double> &es_pre){

    rows = taus_pre.size();

    tau_max = taus_pre[rows-1];
    tau_min = taus_pre[0];
    e_max   = es_pre[rows-1];
    e_min   = es_pre[0];

    // Spline for e(tau)
    acc     = gsl_interp_accel_alloc();
    spline  = gsl_spline_alloc(gsl_interp_cspline, rows);
    gsl_spline_init(spline, taus_pre.data(), es_pre.data(), rows);

    // Spline for tau(e)
    acci    = gsl_interp_accel_alloc();
    splinei = gsl_spline_alloc(gsl_interp_cspline, rows);
    gsl_spline_init(splinei, es_pre.data(), taus_pre.data(), rows);

    // Coefficients for asymptotic series
    b  = 2./sqrt(1-e_max) - a*tau_max;
}

double Evolve::e_from_tau(const double tau) const{

    double e;
    if(tau<tau_min){
        const double coeff = pow(2,559./726) * pow(3,19./48) / pow(19,145./242);
        e = coeff * pow(tau, 19./48);
    }
    else if(tau>tau_max){
        double atau_b = a*tau+b;
        e = 1 - 4./(atau_b*atau_b);
    }
    else{
        e = gsl_spline_eval(spline, tau, acc);
    }
    return e;
}

double Evolve::tau_from_e(const double e) const{

    double tau;
    if(e<e_min){
        const double coeff = 19*pow(19,1181./2299)/(6*pow(2,2173./2299));
        tau = coeff * pow(e,48./19);
    }
    else if(e>e_max){
        tau = (2./sqrt(1-e) - b)/a;
    }
    else{
        tau = gsl_spline_eval(splinei, e, acci);
    }
    return tau;
}

/* 
 * This function is *NOT* thread-safe.
 */ 
BinaryState Evolve::solve_orbit_equations(const BinaryMass &bin_mass, 
                                          const BinaryState &bin_init, 
                                          const optional<PulsarTermPhase> ptphase,
                                          const double delay) const{

    const double &t0    = bin_init.t,
                 &n0    = bin_init.n,
                 &e0    = bin_init.e,
                 &gamma0= bin_init.gamma,
                 &l0    = bin_init.l,
                 
                 t      = t0 + delay,
                  
                 eta    = bin_mass.symmetric_mass_ratio(),
                 M      = bin_mass.mass(),
                 Mchirp = bin_mass.chirp_mass();
    
    double n, e, gamma, l, gamma1, gamma2, gamma3;
    bool merged = false;

    if(e0==0){
        
        const double A     = pow(Mchirp,5./3)/5;
        const double T     = 256*A*pow(n0,8./3)*delay;

        if(T>1){
            merged = true;
            throw runtime_error("tau<0 found in solve_orbit_equations");
        }
        else{
            e     = e0;
            n     = n0 / pow(1-T, 3./8);
            
            if(ptphase){
                l = ptphase->lp;
                gamma = ptphase->gammap;
            }
            else{
                l     = l0 + (pow(n0,-5./3) - pow(n,-5./3))/(160*A);
                gamma = gamma0 + (1./n0 - 1./n) * pow(M,2./3) / (32*A);
            }
        }
    }
    else{
        const double P     = compute_P_coeff(Mchirp,n0,e0);
        const double tau0  = tau_from_e(e0);
        //fprintf(stderr, "P = %e, tau0 = %e, t = %e\n",P,tau0,t);
        double tau   = tau0 - P*delay;

        if(tau<0){
            //fprintf(stderr, "Error: The binary has already merged.\n");
            //e=n=l=gamma=NAN;
            merged = true;
            throw runtime_error("tau<0 found in solve_orbit_equations");
        }
        else{
            
            e = e_from_tau(tau);
            n = n_from_e(n0, e0, e);
            
            if(ptphase){
                l = ptphase->lp;
                gamma = ptphase->gammap;
            }
            else{
                const double alpha     = compute_alpha_coeff(Mchirp, n0, e0);
                const double lbar0     = lbar_from_e(e0);  
                double       lbar      = lbar_from_e(e);         
                l = l0 + (lbar0-lbar)*alpha;
                
                const double beta      = compute_beta_coeff(Mchirp, M, n0, e0);
                const double gbar0     = gbar_from_e(e0);
                const double gbar      = gbar_from_e(e);
                gamma1 = gamma0 + (gbar0-gbar)*beta;

                const double beta2     = compute_beta2_coeff(Mchirp, M, n0, e0);
                const double gbar20    = gbar2_from_e(e0,eta);
                const double gbar2     = gbar2_from_e(e,eta);
                gamma2 = (gbar20-gbar2)*beta2;

                const double beta3     = compute_beta3_coeff(Mchirp, M, n0, e0);
                const double gbar30    = gbar3_from_e(e0,eta);
                const double gbar3     = gbar3_from_e(e,eta);
                gamma3 = (gbar30-gbar3)*beta3;

                gamma  = gamma1 + gamma2 + gamma3;
            }
            
        }
    }

    return merged ? BinaryState()
                  : BinaryState { t, 
                                  bin_init.Omega, bin_init.i, 
                                  n, e, l, gamma };
}

BinaryState Evolve::solve_orbit_equations(const BinaryState &bin_init,
                                          const EvolveCoeffs_t &ev_coeffs, 
                                          const optional<PulsarTermPhase> ptphase,
                                          const double delay) const{

    const double &t0    = bin_init.t,
                 &n0    = bin_init.n,
                 &e0    = bin_init.e,
                 &gamma0= bin_init.gamma,
                 &l0    = bin_init.l,
              
                 &eta   = ev_coeffs.eta,
              
                 t      = t0 + delay;
    
    double n, e, gamma, l, gamma1, gamma2, gamma3;
    bool merged = false;

    if(e0==0){
        
        const double &A    = ev_coeffs.A,
                     &AT   = ev_coeffs.AT,
                     &AG   = ev_coeffs.AG;
                 
        const double T     = AT*delay;

        if(T>1){
            //fprintf(stderr, "Error: The binary has already merged.\n");
            //e=n=l=gamma=NAN;
            merged = true;
            throw runtime_error("tau<0 found in solve_orbit_equations");
        }
        else{
            e     = e0;
            n     = n0 / pow(1-T, 3./8);
            if(ptphase){
                l = ptphase->lp;
                gamma = ptphase->gammap;
            }
            else{
                l     = l0 + (pow(n0,-5./3) - pow(n,-5./3))/(160*A);
                gamma = gamma0 + (1./n0 - 1./n) * AG;
            }
            
        }
    }
    else{
        const double &P       = ev_coeffs.P;
        const double &tau0    = ev_coeffs.tau0;
        
        double tau = tau0 - P*delay;

        if(tau<0){
            //fprintf(stderr, "Error: The binary has already merged.\n");
            //e=n=l=gamma=NAN;
            merged = true;
            throw runtime_error("tau<0 found in solve_orbit_equations");
        }
        else{
            
            e = e_from_tau(tau);
            n = n_from_e(n0, e0, e);
            
            if(ptphase){
                l = ptphase->lp;
                gamma = ptphase->gammap;
            }
            else{
                const double &alpha  = ev_coeffs.alpha;
                const double &lbar0  = ev_coeffs.lbar0;
                const double lbar    = lbar_from_e(e);         
                l = l0 + (lbar0-lbar)*alpha;
                
                const double &beta   = ev_coeffs.beta;
                const double &gbar0  = ev_coeffs.gbar0;
                const double gbar    = gbar_from_e(e);
                gamma1 = gamma0 + (gbar0-gbar)*beta;

                const double &beta2   = ev_coeffs.beta2;
                const double &gbar20  = ev_coeffs.gbar20;
                const double gbar2    = gbar2_from_e(e,eta);
                gamma2 = (gbar20-gbar2)*beta2;

                const double &beta3   = ev_coeffs.beta3;
                const double &gbar30  = ev_coeffs.gbar30;
                const double gbar3    = gbar3_from_e(e,eta);
                gamma3 = (gbar30-gbar3)*beta3;

                gamma  = gamma1 + gamma2 + gamma3;
            }
        }
    }

    return merged ? BinaryState()
                  : BinaryState { t, 
                                  bin_init.Omega, bin_init.i, 
                                  n, e, l, gamma };
}

double Evolve::phase_err(const BinaryMass &bin_mass, 
                         const BinaryState &bin_init, 
                         const double delay) const{
    
    const BinaryState bin_end = this->solve_orbit_equations(bin_mass, bin_init, nullopt, delay);
    
    const double k0 = advance_of_periastron(bin_mass, bin_init);
    
    return bin_end.l+bin_end.gamma - (bin_init.l + bin_init.gamma + (1+k0)*bin_init.n*delay);
}
