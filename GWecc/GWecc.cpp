#include "GWecc.hpp"
#include "eccentric_residuals.hpp"
#include "FeStat.hpp"
#include "antenna_pattern.hpp"

static constexpr double MSun_to_s   = 4.92703806e-6,       // Solar mass in s (geometric units)
                        parsec_to_s = 102927125.0,         // Parsec in s (geometric units)
                        year_to_s   = 365.25*24*3600,
                        day_to_s    = 24*3600;

//BinaryState solve_orbit_equations(const BinaryMass &bin_mass, const BinaryState &bin_init, const double delay);

/*
 * M        is Total Mass               in SolarMass
 * q        is Mass Ratio
 *
 * psi      is Polarization Angle       in rad
 * i        is Inclination              in rad
 *
 * t0       is Epoch                    in MJD
 * Pb0E     is Orbital period at t0     in year     in Earth frame
 * e0       is Eccentricity at t0 
 * l0       is Mean Anomaly at t0       in rad
 * gamma0   is Periastron Angle at t0   in rad
 * 
 * D_GW     is Distance to Binary       in parsec
 * RA_GW    is RA of Binary             in rad
 * DEC_GW   is DEC of Binary            in rad
 *
 * D_P      is Distance to Pulsar       in parsec
 * RA_P     is RA of Pulsar             in rad
 * DEC_P    is DEC of Pulsar            in rad
 *
 * z        is redshift of GW source
 *
 * residuals_method is one of {Anl,Num}
 *    
 * ts       are TOAs                    in MJD        in Earth frame
 */
std::vector<double> eccentric_residuals(const double M, const double q,
                                        const double psi, const double i, 
                                        const double t0, const double Pb0E, const double e0, const double l0, const double gamma0,
                                        const double D_GW, const double RA_GW, const double DEC_GW, 
                                        const double D_P,  const double RA_P,  const double DEC_P, 
                                        const double z,
                                        const ResidualsMethod residuals_method,
                                        const ResidualsTerms residuals_terms,
                                        const std::vector<double> _ts){
    
    const BinaryMass bin_mass(MSun_to_s*M, q);
    
    const double Pb0 = Pb0E * year_to_s / (1+z),
                 n0 = 2*M_PI/Pb0; 
    
    const BinaryState bin_init { day_to_s*t0,
                                 psi, i,
                                 n0, e0, l0, gamma0 };
                     
    const SkyPosition bin_pos {parsec_to_s*D_GW, RA_GW, DEC_GW, z},
                      psr_pos {parsec_to_s*D_P,  RA_P,  DEC_P,  0};
        
    Signal1D tzs(_ts.data(), _ts.size());
    tzs *= day_to_s/(1+z);
    
    Signal1D result = eccentric_residuals(bin_mass, bin_init, bin_pos, psr_pos, 
                                         residuals_method, 
                                         residuals_terms,
                                         tzs);
    result *= (1+z);
    
    return std::vector<double>(std::begin(result), std::end(result));
}

std::vector<std::vector<double>> eccentric_residuals_and_waveform( 
                                        const double M, const double q,
                                        const double psi, const double i, 
                                        const double t0, const double Pb0E, const double e0, const double l0, const double gamma0,
                                        const double D_GW, const double RA_GW, const double DEC_GW, 
                                        const double D_P,  const double RA_P,  const double DEC_P, 
                                        const double z,
                                        const ResidualsTerms residuals_terms,
                                        const std::vector<double> _ts){
    
    const BinaryMass bin_mass(MSun_to_s*M, q);
    
    const double Pb0 = Pb0E * year_to_s / (1+z),
                 n0 = 2*M_PI/Pb0; 
    
    const BinaryState bin_init { day_to_s*t0,
                                 psi, i,
                                 n0, e0, l0, gamma0 };
                     
    const SkyPosition bin_pos {parsec_to_s*D_GW, RA_GW, DEC_GW, z},
                      psr_pos {parsec_to_s*D_P,  RA_P,  DEC_P,  0};
        
    Signal1D tzs(_ts.data(), _ts.size());
    tzs *= day_to_s/(1+z);
    
    auto [R, h] = eccentric_residuals_and_waveform(bin_mass, bin_init, bin_pos, psr_pos, 
                                                residuals_terms,
                                                tzs);
    R *= (1+z);
    
    const std::vector<double> result_R(std::begin(R), std::end(R));
    const std::vector<double> result_h(std::begin(h), std::end(h));

    return {result_R, result_h};
}

std::vector<double> eccentric_waveform( const double M, const double q,
                                        const double psi, const double i, 
                                        const double t0, const double Pb0E, const double e0, const double l0, const double gamma0,
                                        const double D_GW, const double RA_GW, const double DEC_GW, 
                                        const double D_P,  const double RA_P,  const double DEC_P, 
                                        const double z,
                                        const ResidualsTerms residuals_terms,
                                        const std::vector<double> _ts){
    
    const BinaryMass bin_mass(MSun_to_s*M, q);
    
    const double Pb0 = Pb0E * year_to_s / (1+z),
                 n0 = 2*M_PI/Pb0; 
    
    const BinaryState bin_init { day_to_s*t0,
                                 psi, i,
                                 n0, e0, l0, gamma0 };
                     
    const SkyPosition bin_pos {parsec_to_s*D_GW, RA_GW, DEC_GW, z},
                      psr_pos {parsec_to_s*D_P,  RA_P,  DEC_P,  0};
        
    Signal1D tzs(_ts.data(), _ts.size());
    tzs *= day_to_s/(1+z);
    
    Signal1D result = eccentric_waveform(bin_mass, bin_init, bin_pos, psr_pos, 
                                        residuals_terms,
                                        tzs);
    
    return std::vector<double>(std::begin(result), std::end(result));
}


/*
 * M        is Total Mass               in SolarMass
 * q        is Mass Ratio
 *
 * psi      is Polarization Angle       in rad
 * i        is Inclination              in rad
 *
 * t0       is Epoch                    in MJD
 * Pb0E     is Orbital period at t0     in year     in Earth frame
 * e0       is Eccentricity at t0 
 * l0       is Mean Anomaly at t0       in rad
 * gamma0   is Periastron Angle at t0   in rad
 * 
 * D_GW     is Distance to Binary       in parsec
 * delay    is Delay of pulsar term     in parsec
 * z        is redshift of GW source
 *
 * residuals_method is one of {Anl,Num}
 *    
 * ts       are TOAs                    in MJD        in Earth frame
 */
std::vector<std::vector<double> > eccentric_residuals_px(const double M, const double q,
                                                        const double psi, const double i,
                                                        const double t0, const double Pb0E, const double e0, const double l0, const double gamma0,
                                                        const double DGW, const double delay,
                                                        const double z,
                                                        const ResidualsMethod residuals_method,
                                                        const ResidualsTerms residuals_terms,
                                                        const std::vector<double> _ts){

    const BinaryMass bin_mass(MSun_to_s*M, q);

    const double Pb0 = Pb0E * year_to_s / (1+z),
                 n0 = 2*M_PI/Pb0; 

    const BinaryState bin_init { day_to_s*t0,
                                 psi, i,
                                 n0, e0, l0, gamma0 };

    Signal1D tzs(_ts.data(), _ts.size());
    tzs *= day_to_s/(1+z);

    auto [_Rp, _Rx] = eccentric_residuals_px(bin_mass, bin_init,
                                            parsec_to_s*DGW, parsec_to_s*delay,
                                            residuals_method,
                                            residuals_terms,
                                            tzs);
    _Rp *= (1+z);
    _Rx *= (1+z);

    std::vector<double>  Rp(std::begin(_Rp), std::end(_Rp)),
                         Rx(std::begin(_Rx), std::end(_Rx));

    return std::vector<std::vector<double> > {Rp, Rx};
}


std::vector<std::vector<double> > eccentric_waveform_px( const double M, const double q,
                                                        const double psi, const double i,
                                                        const double t0, const double Pb0E, const double e0, const double l0, const double gamma0,
                                                        const double DGW, const double z,
                                                        const std::vector<double> _ts){

    const double Pb0 = Pb0E * year_to_s / (1+z),
                 n0 = 2*M_PI/Pb0;

    const BinaryMass bin_mass(M*MSun_to_s, q);

    const BinaryState bin_init { day_to_s*t0,
                                 psi, i,
                                 n0, e0, l0, gamma0 };

    Signal1D tzs(_ts.data(), _ts.size());
    tzs *= day_to_s/(1+z);

    const auto [_hp, _hx] = eccentric_waveform_px(bin_mass, bin_init, 
                                                 DGW*parsec_to_s,
                                                 tzs);

    std::vector<double>     hp(std::begin(_hp), std::end(_hp)),
                 			hx(std::begin(_hx), std::end(_hx));

    return std::vector<std::vector<double> > {hp, hx};
}

std::vector<double> antenna_pattern(const double RA_GW, const double DEC_GW, const double RA_P, const double DEC_P){
    const auto [cos_th, Fp, Fx] = antenna_pattern(SkyPosition{0, RA_GW, DEC_GW, 0}, SkyPosition{0, RA_P, DEC_P, 0});
    return std::vector<double>{cos_th, Fp, Fx};
}

bool mergeq(const double M, const double q, 
            const double Pb0E, const double e0,
            const double z,
            const double t0, const double max_toa){
    
    const double n0 = 2*M_PI/(Pb0E*year_to_s) * (1+z) ;

    const BinaryState bin_init {t0*day_to_s,
                                0, 0,
                                n0, e0, 0, 0};
    
    const BinaryMass bin_mass {M*MSun_to_s, q};
    
    const double delay = (max_toa-t0)*day_to_s;

    BinaryState bin_last = solve_orbit_equations(bin_mass, bin_init, delay);

    return bin_last.merged;
}

std::vector<std::vector<double> > fe_stat_funcs(const double M, const double q,
                                              const double t0, const double Pb0E, const double e0, const double l0,
                                              const double D_GW, const double RA_GW, const double DEC_GW, 
                                              const double RA_P,  const double DEC_P, 
                                              const double z,
                                              const std::vector<double> _ts,
                                              const ResidualsMethod residuals_method){


    const double Pb0 = Pb0E * year_to_s / (1+z),
                 n0 = 2*M_PI/Pb0;
    
    const SkyPosition bin_pos {D_GW, RA_GW, DEC_GW, z},
                      psr_pos {0,  RA_P,  DEC_P, 0};
    
    const BinaryMass bin_mass(M*MSun_to_s, q);

    const BinaryState bin_init { day_to_s*t0,
                                 0, 0,
                                 n0, e0, l0, 0 };

    Signal1D tzs(_ts.data(), _ts.size());
    tzs *= day_to_s/(1+z);

    std::array<Signal1D,6> A;
    if(residuals_method == ResidualsMethod::Num){
        A = fe_stat_funcs(bin_mass, bin_init, bin_pos, psr_pos, tzs);
    }
    else if(residuals_method == ResidualsMethod::Adb){
        A = fe_stat_funcs_Adb(bin_mass, bin_init, bin_pos, psr_pos, tzs);
    }
    else{
        throw std::invalid_argument("Only Num and Adb methods are supported in fe_stat_funcs.");
    }

    for(auto Ai : A){
        Ai *= (1+z);
    }
    
    std::vector<double> A0(std::begin(A[0]), std::end(A[0])),
                 		A1(std::begin(A[1]), std::end(A[1])),
                 		A2(std::begin(A[2]), std::end(A[2])),
                 		A3(std::begin(A[3]), std::end(A[3])),
                 		A4(std::begin(A[4]), std::end(A[4])),
                 		A5(std::begin(A[5]), std::end(A[5]));

    return std::vector<std::vector<double>> {A0, A1, A2, A3, A4, A5};
}
