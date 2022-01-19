#include "GWecc.hpp"
#include "EccentricResiduals.hpp"
#include "FeStat.hpp"
#include "AntennaPattern.hpp"

static constexpr double MSun_to_s   = 4.92703806e-6,       // Solar mass in s (geometric units)
                        parsec_to_s = 102927125.0,         // Parsec in s (geometric units)
                        year_to_s   = 365.25*24*3600,
                        day_to_s    = 24*3600;

//BinaryState solve_orbit_equations(const BinaryMass &bin_mass, const BinaryState &bin_init, const double delay);

/*
 * M        is Total Mass               in SolarMass
 * q        is Mass Ratio
 *
 * Omega    is Polarization Angle       in rad
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
std::vector<double> EccentricResiduals( const double M, const double q,
                                        const double Omega, const double i, 
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
                                 Omega, i,
                                 n0, e0, l0, gamma0 };
                     
    const SkyPosition bin_pos {parsec_to_s*D_GW, RA_GW, DEC_GW, z},
                      psr_pos {parsec_to_s*D_P,  RA_P,  DEC_P,  0};
        
    Signal1D tzs(_ts.data(), _ts.size());
    tzs *= day_to_s/(1+z);
    
    Signal1D result = EccentricResiduals(bin_mass, bin_init, bin_pos, psr_pos, 
                                         residuals_method, 
                                         residuals_terms,
                                         tzs);
    result *= (1+z);
    
    return std::vector<double>(std::begin(result), std::end(result));
}


/*
 * M        is Total Mass               in SolarMass
 * q        is Mass Ratio
 *
 * Omega    is Polarization Angle       in rad
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
std::vector<std::vector<double> > EccentricResiduals_px(const double M, const double q,
                                                        const double Omega, const double i,
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
                                 Omega, i,
                                 n0, e0, l0, gamma0 };

    Signal1D tzs(_ts.data(), _ts.size());
    tzs *= day_to_s/(1+z);

    Signal1D _Rp, _Rx;
    std::tie(_Rp, _Rx) = EccentricResiduals_px(bin_mass, bin_init,
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


std::vector<std::vector<double> > EccentricWaveform_px( const double M, const double q,
                                                        const double Omega, const double i,
                                                        const double t0, const double Pb0E, const double e0, const double l0, const double gamma0,
                                                        const double DGW, const double z,
                                                        const std::vector<double> _ts){


    const double Pb0 = Pb0E * year_to_s / (1+z),
                 n0 = 2*M_PI/Pb0;

    const BinaryMass bin_mass(M*MSun_to_s, q);

    const BinaryState bin_init { day_to_s*t0,
                                 Omega, i,
                                 n0, e0, l0, gamma0 };

    Signal1D tzs(_ts.data(), _ts.size());
    tzs *= day_to_s/(1+z);

    Signal1D _hp, _hx;
    std::tie(_hp, _hx) = EccentricWaveform_px( bin_mass, bin_init, 
                                                     DGW*parsec_to_s,
                                                     tzs);

    std::vector<double>     hp(std::begin(_hp), std::end(_hp)),
                 			hx(std::begin(_hx), std::end(_hx));

    return std::vector<std::vector<double> > {hp, hx};
}

std::vector<double> AntennaPattern(const double RA_GW, const double DEC_GW, const double RA_P, const double DEC_P){
    double cos_th, Fp, Fx;
    std::tie(cos_th, Fp, Fx) = AntennaPattern(SkyPosition{0, RA_GW, DEC_GW, 0}, SkyPosition{0, RA_P, DEC_P, 0});
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

std::vector<std::vector<double> > FeStatFuncs(const double M, const double q,
                                              const double t0, const double Pb0E, const double e0, 
                                              const double D_GW, const double RA_GW, const double DEC_GW, 
                                              const double RA_P,  const double DEC_P, 
                                              const double z,
                                              const std::vector<double> _ts){


    const double Pb0 = Pb0E * year_to_s / (1+z),
                 n0 = 2*M_PI/Pb0;
    
    const SkyPosition bin_pos {D_GW, RA_GW, DEC_GW, z},
                      psr_pos {0,  RA_P,  DEC_P, 0};
    
    const BinaryMass bin_mass(M*MSun_to_s, q);

    const BinaryState bin_init { day_to_s*t0,
                                 0, 0,
                                 n0, e0, 0, 0 };

    Signal1D tzs(_ts.data(), _ts.size());
    tzs *= day_to_s/(1+z);

    std::array<Signal1D,10> A = FeStatFuncs(bin_mass, bin_init, bin_pos, psr_pos, tzs);
        
    std::vector<double>     A0(std::begin(A[0]), std::end(A[0])),
                 			A1(std::begin(A[1]), std::end(A[1])),
                 			A2(std::begin(A[2]), std::end(A[2])),
                 			A3(std::begin(A[3]), std::end(A[3])),
                 			A4(std::begin(A[4]), std::end(A[4])),
                 			A5(std::begin(A[5]), std::end(A[5])),
                 			A6(std::begin(A[6]), std::end(A[6])),
                 			A7(std::begin(A[7]), std::end(A[7])),
                 			A8(std::begin(A[8]), std::end(A[8])),
                            A9(std::begin(A[9]), std::end(A[9]));

    return std::vector<std::vector<double> > {A0, A1, A2, A3, A4, A5, A6, A7, A8, A9};
}

/*
std::vector<std::vector<double> > FeStatFuncs_h(const double M, const double q,
                                              const double t0, const double Pb0E, const double e0, const double l0, const double gamma0,
                                              const double RA_GW, const double DEC_GW, 
                                              const double RA_P,  const double DEC_P, 
                                              const double z,
                                              const std::vector<double> _ts){


    const double Pb0 = Pb0E * year_to_s / (1+z),
                 n0 = 2*M_PI/Pb0;
    
    const SkyPosition bin_pos {0, RA_GW, DEC_GW, 0},
                      psr_pos {0,  RA_P,  DEC_P, 0};
    
    const BinaryMass bin_mass(M*MSun_to_s, q);

    const BinaryState bin_init { day_to_s*t0,
                                 0, 0,
                                 n0, e0, l0, gamma0 };

    Signal1D tzs(_ts.data(), _ts.size());
    tzs *= day_to_s/(1+z);

    std::array<Signal1D,5> B = FeStatFuncs_h(bin_mass, bin_init, bin_pos, psr_pos, tzs);
        
    std::vector<double>     B0(std::begin(B[0]), std::end(B[0])),
                 			B1(std::begin(B[1]), std::end(B[1])),
                 			B2(std::begin(B[2]), std::end(B[2])),
                 			B3(std::begin(B[3]), std::end(B[3])),
                 			B4(std::begin(B[4]), std::end(B[4]));

    return std::vector<std::vector<double> > {B0, B1, B2, B3, B4};
}
*/