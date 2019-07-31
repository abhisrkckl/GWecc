#include "EccentricResiduals_py.hpp"
#include "EccentricResiduals.hpp"
#include "AntennaPattern.hpp"

static constexpr double MSun_geom	= 4.92703806e-6;	// Solar mass in s (geometric units)
static constexpr double parsec_geom 	= 102927125.0;		// Parsec in s (geometric units)


/*
 * M		is Total Mass			in SolarMass
 * q		is Mass Ratio
 *
 * Omega	is Polarization Angle 		in rad
 * i		is Inclination 			in rad
 *
 * t0		is Epoch			in s
 * n0		is Mean Motion at t0		in rad/s
 * e0		is Eccentricity at t0
 * l0		is Mean Anomaly at t0		in rad
 * gamma0	is Periastron Angle at t0	in rad
 * 
 * D_GW		is Distance to Binary		in parsec
 * RA_GW	is RA of Binary			in rad
 * DEC_GW	is DEC of Binary		in rad
 *
 * D_P		is Distance to Pulsar		in parsec
 * RA_P		is RA of Pulsar			in rad
 * DEC_P	is DEC of Pulsar		in rad
 *
 * residuals_method	is one of	{Anl,FFT,Num}
 *	
 * ts 		is the time instants		in s
 */
std::vector<double> EccentricResiduals( const double M, const double q,
					const double Omega, const double i, 
					const double t0, const double n0, const double e0, const double l0, const double gamma0,
					const double D_GW, const double RA_GW, const double DEC_GW, 
					const double D_P,  const double RA_P,  const double DEC_P, 
					const ResidualsMethod residuals_method,
					const ResidualsTerms residuals_terms,
					const std::vector<double> _ts){
	
	const BinaryMass bin_mass(MSun_geom*M, q);
	
	const BinaryState bin_init { t0,
				     Omega, i,
				     n0, e0, l0, gamma0 };
				     
	const SkyPosition bin_pos {parsec_geom*D_GW, RA_GW, DEC_GW},
			  psr_pos {parsec_geom*D_P,  RA_P,  DEC_P};
		
	const Signal1D ts(_ts.data(), _ts.size());
	
	/*
	    This should be decided based on residuals_method once EccentricResiduals_AB_FFT and EccentricResiduals_AB_Num are implemented.
	    Also, passing ResidualsTerms should be implemented.
	*/
	const Signal1D result = EccentricResiduals(bin_mass, bin_init, bin_pos, psr_pos, 
					      	   residuals_method, 
					      	   residuals_terms,
					      	   ts);
	
	return std::vector<double>(std::begin(result), std::end(result));
}

std::vector<std::vector<double> > EccentricResiduals_px(const double M, const double q,
                                        	      	const double Omega, const double i,
                                        	      	const double t0, const double n0, const double e0, const double l0, const double gamma0,
                                        	      	const double DGW,
	                                             	const ResidualsMethod residuals_method,
                                        	     	const std::vector<double> _ts){

        const BinaryMass bin_mass(MSun_geom*M, q);

        const BinaryState bin_init { t0,
                                     Omega, i,
                                     n0, e0, l0, gamma0 };

        const Signal1D ts(_ts.data(), _ts.size());

        /*
            This should be decided based on residuals_method once EccentricResiduals_AB_FFT and EccentricResiduals_AB_Num are implemented.
            Also, passing ResidualsTerms should be implemented.
        */
	Signal1D _Rp, _Rx;
	std::tie(_Rp, _Rx) = EccentricResiduals_px(bin_mass, bin_init, 
                                              	   residuals_method,
                                              	   DGW,
                                              	   ts);

	std::vector<double> 	Rp(std::begin(_Rp), std::end(_Rp)),
			 	Rx(std::begin(_Rx), std::end(_Rx));

        return std::vector<std::vector<double> > {Rp, Rx};
}

std::vector<std::vector<double> > EccentricWaveform_px(	const double M, const double q,
                                        	      	const double Omega, const double i,
                                        	      	const double t0, const double n0, const double e0, const double l0, const double gamma0,
                                        	      	const double DGW,
                                        	     	const std::vector<double> _ts){

        const BinaryMass bin_mass(MSun_geom*M, q);

        const BinaryState bin_init { t0,
                                     Omega, i,
                                     n0, e0, l0, gamma0 };

        const Signal1D ts(_ts.data(), _ts.size());

        /*
            This should be decided based on residuals_method once EccentricResiduals_AB_FFT and EccentricResiduals_AB_Num are implemented.
            Also, passing ResidualsTerms should be implemented.
        */
	Signal1D _hp, _hx;
	std::tie(_hp, _hx) = EccentricWaveform_px( bin_mass, bin_init, 
                                              	   DGW,
                                              	   ts);

	std::vector<double> 	hp(std::begin(_hp), std::end(_hp)),
			 	hx(std::begin(_hx), std::end(_hx));

        return std::vector<std::vector<double> > {hp, hx};
}

std::vector<double> AntennaPattern(const double RA_GW, const double DEC_GW, const double RA_P, const double DEC_P){
    const auto [cos_th, Fp, Fx] = AntennaPattern(SkyPosition{0, RA_GW, DEC_GW}, SkyPosition{0, RA_P, DEC_P});
    return std::vector<double>{cos_th, Fp, Fx};
}
