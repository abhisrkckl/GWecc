#include "NumericalWaveform.hpp"
#include "OrbitalEvolution.hpp"
#include "mikkola.h"
#include "PN.hpp"
#include "ipow.hpp"

/*
WaveformInterpolator::WaveformInterpolator(const Signal1D &tsamples, 
					   const Signal1D &hps, 
					   const Signal1D &hxs){
		
		const size_t interp_length = tsamples.size();
		
		tmin = tsamples[0];
		tmax = tsamples[interp_length-1];
		
		acc_p = gsl_interp_accel_alloc();
		acc_x = gsl_interp_accel_alloc();
		
		spline_p = gsl_spline_alloc(gsl_interp_cspline, interp_length);
		spline_x = gsl_spline_alloc(gsl_interp_cspline, interp_length);
		
		gsl_spline_init(spline_p, &tsamples[0], &hps[0], interp_length);
		gsl_spline_init(spline_x, &tsamples[0], &hxs[0], interp_length);
}

WaveformInterpolator::~WaveformInterpolator(){
	gsl_interp_accel_free(acc_p);
	gsl_interp_accel_free(acc_x);
		
	gsl_spline_free(spline_p);
	gsl_spline_free(spline_x);
}

double WaveformInterpolator::Rp(double t){
	return gsl_spline_eval_integ(spline_p, tmin, t, acc_p);
	//return gsl_spline_eval(spline_p, t, acc_p);
}

double WaveformInterpolator::Rx(double t){
	return gsl_spline_eval_integ(spline_x, tmin, t, acc_x);
	//return gsl_spline_eval(spline_x, t, acc_x);
}


std::tuple<Signal1D,Signal1D> waveform_hpx_Num(	const BinaryMass &bin_mass, 
			  			const BinaryState &bin_init, 
			 			const Signal1D &ts){
			
	const size_t length = ts.size();
	
	Signal1D hp(length), hx(length);
		
	for(size_t i=0; i<length; i++) {
		
		const BinaryState bin_now = solve_orbit_equations(bin_mass, bin_init, ts[i]-bin_init.t);
		
		const double 	//&n = bin_now.n,
				&e = bin_now.e,
				&g = bin_now.gamma,
				&l = bin_now.l,
				
				u  = MIKKOLA(l,e),
				L  = l+g,
				
				su = sin(u),
				cu = cos(u),
				
				esu= e*su,
				ecu= e*cu;
		
		const double	//x    = PN_param_x(bin_mass, bin_now), 
				k    = advance_of_periastron(bin_mass, bin_now), 
				ephi = angular_eccentricity(bin_mass, bin_now),
		
				betaphi = (1-sqrt(1-ephi*ephi))/ephi,
				v_u	= 2*atan2(betaphi*su, 1-betaphi*cu),
				v_l	= v_u + esu,
				
				W 	= (1+k)*v_l,
				phi	= L+W,
				
				si	= sin(bin_now.i),
				ci	= cos(bin_now.i),
				
				s2phi	= sin(2*phi),
				c2phi	= cos(2*phi),
				
				OTS	= sqrt(1-e*e);
				
		const double	h0A = 1./ipow(1-ecu,2) * (	- (ci*ci+1)*(2*OTS*esu)           	*s2phi
								+ (ci*ci+1)*(2*e*e - ecu*ecu + ecu - 2)	*c2phi
								+ (si*si)  *(1-ecu)*ecu
                           				 ),
                           	h0B = 1./ipow(1-ecu,2) * (	 	    (2*OTS*esu)			*c2phi
                           					+ 	    (2*e*e - ecu*ecu + ecu - 2)*s2phi
							 )*2*ci;
				
		const double 	c2Om = cos(2*bin_now.Omega),
				s2Om = sin(2*bin_now.Omega);
				
		hp[i] = c2Om*h0A - s2Om*h0B;
		hx[i] = c2Om*h0B + s2Om*h0A;
	}
	
	return std::make_tuple(hp,hx);
	
}*/

