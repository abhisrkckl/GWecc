#include "EccentricResiduals.hpp"
#include "OrbitalEvolution.hpp"
#include "AntennaPattern.hpp"
#include "FourierWaveform.hpp"
#include "PN.hpp"
#include <iostream>

Signal1D EccentricResiduals_Anl(const BinaryMass &bin_mass,
			    	const BinaryState &bin_init,
				const SkyPosition &bin_pos,
				const SkyPosition &psr_pos,
			    	const ResidualsTerms residuals_terms,
			    	const Signal1D &ts){

	const auto [cosmu, Fp, Fx] = AntennaPattern(bin_pos, psr_pos);	
	
	if(residuals_terms==ResidualsTerms::Earth){
		
		const auto [RpE, RxE] = EccentricResiduals_px_Anl(bin_mass, bin_init, 1, ts);
		
		//const auto H0E = GWAmplitude(bin_mass, bin_init, bin_pos.DL);
		
		return -(Fp*RpE + Fx*RxE);
	}
	else if(residuals_terms==ResidualsTerms::Pulsar){
		
		const auto delay = -psr_pos.DL*(1-cosmu);
		const auto bin_psrterm = solve_orbit_equations(bin_mass, bin_init, delay);
		
		const auto [RpP, RxP] = EccentricResiduals_px_Anl(bin_mass, bin_psrterm, 1, ts + delay);
		
		//const auto H0P = GWAmplitude(bin_mass, bin_psrterm, bin_pos.DL);
		
		return (Fp*RpP + Fx*RxP);
	}
	else{
		const auto delay = -psr_pos.DL*(1-cosmu);
		const auto bin_psrterm = solve_orbit_equations(bin_mass, bin_init, delay);
		
		const auto [RpP, RxP] = EccentricResiduals_px_Anl(bin_mass, bin_psrterm, 1, ts + delay);
		const auto [RpE, RxE] = EccentricResiduals_px_Anl(bin_mass, bin_init, 1, ts);
		
		//const auto H0E = GWAmplitude(bin_mass, bin_init, bin_pos.DL),
		//	   H0P = GWAmplitude(bin_mass, bin_psrterm, bin_pos.DL);
		
		return    (Fp*RpP + Fx*RxP) 
			- (Fp*RpE + Fx*RxE);
	}
}

std::tuple<Signal1D, Signal1D> EccentricResiduals_px_Anl(const BinaryMass &bin_mass,
							 const BinaryState &bin_init,
							 const double DGW,
							 const Signal1D &ts){
	
	const auto length = ts.size();
	
	const EvolveCoeffs_t ev_coeffs = compute_evolve_coeffs(bin_mass, bin_init);
	
	const auto	cos2Omega = cos(2*bin_init.Omega),
			sin2Omega = sin(2*bin_init.Omega);
	
	Signal1D Rp(length), Rx(length);
	
	//const auto	l0	= bin_init.l 		      -       bin_init.n*bin_init.t,
	//		lambda0	= bin_init.l + bin_init.gamma - (1+k)*bin_init.n*bin_init.t;
	
	for(size_t i=0; i<length; i++){
	
		const auto bin_now = solve_orbit_equations(bin_init, ev_coeffs, ts[i]-bin_init.t);
		
		const auto k = advance_of_periastron(bin_mass, bin_now);
	
		const auto res_coeffs  = FourierResidualCoeffs(bin_mass, bin_now, 0);
		
		const auto	nt 	= bin_init.n * (ts[i]-bin_init.t),
				l	=       nt + bin_init.l,
				lambda	= (1+k)*nt + bin_init.l + bin_init.gamma;
		
		const auto 	[RA, RB]= FourierResidual_pt(res_coeffs, l, lambda);
		
		const auto 	H0	= GWAmplitude(bin_mass, bin_init, DGW);
		
		Rp[i] = H0*(cos2Omega*RA - sin2Omega*RB);
		Rx[i] = H0*(cos2Omega*RB + sin2Omega*RA);
	}
	
	return std::make_tuple(Rp, Rx);
}

