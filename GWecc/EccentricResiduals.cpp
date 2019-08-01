#include "EccentricResiduals.hpp"
#include "OrbitalEvolution.hpp"
#include "AntennaPattern.hpp"
#include "PN.hpp"

Signal1D EccentricResiduals(const BinaryMass &bin_mass,
			    const BinaryState &bin_init,
			    const SkyPosition &bin_pos,
			    const SkyPosition &psr_pos,
			    const ResidualsMethod residuals_method,
			    const ResidualsTerms residuals_terms,
			    const Signal1D &ts){

	const auto EccentricResiduals_fn =	 (residuals_method==ResidualsMethod::Anl) ? EccentricResiduals_Anl
											  : EccentricResiduals_Num;
								
	return EccentricResiduals_fn(bin_mass,
			    	     bin_init,
				     bin_pos,
				     psr_pos,
				     residuals_terms,
				     ts);
}

std::tuple<Signal1D, Signal1D> EccentricResiduals_px(const BinaryMass &bin_mass,
                                                     const BinaryState &bin_init,
						     const ResidualsMethod residuals_method,
						     const double DGW,
                                                     const Signal1D &ts){

	const auto EccentricResiduals_px_fn = (residuals_method==ResidualsMethod::Anl) ? EccentricResiduals_px_Anl
                                                                                   : EccentricResiduals_px_Num;

	return EccentricResiduals_px_fn(bin_mass, bin_init, DGW, ts);
}

double GWAmplitude(const BinaryMass &bin_mass,
		   const BinaryState &bin_init,
		   const double DGW){

	const auto eta = bin_mass.symmetric_mass_ratio(),
		   x   = PN_param_x(bin_mass, bin_init);
	
	return bin_mass.mass()*eta/DGW  *  x;
}
