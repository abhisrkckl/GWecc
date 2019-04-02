
#include "AntennaPattern.hpp"
#include <Eigen/Dense>
#include <cmath>

std::tuple<double,double,double> AntennaPattern(const SkyPosition &bin_pos, const SkyPosition &psr_pos){

	const auto 	sin_th_gw = cos(bin_pos.DEC),
			cos_th_gw = sin(bin_pos.DEC),
			sin_ph_gw = sin(bin_pos.RA),
			cos_ph_gw = cos(bin_pos.RA),
			
			sin_th_ps = cos(psr_pos.DEC),
			cos_th_ps = sin(psr_pos.DEC),
			sin_ph_ps = sin(psr_pos.RA),
			cos_ph_ps = cos(psr_pos.RA);
	
	// Equations 4-6 of https://arxiv.org/pdf/1204.4218.pdf
	const Eigen::Vector3d 	Rhat( sin_th_gw*cos_ph_gw,  sin_th_gw*sin_ph_gw, cos_th_gw),	// Line of sight to GW source
				mhat(-sin_ph_gw,            cos_ph_gw,           0        ),	// Direction of increasing RA of GW source
				nhat(-cos_th_gw*cos_ph_gw, -cos_th_gw*sin_ph_gw, sin_th_gw),	// Direction of increasing DEC of GW source
	
				phat( sin_th_ps*cos_ph_ps,  sin_th_ps*sin_ph_ps, cos_th_ps);	// Line of sight to pulsar	
	
	// Dot products
	const auto 	mhat_phat = mhat.dot(phat),
			nhat_phat = nhat.dot(phat),
			Rhat_phat = Rhat.dot(phat);
	
	// Equations 9-10 of https://arxiv.org/pdf/1204.4218.pdf
	double Fp, Fx;
	if(Rhat_phat != 1){
		Fp   = 0.5 * (mhat_phat*mhat_phat - nhat_phat*nhat_phat) / (1-Rhat_phat);
		Fx   = mhat_phat * nhat_phat / (1-Rhat_phat);
	}
	else{
		Fp=Fx = NAN;
	}

	return std::make_tuple(Rhat_phat, Fp, Fx);
}
