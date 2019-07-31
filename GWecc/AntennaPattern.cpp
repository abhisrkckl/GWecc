
#include "AntennaPattern.hpp"
#include <Eigen/Dense>
#include <cmath>

/*
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
}*/

std::tuple<double,double,double> AntennaPattern(const SkyPosition &bin_pos, const SkyPosition &psr_pos){
    
    const double &lambda_p = psr_pos.RA,
                 &beta_p   = psr_pos.DEC,
                 &lambda   = bin_pos.RA,
                 &beta     = bin_pos.DEC;
    
    // Pulsar vector
    const double n1 = cos(lambda_p)*cos(beta_p),
                 n2 = sin(lambda_p)*cos(beta_p),
                 n3 = sin(beta_p);
    
    const double cos_theta = cos(beta)*cos(beta_p)*cos(lambda-lambda_p) + sin(beta)*sin(beta_p);
    
    
    double Fp = 0, Fx = 0;
    if(cos_theta!=1){
        const double    e11p = pow(sin(lambda),2)-pow(cos(lambda),2)*pow(sin(beta),2),
                        e21p = -sin(lambda)*cos(lambda)*(pow(sin(beta),2)+1),
                        e31p = cos(lambda)*sin(beta)*cos(beta),
                        e12p = -sin(lambda)*cos(lambda)*(pow(sin(beta),2)+1),
                        e22p = pow(cos(lambda),2)-pow(sin(lambda),2)*pow(sin(beta),2),
                        e32p = sin(lambda)*sin(beta)*cos(beta),
                        e13p = cos(lambda)*sin(beta)*cos(beta),
                        e23p = sin(lambda)*sin(beta)*cos(beta),
                        e33p = -pow(cos(beta),2);
    
        Fp = (n1*(n1*e11p+n2*e12p+n3*e13p)+
                n2*(n1*e21p+n2*e22p+n3*e23p)+
                n3*(n1*e31p+n2*e32p+n3*e33p));
        Fp *= 0.5 / (1-cos_theta);
        
        /************/
        
        const double    e11c = sin(2*lambda)*sin(beta),
                        e21c = -cos(2*lambda)*sin(beta),
                        e31c = -sin(lambda)*cos(beta),
                        e12c = -cos(2*lambda)*sin(beta),
                        e22c = -sin(2*lambda)*sin(beta),
                        e32c = cos(lambda)*cos(beta),
                        e13c = -sin(lambda)*cos(beta),
                        e23c = cos(lambda)*cos(beta),
                        e33c  = 0;

        Fx = (n1*(n1*e11c+n2*e12c+n3*e13c)+
              n2*(n1*e21c+n2*e22c+n3*e23c)+
              n3*(n1*e31c+n2*e32c+n3*e33c));
        Fx *= 0.5 / (1-cos_theta);
        
    }
    
    return std::make_tuple(cos_theta, Fp, Fx);
}
