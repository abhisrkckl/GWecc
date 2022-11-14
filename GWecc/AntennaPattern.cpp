#include "AntennaPattern.hpp"
#include <Eigen/Dense>
#include <cmath>
#include "ipow.hpp"

AntennaPattern antenna_pattern(const SkyPosition &bin_pos, const SkyPosition &psr_pos){
    
    const double &lambda_p = psr_pos.RA,
                 &beta_p   = psr_pos.DEC,
                 &lambda   = bin_pos.RA,
                 &beta     = bin_pos.DEC;
    
    // Pulsar vector
    const double n1 = cos(lambda_p)*cos(beta_p),
                 n2 = sin(lambda_p)*cos(beta_p),
                 n3 = sin(beta_p);
    
    const double cosmu = cos(beta)*cos(beta_p)*cos(lambda-lambda_p) + sin(beta)*sin(beta_p);
    
    
    double Fp = 0, Fx = 0;
    if(cosmu!=1){
        const double    e11p = ipow(sin(lambda),2)-ipow(cos(lambda),2)*ipow(sin(beta),2),
                        e21p = -sin(lambda)*cos(lambda)*(ipow(sin(beta),2)+1),
                        e31p = cos(lambda)*sin(beta)*cos(beta),
                        e12p = -sin(lambda)*cos(lambda)*(ipow(sin(beta),2)+1),
                        e22p = ipow(cos(lambda),2)-ipow(sin(lambda),2)*ipow(sin(beta),2),
                        e32p = sin(lambda)*sin(beta)*cos(beta),
                        e13p = cos(lambda)*sin(beta)*cos(beta),
                        e23p = sin(lambda)*sin(beta)*cos(beta),
                        e33p = -ipow(cos(beta),2);
    
        Fp = (n1*(n1*e11p+n2*e12p+n3*e13p)+
              n2*(n1*e21p+n2*e22p+n3*e23p)+
              n3*(n1*e31p+n2*e32p+n3*e33p));
        Fp *= 0.5 / (1-cosmu);
        
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
        Fx *= 0.5 / (1-cosmu);
        
    }
    
    return {cosmu, Fp, Fx};
}
