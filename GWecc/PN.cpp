
#include <cmath>
#include "PN.hpp"
#include "ipow.hpp"

double PN_param_xi(const BinaryMass &bin_mass,
                   const BinaryState &bin_state){
    return bin_mass.mass() * bin_state.n;
}

double PN_param_x(const BinaryMass &bin_mass,
                  const BinaryState &bin_state){

    const double xi = bin_mass.mass() * bin_state.n;
    const double k  = advance_of_periastron(bin_mass, bin_state);
    
    return cbrt(ipow((1+k)*xi, 2));
}

double advance_of_periastron(const BinaryMass &binmass, const BinaryState &binstate){
    const double xi  = cbrt( ipow(binmass.mass()*binstate.n, 2) ),
                 e   = binstate.e,
                 OTS = sqrt(1-e*e),
                 eta = binmass.symmetric_mass_ratio();
    
    return       3*xi/ipow(OTS,2) 
                    + ((78 + ipow(e,2)*(51 - 26*eta) - 28*eta)*ipow(xi,2))/(4.*ipow(OTS,4)) 
                    + ((18240 - 25376*eta + 896*ipow(eta,2) + ipow(e,4)*(2496 - 1760*eta + 1040*ipow(eta,2)) 
                        + (1920 + ipow(e,2)*(3840 - 1536*eta) - 768*eta)*OTS + 492*eta*ipow(M_PI,2) 
                        + ipow(e,2)*(28128 - 27840*eta + 5120*ipow(eta,2) 
                        + 123*eta*ipow(M_PI,2)))*ipow(xi,3))/(128.*ipow(OTS,6));
}

double angular_eccentricity(const BinaryMass &bin_mass, 
                            const BinaryState &bin_state){
    
    const double x   = PN_param_x(bin_mass, bin_state),
                 e     = bin_state.e,
                 OTS = sqrt(1 - e*e),
                 eta = bin_mass.symmetric_mass_ratio();
    //constexpr double PI2 = M_PI*M_PI;
    
    return  e * (1     + x   *(4 - eta) 
                + x*x *(4*(-12*(26 + 15*OTS) + eta*(17 + 72*OTS + eta)) + e*e *(1152 + eta*(-659 + 41*eta)))/(96*(-1 + e*e)));
            /*
        + (x**3 *(-70*e**4*(-12288 + eta*(11233 + 5*eta*(-383 + 3*eta))) 
        + 20*(1344*(54 + 65*ots) + eta*(861*(1 + ots)*pi2 + 4*(-33431 - 29960*ots - 3458*eta + 3360*ots*eta)))
        + 3*e**2*(-8960*(76 + 35*ots) + eta*(-1435*pi2 + 4*(8*(12983 + 5040*ots) + 35*eta*(-1319 + 32*ots + 45*eta))))))/(26880*(-1 + e**2)**2)))    */
    
    //return  e*(1 + x*(4-eta));
}
