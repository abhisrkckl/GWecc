#include "OrbitalEvolution.hpp"
#include "PN.hpp"
#include "ipow.hpp"
#include "mikkola.h"
#include <array>
#include <cmath>

std::array<double,2> coords(const BinaryMass &bin_mass, const BinaryState &bin_now){
    
    const auto &n       = bin_now.n,
               &e       = bin_now.e,
               &g       = bin_now.gamma,
               &l       = bin_now.l,
               
               M        = bin_mass.mass(),
               eta      = bin_mass.symmetric_mass_ratio(), 
                    
               u        = mikkola(l,e),
               L        = l+g,
               
               su       = sin(u),
               cu       = cos(u),
                       
               esu      = e*su,
               ecu      = e*cu,
               
               x        = PN_param_x(bin_mass, bin_now),
               sqrtx    = sqrt(x),
               k        = advance_of_periastron(bin_mass, bin_now),
               ephi     = angular_eccentricity(bin_mass, bin_now),
               
               betaphi  = (ephi>1e-15)  ? (1-sqrt(1-ephi*ephi))/ephi 
                                        : (e/2. + ipow(e,3)/8. + ipow(e,5)/16.),
               v_u      = 2*atan2(betaphi*su, 1-betaphi*cu),
               v_l      = v_u + esu,
                    
               W        = (1+k)*v_l,
               phi      = L+W,
               
               a        = cbrt(M/(n*n)),
               ar       = a*(1-x/3*(9-eta)),
               er       = e*(1+x/2*(8-3*eta)),
               r        = ar*(1-er*cu);
    
    return {r,phi};
}
