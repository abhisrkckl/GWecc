#ifndef _PN_hpp_
#define _PN_hpp_ 1

#include "binary.hpp"

double pn_param_xi(const BinaryMass &bin_mass,
                   const BinaryState &bin_state);

double advance_of_periastron(const BinaryMass &bin_mass,
                             const BinaryState &bin_state);

double pn_param_x(const BinaryMass &bin_mass,
                  const BinaryState &bin_state);

double angular_eccentricity(const BinaryMass &bin_mass, 
                            const BinaryState &bin_state);

#endif
