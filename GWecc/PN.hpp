#ifndef _PN_hpp_
#define _PN_hpp_ 1

#include "Binary.hpp"

double PN_param_xi(const BinaryMass &bin_mass,
                   const BinaryState &bin_state);

double advance_of_periastron(const BinaryMass &bin_mass,
                             const BinaryState &bin_state);

double PN_param_x(const BinaryMass &bin_mass,
                  const BinaryState &bin_state);

double angular_eccentricity(const BinaryMass &bin_mass, 
                            const BinaryState &bin_state);

#endif
