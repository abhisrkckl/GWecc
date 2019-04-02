#ifndef _AntennaPattern_hpp_
#define _AntennaPattern_hpp_ 1

#include <tuple>
#include "Binary.hpp"

/*
 * Returns:
 *	[cos(eta), Fp, Fx]
 */
std::tuple<double,double,double> AntennaPattern(const SkyPosition &bin_pos, const SkyPosition &psr_pos);

#endif
