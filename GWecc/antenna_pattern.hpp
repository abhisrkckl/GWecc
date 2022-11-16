#ifndef _AntennaPattern_hpp_
#define _AntennaPattern_hpp_ 1

#include <tuple>
#include "binary.hpp"

struct AntennaPattern{
    double cosmu, Fp, Fx;
};

AntennaPattern antenna_pattern(const SkyPosition &bin_pos, const SkyPosition &psr_pos);

double pulsar_term_delay(const double cosmu, const double Dp, const double z);

#endif
