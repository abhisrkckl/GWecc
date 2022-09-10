#ifndef _FeStat_hpp_
#define _FeStat_hpp_

#include "NumericalWaveform.hpp"
#include <array>

std::array<Signal1D, 6> FeStatFuncs(const BinaryMass &bin_mass,
                                    const BinaryState &bin_init,
                                    const SkyPosition &bin_pos,
                                    const SkyPosition &psr_pos,
                                    const Signal1D &ts);

std::array<Signal1D, 6> FeStatFuncs_Adb(const BinaryMass &bin_mass,
                                        const BinaryState &bin_init,
                                        const SkyPosition &bin_pos,
                                        const SkyPosition &psr_pos,
                                        const Signal1D &ts);

#endif
