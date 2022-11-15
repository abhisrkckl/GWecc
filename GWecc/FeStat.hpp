#ifndef _FeStat_hpp_
#define _FeStat_hpp_

#include <array>
#include "numerical_residuals.hpp"

std::array<Signal1D, 6> fe_stat_funcs(const BinaryMass &bin_mass,
                                    const BinaryState &bin_init,
                                    const SkyPosition &bin_pos,
                                    const SkyPosition &psr_pos,
                                    const Signal1D &ts);

std::array<Signal1D, 6> fe_stat_funcs_Adb(const BinaryMass &bin_mass,
                                        const BinaryState &bin_init,
                                        const SkyPosition &bin_pos,
                                        const SkyPosition &psr_pos,
                                        const Signal1D &ts);

#endif
