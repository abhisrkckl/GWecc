#ifndef _EccentricResiduals_py_hpp_
#define _EccentricResiduals_py_hpp_

#include <vector>
#include "ResidualsOptions.hpp"
#include "Binary.hpp"

/*
 * Computes R(t)
 */
std::vector<double> eccentric_residuals( const double M, const double q,
                    const double psi, const double i, 
                    const double t0, const double Pb0E, const double e0, const double l0, const double gamma0,
                    const double D_GW, const double RA_GW, const double DEC_GW, 
                    const double D_P,  const double RA_P,  const double DEC_P, 
                    const double z,
                    const ResidualsMethod residuals_method,
                    const ResidualsTerms residuals_terms,
                    const std::vector<double> ts);

std::vector<double> eccentric_waveform( const double M, const double q,
                    const double psi, const double i, 
                    const double t0, const double Pb0E, const double e0, const double l0, const double gamma0,
                    const double D_GW, const double RA_GW, const double DEC_GW, 
                    const double D_P,  const double RA_P,  const double DEC_P, 
                    const double z,
                    const ResidualsTerms residuals_terms,
                    const std::vector<double> ts);

std::vector<std::vector<double>> eccentric_residuals_and_waveform(const double M, const double q,
                    const double psi, const double i, 
                    const double t0, const double Pb0E, const double e0, const double l0, const double gamma0,
                    const double D_GW, const double RA_GW, const double DEC_GW, 
                    const double D_P,  const double RA_P,  const double DEC_P, 
                    const double z,
                    const ResidualsTerms residuals_terms,
                    const std::vector<double> ts);

std::vector<std::vector<double> > eccentric_residuals_px(const double M, const double q,
                                                        const double psi, const double i,
                                                        const double t0, const double n0, const double e0, const double l0, const double gamma0,
                                                        const double DGW, const double delay,
                                                        const double z,
                                                        const ResidualsMethod residuals_method,
                                                        const ResidualsTerms residuals_terms,
                                                        const std::vector<double> _ts);

std::vector<std::vector<double> > eccentric_waveform_px( const double M, const double q,
                                                        const double psi, const double i,
                                                        const double t0, const double n0, const double e0, const double l0, const double gamma0,
                                                        const double DGW, const double z,
                                                        const std::vector<double> _ts);

/* std::vector<double> eccentric_waveform_fn(const double M, const double q,
                                         const double psi, const double i,
                                         const double n, const double e, const double l, const double gamma,
                                         const double DGW);  */

BinaryState solve_orbit_equations(const BinaryMass &bin_mass, const BinaryState &bin_init, const double delay);

std::vector<double> antenna_pattern(const double RA_GW, const double DEC_GW, const double RA_P, const double DEC_P);

bool mergeq(const double M, const double q, 
            const double Pb0E, const double e0, 
            const double z,
            const double t0, const double max_toa);

std::vector<std::vector<double> > fe_stat_funcs(const double M, const double q,
                                              const double t0, const double Pb0E, const double e0, const double l0,
                                              const double D_GW, const double RA_GW, const double DEC_GW, 
                                              const double RA_P,  const double DEC_P, 
                                              const double z,
                                              const std::vector<double> _ts,
                                              const ResidualsMethod residuals_method);

#endif
