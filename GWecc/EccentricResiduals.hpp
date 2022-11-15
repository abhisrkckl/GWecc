#ifndef _eccentric_residuals_hpp_
#define _eccentric_residuals_hpp_ 1

#include <valarray>
#include <tuple>
#include "binary.hpp"
#include "residual_options.hpp"

typedef std::valarray<double> Signal1D;

/*
 * Computes H0
 */
double gw_amplitude(const BinaryMass &bin_mass,
                    const BinaryState &bin_init,
                    const double DGW);


/*
 * This is to abstract over different methods of computing RA and RB.
 */
//typedef std::tuple<Signal1D, Signal1D> (*)(const BinaryMass&, const BinaryState&, const Signal1D&);
typedef Signal1D (*eccentric_residuals_func_t) ( const BinaryMass &,
                                                const BinaryState &,
                                                const SkyPosition &,
                                                const SkyPosition &,
                                                const ResidualsTerms,
                                                const Signal1D &);

/*
 * Computes R(t)
 */
Signal1D eccentric_residuals(const BinaryMass &bin_mass,
                            const BinaryState &bin_init,
                            const SkyPosition &bin_pos,
                            const SkyPosition &psr_pos,
                            const ResidualsMethod residuals_method,
                            const ResidualsTerms residuals_terms,
                            const Signal1D &ts);

/*
 * Computes h(t)
 */
Signal1D eccentric_waveform( const BinaryMass &bin_mass,
                            const BinaryState &bin_init,
                            const SkyPosition &bin_pos,
                            const SkyPosition &psr_pos,
                            const ResidualsTerms residuals_terms,
                            const Signal1D &ts);

/*
 * Computes R(t) and h(t)
 */
std::tuple<Signal1D,Signal1D> eccentric_residuals_and_waveform(const BinaryMass &bin_mass,
                                                            const BinaryState &bin_init,
                                                            const SkyPosition &bin_pos,
                                                            const SkyPosition &psr_pos,
                                                            const ResidualsTerms residuals_terms,
                                                            const Signal1D &ts);

std::tuple<Signal1D, Signal1D> eccentric_residuals_px(const BinaryMass &bin_mass,
                                                     const BinaryState &bin_init,
                                                     const double DGW, const double delay,
                                                     const ResidualsMethod residuals_method,
                                                     const ResidualsTerms residuals_terms,
                                                     const Signal1D &ts);

std::tuple<Signal1D, Signal1D> eccentric_waveform_px( const BinaryMass &bin_mass,
                                                     const BinaryState &bin_init,
                                                     const double DGW,
                                                     const Signal1D &ts);

/*
 * Compute Rp and Rx using analytic expressions. Valid for low eccentricities (e<0.3).
 * Based on Boetzel et al. 2017
 */
Signal1D eccentric_residuals_Anl(const BinaryMass &bin_mass,
                                const BinaryState &bin_init,
                                const SkyPosition &bin_pos,
                                const SkyPosition &psr_pos,
                                const ResidualsTerms residuals_terms,
                                const Signal1D &ts);
                    
std::tuple<Signal1D, Signal1D> eccentric_residuals_px_Anl(const BinaryMass &bin_mass,
                                                         const BinaryState &bin_init,
                                                         const double DGW,
                                                         const Signal1D &ts);


Signal1D eccentric_residuals_Adb(const BinaryMass &bin_mass,
                                const BinaryState &bin_init,
                                const SkyPosition &bin_pos,
                                const SkyPosition &psr_pos,
                                const ResidualsTerms residuals_terms,
                                const Signal1D &ts);

std::tuple<Signal1D, Signal1D> eccentric_residuals_px_Adb(const BinaryMass &bin_mass,
                                                         const BinaryState &bin_init,
                                                         const double DGW,
                                                         const Signal1D &ts);

Signal1D eccentric_residuals_PM(const BinaryMass &bin_mass,
                               const BinaryState &bin_init,
                               const SkyPosition &bin_pos,
                               const SkyPosition &psr_pos,
                               const ResidualsTerms residuals_terms,
                               const Signal1D &ts);

std::tuple<Signal1D, Signal1D> eccentric_residuals_px_PM(const BinaryMass &bin_mass,
                                                        const BinaryState &bin_init,
                                                        const double DGW,
                                                        const Signal1D &ts);

/*
 * Compute Rp and Rx using FFT. Valid for all eccentricities.
 * Based on Tessmer & Gopakumar 2006
 */
/*
std::tuple<Signal1D, Signal1D> eccentric_residuals_px_FFT(const BinaryMass &bin_mass,
                                                         const BinaryState &bin_init,
                                                         const Signal1D &ts);
*/

/*
 * Compute Rp and Rx using numerical integration. Valid for all eccentricities.
 */
Signal1D eccentric_residuals_Num(const BinaryMass &bin_mass,
                                const BinaryState &bin_init,
                                const SkyPosition &bin_pos,
                                const SkyPosition &psr_pos,
                                const ResidualsTerms residuals_terms,
                                const Signal1D &ts);
                    
Signal1D eccentric_residuals_fn_Num(const BinaryMass &bin_mass,
                                   const BinaryState &bin_init,
                                   const double Fp, const double Fx, const double DGW,
                                   const Signal1D &ts);

std::tuple<Signal1D, Signal1D> eccentric_residuals_px_Num(const BinaryMass &bin_mass,
                                                         const BinaryState &bin_init,
                                                         const double DGW,
                                                         const Signal1D &ts);

double eccentric_residuals_fn_pt(double t, void *_params);

#endif
