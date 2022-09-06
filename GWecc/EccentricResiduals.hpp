#ifndef _EccentricResiduals_hpp_
#define _EccentricResiduals_hpp_ 1

#include <valarray>
#include <tuple>
#include "Binary.hpp"
#include "ResidualsOptions.hpp"

typedef std::valarray<double> Signal1D;

/*
 * Computes H0
 */
double GWAmplitude(const BinaryMass &bin_mass,
                   const BinaryState &bin_init,
                   const double DGW);


/*
 * This is to abstract over different methods of computing RA and RB.
 */
//typedef std::tuple<Signal1D, Signal1D> (*)(const BinaryMass&, const BinaryState&, const Signal1D&);
typedef Signal1D (*EccentricResiduals_func_t) ( const BinaryMass &,
                                                const BinaryState &,
                                                const SkyPosition &,
                                                const SkyPosition &,
                                                const ResidualsTerms,
                                                const Signal1D &);

/*
 * Computes R(t)
 */
Signal1D EccentricResiduals(const BinaryMass &bin_mass,
                            const BinaryState &bin_init,
                            const SkyPosition &bin_pos,
                            const SkyPosition &psr_pos,
                            const ResidualsMethod residuals_method,
                            const ResidualsTerms residuals_terms,
                            const Signal1D &ts);

/*
 * Computes h(t)
 */
Signal1D EccentricWaveform( const BinaryMass &bin_mass,
                            const BinaryState &bin_init,
                            const SkyPosition &bin_pos,
                            const SkyPosition &psr_pos,
                            const ResidualsTerms residuals_terms,
                            const Signal1D &ts);

std::tuple<Signal1D, Signal1D> EccentricResiduals_px(const BinaryMass &bin_mass,
                                                     const BinaryState &bin_init,
                                                     const double DGW, const double delay,
                                                     const ResidualsMethod residuals_method,
                                                     const ResidualsTerms residuals_terms,
                                                     const Signal1D &ts);

std::tuple<Signal1D, Signal1D> EccentricWaveform_px( const BinaryMass &bin_mass,
                                                     const BinaryState &bin_init,
                                                     const double DGW,
                                                     const Signal1D &ts);

/*
 * Compute Rp and Rx using analytic expressions. Valid for low eccentricities (e<0.3).
 * Based on Boetzel et al. 2017
 */
Signal1D EccentricResiduals_Anl(const BinaryMass &bin_mass,
                                const BinaryState &bin_init,
                                const SkyPosition &bin_pos,
                                const SkyPosition &psr_pos,
                                const ResidualsTerms residuals_terms,
                                const Signal1D &ts);
                    
std::tuple<Signal1D, Signal1D> EccentricResiduals_px_Anl(const BinaryMass &bin_mass,
                                                         const BinaryState &bin_init,
                                                         const double DGW,
                                                         const Signal1D &ts);


Signal1D EccentricResiduals_Adb(const BinaryMass &bin_mass,
                                const BinaryState &bin_init,
                                const SkyPosition &bin_pos,
                                const SkyPosition &psr_pos,
                                const ResidualsTerms residuals_terms,
                                const Signal1D &ts);

std::tuple<Signal1D, Signal1D> EccentricResiduals_px_Adb(const BinaryMass &bin_mass,
                                                         const BinaryState &bin_init,
                                                         const double DGW,
                                                         const Signal1D &ts);

Signal1D EccentricResiduals_PM(const BinaryMass &bin_mass,
                               const BinaryState &bin_init,
                               const SkyPosition &bin_pos,
                               const SkyPosition &psr_pos,
                               const ResidualsTerms residuals_terms,
                               const Signal1D &ts);

std::tuple<Signal1D, Signal1D> EccentricResiduals_px_PM(const BinaryMass &bin_mass,
                                                        const BinaryState &bin_init,
                                                        const double DGW,
                                                        const Signal1D &ts);

/*
 * Compute Rp and Rx using FFT. Valid for all eccentricities.
 * Based on Tessmer & Gopakumar 2006
 */
/*
std::tuple<Signal1D, Signal1D> EccentricResiduals_px_FFT(const BinaryMass &bin_mass,
                                                         const BinaryState &bin_init,
                                                         const Signal1D &ts);
*/

/*
 * Compute Rp and Rx using numerical integration. Valid for all eccentricities.
 */
Signal1D EccentricResiduals_Num(const BinaryMass &bin_mass,
                                const BinaryState &bin_init,
                                const SkyPosition &bin_pos,
                                const SkyPosition &psr_pos,
                                const ResidualsTerms residuals_terms,
                                const Signal1D &ts);
                    
Signal1D EccentricResiduals_fn_Num(const BinaryMass &bin_mass,
                                   const BinaryState &bin_init,
                                   const double Fp, const double Fx, const double DGW,
                                   const Signal1D &ts);

std::tuple<Signal1D, Signal1D> EccentricResiduals_px_Num(const BinaryMass &bin_mass,
                                                         const BinaryState &bin_init,
                                                         const double DGW,
                                                         const Signal1D &ts);

double EccentricResiduals_fn_pt(double t, void *_params);

#endif
