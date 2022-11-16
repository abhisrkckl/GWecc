#ifndef _eccentric_residuals_hpp_
#define _eccentric_residuals_hpp_ 1

#include <valarray>
#include <tuple>
#include <array>

#include "binary.hpp"

enum class ResidualsMethod {PC, Adb, Num, PM};
enum class ResidualsTerms {Earth=1, Pulsar=2, Both=3};

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
/*typedef Signal1D (*eccentric_residuals_func_t)(const BinaryMass &,
                                               const BinaryState &,
                                               const SkyPosition &,
                                               const SkyPosition &,
                                               const ResidualsTerms,
                                               const Signal1D &);*/
template <typename ResidualFunction>
double get_total_residual(const ResidualFunction& res, double t, double delay, ResidualsTerms residuals_terms){
    double R = 0;
    if(residuals_terms == ResidualsTerms::Earth || residuals_terms == ResidualsTerms::Both) R -= res(t);
    if(residuals_terms == ResidualsTerms::Pulsar|| residuals_terms == ResidualsTerms::Both) R += res(t-delay);
    return R;
}

template <typename WfResFunction>
std::array<double,2> get_total_residual_and_waveform(const WfResFunction& wfres, double t, double delay, ResidualsTerms residuals_terms){
    double h=0, R=0;
    if(residuals_terms == ResidualsTerms::Earth || residuals_terms == ResidualsTerms::Both){
        const auto [RE, hE] = wfres(t);
        R -= RE;
        h -= hE;
    }
    if(residuals_terms == ResidualsTerms::Pulsar|| residuals_terms == ResidualsTerms::Both){
        const auto [RP, hP] = wfres(t-delay);
        R += RP;
        h += hP;
    }
    return {R, h};
}

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
std::tuple<Signal1D,Signal1D> adiabatic_residuals_and_waveform(const BinaryMass &bin_mass,
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

Signal1D adiabatic_residuals(const BinaryMass &bin_mass,
                             const BinaryState &bin_init,
                             const SkyPosition &bin_pos,
                             const SkyPosition &psr_pos,
                             const ResidualsTerms residuals_terms,
                             const Signal1D &ts);

std::tuple<Signal1D, Signal1D> adiabatic_residuals_px(const BinaryMass &bin_mass,
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
