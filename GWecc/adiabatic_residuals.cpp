#include <iostream>
#include <stdexcept>
#include "orbital_evolution.hpp"
#include "eccentric_residuals.hpp"
#include "antenna_pattern.hpp"
#include "mikkola.h"
#include "post_newtonian.hpp"
#include "ipow.hpp"
#include "waveform_vars.hpp"

auto adiabatic_residual_px(const BinaryMass &bin_mass,
                           const BinaryState &bin_init,
                           const double DGW,
                           const EvolveCoeffs& ev_coeffs,
                           const double t){

    const auto bin_now = solve_orbit_equations(bin_init, ev_coeffs, t-bin_init.t);

    if(bin_now.merged){
        throw std::runtime_error("The binary has already merged.");
    }
    
    const auto [scu, phi, w] = get_sincosu_phi_omega(bin_mass, bin_now);
    const auto A012 = get_residual_A012(bin_now.e, scu, w);
    const auto a012 = get_residual_a012(ev_coeffs.cosi);
    const auto H0 = gw_amplitude(bin_mass, bin_now, DGW);
    const auto S = H0/bin_now.n;
    const auto sAB = get_residual_sAB(A012, a012, S);
    const SinCos sc2psi{bin_now.psi};
    const auto spx = get_residual_spx(sAB, sc2psi);

    return spx;
}

auto adiabatic_residual_and_waveform_px(const BinaryMass &bin_mass,
                                        const BinaryState &bin_init,
                                        const double DGW,
                                        const EvolveCoeffs& ev_coeffs,
                                        const double t){

    const auto bin_now = solve_orbit_equations(bin_init, ev_coeffs, t-bin_init.t);

    if(bin_now.merged){
        throw std::runtime_error("The binary has already merged.");
    }

    const auto [scu, phi, w] = get_sincosu_phi_omega(bin_mass, bin_now);
    const auto A012 = get_residual_A012(bin_now.e, scu, w);
    const auto a012 = get_residual_a012(ev_coeffs.cosi);
    const auto H0 = gw_amplitude(bin_mass, bin_now, DGW);
    const auto S = H0/bin_now.n;
    const auto sAB = get_residual_sAB(A012, a012, S);
    const SinCos sc2psi{bin_now.psi};
    const auto spx = get_residual_spx(sAB, sc2psi);

    const auto hA012 = get_waveform_A012(bin_now.e, scu, phi);
    const auto hAB = get_residual_sAB(hA012, a012, H0);
    const auto hpx = get_residual_spx(hAB, sc2psi);

    const auto& [sp, sx] = spx;
    const auto& [hp, hx] = hpx;

    return std::make_tuple(sp, sx, hp, hx);
}

double adiabatic_residual(const BinaryMass &bin_mass,
                          const BinaryState &bin_init,
                          const AntennaPattern Fpx, const double DGW,
                          const EvolveCoeffs& ev_coeffs,
                          const double t) {
    const auto [sp, sx] = adiabatic_residual_px(bin_mass, bin_init, DGW, ev_coeffs, t);
    const auto [_, Fp, Fx] = Fpx;
    return Fp*sp + Fx*sx;
}

Signal1D adiabatic_residuals(const BinaryMass &bin_mass,
                             const BinaryState &bin_init,
                             const SkyPosition &bin_pos,
                             const SkyPosition &psr_pos,
                             const ResidualsTerms residuals_terms,
                             const Signal1D &ts){
    
    const auto Fpx = antenna_pattern(bin_pos, psr_pos);
    const auto DGW = bin_pos.DL;
    const auto ev_coeffs = compute_evolve_coeffs(bin_mass, bin_init);

    const auto res = [&](double t){
        return adiabatic_residual(bin_mass, bin_init, Fpx, DGW, ev_coeffs, t);
    };

    const auto nts = ts.size();
    Signal1D Rs(nts);
    for(unsigned idx=0; idx<nts; idx++){
        const auto& t = ts[idx];

        if(residuals_terms == ResidualsTerms::Earth){
            Rs[idx] = -res(t);
        }
        else{
            const auto delay = -psr_pos.DL*(1-Fpx.cosmu) / (1+bin_pos.z);
            if(residuals_terms == ResidualsTerms::Pulsar){
                Rs[idx] = res(t + delay);
            }
            else{
                Rs[idx] = res(t+delay) - res(t);
            }
        }
    }

    return Rs;
}

std::tuple<Signal1D,Signal1D> eccentric_residuals_and_waveform(const BinaryMass &bin_mass,
                                                            const BinaryState &bin_init,
                                                            const SkyPosition &bin_pos,
                                                            const SkyPosition &psr_pos,
                                                            const ResidualsTerms residuals_terms,
                                                            const Signal1D &ts){
    
    const auto [cosmu, Fp, Fx] = antenna_pattern(bin_pos, psr_pos);
    const auto DGW = bin_pos.DL;
    const auto ev_coeffs = compute_evolve_coeffs(bin_mass, bin_init);

    const auto wfres = [&](double t){
        return adiabatic_residual_and_waveform_px(bin_mass, bin_init, DGW, ev_coeffs, t);
    };

    const auto nts = ts.size();
    Signal1D Rs(nts), hs(nts);
    for(unsigned idx=0; idx<nts; idx++){
        const auto& t = ts[idx];

        if(residuals_terms == ResidualsTerms::Earth){
            const auto [sp, sx, hp, hx] = wfres(t);
            Rs[idx] = -(Fp*sp + Fx*sx);
            hs[idx] = -(Fp*hp + Fx*hx);
        }
        else{
            const auto delay = -psr_pos.DL*(1-cosmu) / (1+bin_pos.z);
            if(residuals_terms == ResidualsTerms::Pulsar){
                const auto [sp, sx, hp, hx] = wfres(t+delay);
                Rs[idx] = (Fp*sp + Fx*sx);
                hs[idx] = (Fp*hp + Fx*hx);
            }
            else{
                const auto [spP, sxP, hpP, hxP] = wfres(t+delay);
                const auto [spE, sxE, hpE, hxE] = wfres(t);
                Rs[idx] = Fp*(spP-spE) + Fx*(sxP-sxE);
                hs[idx] = Fp*(hpP-hpE) + Fx*(hxP-hxE);
            }
        }
    }

    return {Rs, hs};
}

std::tuple<Signal1D, Signal1D> adiabatic_residuals_px(const BinaryMass &bin_mass,
                                                      const BinaryState &bin_init,
                                                      const double DGW,
                                                      const Signal1D &ts){

    const auto ev_coeffs = compute_evolve_coeffs(bin_mass, bin_init);

    const auto nts = ts.size();
    Signal1D Rps(nts), Rxs(nts);
    for(unsigned idx=0; idx<nts; idx++){
        const auto [spE, sxE] = adiabatic_residual_px(bin_mass, bin_init, DGW, ev_coeffs, ts[idx]);
        Rps[idx] = spE;
        Rxs[idx] = sxE;
    }

    return std::make_tuple(Rps, Rxs);
}