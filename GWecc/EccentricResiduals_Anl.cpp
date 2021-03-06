#include "EccentricResiduals.hpp"
#include "OrbitalEvolution.hpp"
#include "AntennaPattern.hpp"
#include "FourierWaveform.hpp"
#include "PN.hpp"
#include <iostream>

Signal1D EccentricResiduals_Anl(const BinaryMass &bin_mass,
                                const BinaryState &bin_init,
                                const SkyPosition &bin_pos,
                                const SkyPosition &psr_pos,
                                const ResidualsTerms residuals_terms,
                                const Signal1D &ts){

    double cosmu, Fp, Fx;
    std::tie(cosmu, Fp, Fx) = AntennaPattern(bin_pos, psr_pos);
    //const auto [cosmu, Fp, Fx] = AntennaPattern(bin_pos, psr_pos);
    
    if(residuals_terms==ResidualsTerms::Earth){
        
        Signal1D RpE, RxE;
        std::tie(RpE, RxE) = EccentricResiduals_px_Anl(bin_mass, bin_init, bin_pos.DL, ts);
        //const auto [RpE, RxE] = EccentricResiduals_px_Anl(bin_mass, bin_init, bin_pos.DL, ts);
        
        return -(Fp*RpE + Fx*RxE);
    }
    else if(residuals_terms==ResidualsTerms::Pulsar){
        
        const auto delay = -psr_pos.DL*(1-cosmu) / (1+bin_pos.z);
        //const auto bin_psrterm = solve_orbit_equations(bin_mass, bin_init, delay);
        
        Signal1D RpP, RxP;
        std::tie(RpP, RxP) = EccentricResiduals_px_Anl(bin_mass, bin_init, bin_pos.DL, ts + delay);
        //const auto [RpP, RxP] = EccentricResiduals_px_Anl(bin_mass, bin_init, bin_pos.DL, ts + delay);
        
        return (Fp*RpP + Fx*RxP);
    }
    else{
        const auto delay = -psr_pos.DL*(1-cosmu) / (1+bin_pos.z);
        //const auto bin_psrterm = solve_orbit_equations(bin_mass, bin_init, delay);
        
        Signal1D RpP, RxP, RpE, RxE;
        std::tie(RpP, RxP) = EccentricResiduals_px_Anl(bin_mass, bin_init, bin_pos.DL, ts + delay);
        std::tie(RpE, RxE) = EccentricResiduals_px_Anl(bin_mass, bin_init, bin_pos.DL, ts);
        //const auto [RpP, RxP] = EccentricResiduals_px_Anl(bin_mass, bin_init, bin_pos.DL, ts + delay);
        //const auto [RpE, RxE] = EccentricResiduals_px_Anl(bin_mass, bin_init, bin_pos.DL, ts);
        
        return    (Fp*RpP + Fx*RxP) 
                - (Fp*RpE + Fx*RxE);
    }
}

std::tuple<Signal1D, Signal1D> EccentricResiduals_px_Anl(const BinaryMass &bin_mass,
                                                         const BinaryState &bin_init,
                                                         const double DGW,
                                                         const Signal1D &ts){
    
    const auto length = ts.size();
    
    const EvolveCoeffs_t ev_coeffs = compute_evolve_coeffs(bin_mass, bin_init);
    
    const auto cos2Omega = cos(2*bin_init.Omega),
               sin2Omega = sin(2*bin_init.Omega);
    
    Signal1D Rp(length), Rx(length);
    
    //const auto    l0    = bin_init.l               -       bin_init.n*bin_init.t,
    //              lambda0    = bin_init.l + bin_init.gamma - (1+k)*bin_init.n*bin_init.t;
    
    for(size_t i=0; i<length; i++){
    
        const auto bin_now = solve_orbit_equations(bin_init, ev_coeffs, ts[i]-bin_init.t);
        
        //const auto k = advance_of_periastron(bin_mass, bin_now);
    
        const auto res_coeffs  = FourierResidualCoeffs(bin_mass, bin_now, 0);
        
        const auto //nt     = bin_now.n,                //bin_init.n * (ts[i]-bin_init.t),
                   l      = bin_now.l,                  //      nt + bin_init.l,
                   lambda = l+bin_now.gamma;            //(1+k)*nt + bin_init.l + bin_init.gamma;
        
        double RA, RB;
        std::tie(RA, RB) = FourierResidual_pt(res_coeffs, l, lambda);
        //const auto [RA, RB] = FourierResidual_pt(res_coeffs, l, lambda);

        const auto H0 = GWAmplitude(bin_mass, bin_now, DGW);
        
        Rp[i] = H0*(cos2Omega*RA - sin2Omega*RB);
        Rx[i] = H0*(cos2Omega*RB + sin2Omega*RA);
    }
    
    return std::make_tuple(Rp, Rx);
}

