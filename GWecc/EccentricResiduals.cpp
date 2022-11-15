#include <stdexcept>
#include "EccentricResiduals.hpp"
#include "OrbitalEvolution.hpp"
#include "antenna_pattern.hpp"
#include "post_newtonian.hpp"

auto choose_eccentric_residuals_fn(const ResidualsMethod residuals_method){
    switch (residuals_method){
        break; case ResidualsMethod::Anl: return eccentric_residuals_Anl;
        break; case ResidualsMethod::Adb: return eccentric_residuals_Adb;
        break; case ResidualsMethod::Num: return eccentric_residuals_Num;
        break; case ResidualsMethod::PM : return eccentric_residuals_PM;
    }
}

auto choose_eccentric_residuals_px_fn(const ResidualsMethod residuals_method){
    switch (residuals_method){
        break; case ResidualsMethod::Anl: return eccentric_residuals_px_Anl;
        break; case ResidualsMethod::Adb: return eccentric_residuals_px_Adb;
        break; case ResidualsMethod::Num: return eccentric_residuals_px_Num;
        break; case ResidualsMethod::PM : return eccentric_residuals_px_PM;
    }
}

Signal1D eccentric_residuals(const BinaryMass &bin_mass,
                            const BinaryState &bin_init,
                            const SkyPosition &bin_pos,
                            const SkyPosition &psr_pos,
                            const ResidualsMethod residuals_method,
                            const ResidualsTerms residuals_terms,
                            const Signal1D &ts){

    const auto eccentric_residuals_fn = choose_eccentric_residuals_fn(residuals_method);

    return eccentric_residuals_fn(bin_mass,
                                 bin_init,
                                 bin_pos,
                                 psr_pos,
                                 residuals_terms,
                                 ts);
}

std::tuple<Signal1D, Signal1D> eccentric_residuals_px(const BinaryMass &bin_mass,
                                                     const BinaryState &bin_init,
                                                     const double DGW, const double delay,
                                                     const ResidualsMethod residuals_method,
                                                     const ResidualsTerms residuals_terms,
                                                     const Signal1D &ts){

    
    const auto eccentric_residuals_px_fn = choose_eccentric_residuals_px_fn(residuals_method);

    if(residuals_terms == ResidualsTerms::Earth){
        const auto [RpE,RxE] = eccentric_residuals_px_fn(bin_mass, bin_init, DGW, ts);
        return std::make_tuple(-RpE, -RxE);
    }
    else if(residuals_terms == ResidualsTerms::Pulsar){    
        return eccentric_residuals_px_fn(bin_mass, bin_init, DGW, ts + delay);
    }
    else{
        const auto [RpE,RxE] = eccentric_residuals_px_fn(bin_mass, bin_init, DGW, ts);
        const auto [RpP,RxP] = eccentric_residuals_px_fn(bin_mass, bin_init, DGW, ts + delay);
        
        const auto Rp = RpP-RpE,
                   Rx = RxP-RxE;
        
        return std::make_tuple(Rp, Rx);
    }
}

double gw_amplitude(const BinaryMass &bin_mass,
                   const BinaryState &bin_init,
                   const double DGW){

    if(DGW<=0){
        throw std::invalid_argument("Invalid luminosity distance value in GWAmplitude.");
    }

    const auto eta = bin_mass.symmetric_mass_ratio(),
               x   = pn_param_x(bin_mass, bin_init);

    return bin_mass.mass()*eta/DGW * x;
}
