#include "EccentricResiduals.hpp"
#include "OrbitalEvolution.hpp"
#include "AntennaPattern.hpp"
#include "PN.hpp"
#include <stdexcept>

auto choose_EccentricResiduals_fn(const ResidualsMethod residuals_method){
    auto fn = EccentricResiduals_Num;
    if(residuals_method == ResidualsMethod::Anl){
        fn = EccentricResiduals_Anl;
    }
    else if(residuals_method == ResidualsMethod::Adb){
        fn = EccentricResiduals_Adb;
    }
    return fn;
}

auto choose_EccentricResiduals_px_fn(const ResidualsMethod residuals_method){
    auto fn = EccentricResiduals_px_Num;
    if(residuals_method == ResidualsMethod::Anl){
        fn = EccentricResiduals_px_Anl;
    }
    else if(residuals_method == ResidualsMethod::Adb){
        fn = EccentricResiduals_px_Adb;
    }
    return fn;
}

Signal1D EccentricResiduals(const BinaryMass &bin_mass,
                            const BinaryState &bin_init,
                            const SkyPosition &bin_pos,
                            const SkyPosition &psr_pos,
                            const ResidualsMethod residuals_method,
                            const ResidualsTerms residuals_terms,
                            const Signal1D &ts){

    const auto EccentricResiduals_fn = choose_EccentricResiduals_fn(residuals_method);

    return EccentricResiduals_fn(bin_mass,
                                 bin_init,
                                 bin_pos,
                                 psr_pos,
                                 residuals_terms,
                                 ts);
}

std::tuple<Signal1D, Signal1D> EccentricResiduals_px(const BinaryMass &bin_mass,
                                                     const BinaryState &bin_init,
                                                     const double DGW, const double delay,
                                                     const ResidualsMethod residuals_method,
                                                     const ResidualsTerms residuals_terms,
                                                     const Signal1D &ts){

    
    const auto EccentricResiduals_px_fn = (residuals_method==ResidualsMethod::Anl) ? EccentricResiduals_px_Anl
                                                                                   : EccentricResiduals_px_Num;

    if(residuals_terms == ResidualsTerms::Earth){
        Signal1D RpE, RxE;
        std::tie(RpE,RxE) = EccentricResiduals_px_fn(bin_mass, bin_init, DGW, ts);
        //auto [RpE,RxE] = EccentricResiduals_px_fn(bin_mass, bin_init, DGW, ts);
        RpE *= -1;
        RxE *= -1;
        return std::make_tuple(RpE,RxE);
    }
    else if(residuals_terms == ResidualsTerms::Pulsar){    
        return EccentricResiduals_px_fn(bin_mass, bin_init, DGW, ts + delay);
    }
    else{
        Signal1D RpE, RxE, RpP, RxP;
        std::tie(RpE,RxE) = EccentricResiduals_px_fn(bin_mass, bin_init, DGW, ts);
        std::tie(RpP,RxP) = EccentricResiduals_px_fn(bin_mass, bin_init, DGW, ts + delay);
        //const auto [RpE,RxE] = EccentricResiduals_px_fn(bin_mass, bin_init, DGW, ts);
        //const auto [RpP,RxP] = EccentricResiduals_px_fn(bin_mass, bin_init, DGW, ts + delay);
        
        const auto Rp = RpP-RpE,
                   Rx = RxP-RxE;
        
        return std::make_tuple(Rp, Rx);
    }
}

double GWAmplitude(const BinaryMass &bin_mass,
                   const BinaryState &bin_init,
                   const double DGW){

    if(DGW<=0){
        throw std::invalid_argument("Invalid luminosity distance value in GWAmplitude.");
    }

    const auto eta = bin_mass.symmetric_mass_ratio(),
               x   = PN_param_x(bin_mass, bin_init);

    return bin_mass.mass()*eta/DGW * x;
}
