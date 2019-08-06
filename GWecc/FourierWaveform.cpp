#include "FourierWaveform.hpp"
#include <cmath>
#include "ipow.hpp"
#include "PN.hpp"

/*
Fourier2DCoeffs_t FourierWaveformCoeffs(const BinaryMass &binmass, const BinaryState &binstate, const double beta){
    
    / *
     * p => A
     * x => B
     * /
    
    Array_pq    ap, bp, cp, dp, 
            ax, bx, cx, dx;
    
    const double     C1beta = cos(beta),
            C2beta = cos(2*beta),
            C3beta = cos(3*beta),
            C4beta = cos(4*beta),
            S1beta = sin(beta),
            S2beta = sin(2*beta),
            S3beta = sin(3*beta),
            S4beta = sin(4*beta),
            Ci     = cos(binstate.i),
            Si     = sin(binstate.i);
    
    const double    e    = binstate.e,
            eta  = binmass.symmetric_mass_ratio(),
            delta= binmass.differential_mass_ratio(),
            xi   = PN_param_xi(binmass, binstate),
            k    = advance_of_periastron(binmass, binstate),
            x    = pow((1+k)*xi, 3./2);
    
    / *
     * These files provide the expressions for ap,...,dx.
     * They are created using Mathematica and a sed script.
     * Ideally not to be edited by humans.
     * 
     * These files refer to the following functions:
     *    sqrt, ipow
     * and the following variables:
     *    C1beta,...,C4beta
     *    S1beta,...,S4beta,
     *    Ci, Si,
     *    e, x, eta, delta
     * /
    #include "ap.incl"
    #include "bp.incl"
    #include "cp.incl"
    #include "dp.incl"
    #include "ax.incl"
    #include "bx.incl"
    #include "cx.incl"
    #include "dx.incl"
    
    return Fourier2DCoeffs_t{ ap, bp, cp, dp,
                  ax, bx, cx, dx };
}*/

Fourier2DCoeffs_t FourierResidualCoeffs(const BinaryMass &binmass, const BinaryState &binstate, const double beta){
    
    /*
     * p => A
     * x => B
     */
    
    const double C1beta = cos(beta),
                 C2beta = cos(2*beta),
                 C3beta = cos(3*beta),
                 C4beta = cos(4*beta),
                 S1beta = sin(beta),
                 S2beta = sin(2*beta),
                 S3beta = sin(3*beta),
                 S4beta = sin(4*beta),
                 Ci     = cos(binstate.i),
                 Si     = sin(binstate.i);
    
    const double &e    = binstate.e,
                 n     = binstate.n,
                 eta   = binmass.symmetric_mass_ratio(),
                 delta = binmass.differential_mass_ratio(),
                 xi    = PN_param_xi(binmass, binstate),
                 k     = advance_of_periastron(binmass, binstate),
                 x     = pow((1+k)*xi, 3./2);
    
    Array_pq Ap, Bp, Cp, Dp, 
             Ax, Bx, Cx, Dx;
    
    /* These files provide the expressions for ap,...,dx.
     * They are created using Mathematica and a sed script.
     * Ideally not to be edited by humans.
     * 
     * These files refer to the following functions:
     *    sqrt, ipow
     * and the following variables:
     *    C1beta,...,C4beta
     *    S1beta,...,S4beta,
     *    Ci, Si,
     *    e, x, eta, delta
     */
    #include "ResCoeffs.incl"
    
    return Fourier2DCoeffs_t{ Ap, Bp, Cp, Dp, 
                              Ax, Bx, Cx, Dx };
}

std::tuple<double,double> FourierResidual_pt(const Fourier2DCoeffs_t &rc, const double l, const double lambda){
    
    Array_pq CC, CS, SC, SS;
    
    std::array<double,qmax+1> cosqlambda, sinqlambda;
    for(size_t q=0; q<=qmax; q++){
        cosqlambda[q] = cos(q*lambda);
        sinqlambda[q] = sin(q*lambda);
    }
    for(size_t p=0; p<=pmax; p++){

        const double sinpl = sin(p*l),
                     cospl = cos(p*l);
            
        for(size_t q=0; q<=qmax; q++){
            
            CC(p,q) = cospl * cosqlambda[q],
            CS(p,q) = cospl * sinqlambda[q],
            SC(p,q) = sinpl * cosqlambda[q],
            SS(p,q) = sinpl * sinqlambda[q];
        }
    }
    
    const double RA_pt = (CC*rc.Ap + SC*rc.Bp + CS*rc.Cp + SS*rc.Dp).sum(), 
                 RB_pt = (CC*rc.Ax + SC*rc.Bx + CS*rc.Cx + SS*rc.Dx).sum();
            
    return std::make_tuple(RA_pt, RB_pt);
}

