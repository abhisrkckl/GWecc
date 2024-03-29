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
    /* Begin auto-generated code */
    Ap(0,0) = 0;
    Bp(0,0) = 0;
    Cp(0,0) = 0;
    Dp(0,0) = 0;
    Ax(0,0) = 0;
    Bx(0,0) = 0;
    Cx(0,0) = 0;
    Dx(0,0) = 0;
    Ap(0,1) = -((-737280 + 884736*ipow(e,2) - 117504*ipow(e,4) + 10496*ipow(e,6) + 4791*ipow(e,8) + ipow(Ci,2)*(-147456 - 294912*ipow(e,2) + 94464*ipow(e,4) + 9472*ipow(e,6) + 7307*ipow(e,8)))*S1beta*Si*sqrt(x)*delta)/(589824.*(1 + k)*n);
    Bp(0,1) = 0;
    Cp(0,1) = (C1beta*(-737280 + 884736*ipow(e,2) - 117504*ipow(e,4) + 10496*ipow(e,6) + 4791*ipow(e,8) + ipow(Ci,2)*(-147456 - 294912*ipow(e,2) + 94464*ipow(e,4) + 9472*ipow(e,6) + 7307*ipow(e,8)))*Si*sqrt(x)*delta)/(589824.*(1 + k)*n);
    Dp(0,1) = 0;
    Ax(0,1) = -((C1beta*Ci*(-1.5 + ipow(e,2) - (5*ipow(e,4))/128. + (13*ipow(e,6))/384. + (6049*ipow(e,8))/294912.)*Si*sqrt(x)*delta)/((1 + k)*n));
    Bx(0,1) = 0;
    Cx(0,1) = -(Ci*(-442368 + 294912*ipow(e,2) - 11520*ipow(e,4) + 9984*ipow(e,6) + 6049*ipow(e,8))*S1beta*Si*sqrt(x)*delta)/(294912.*(1 + k)*n);
    Dx(0,1) = 0;
    Ap(0,2) = -((1 + ipow(Ci,2))*(-2 + 5*ipow(e,2) - (23*ipow(e,4))/8. + (65*ipow(e,6))/144. - (85*ipow(e,8))/1152.)*S2beta + S2beta*x*(2.6666666666666665 + (83*ipow(e,4)*(1 - 3*eta))/6. - 8*eta + (40*ipow(e,2)*(-1 + 3*eta))/3. + (110*ipow(e,6)*(-1 + 3*eta))/27. - (455*ipow(e,8)*(-1 + 3*eta))/864. + ((1 + ipow(Ci,2))*(-1152*ipow(e,2)*(-315 + 151*eta + 4*ipow(Si,2)*(-1 + 3*eta)) - 56*ipow(e,6)*(-909 + 925*eta + 80*ipow(Si,2)*(-1 + 3*eta)) + 144*ipow(e,4)*(-1101 + 1309*eta + 86*ipow(Si,2)*(-1 + 3*eta)) + ipow(e,8)*(13176 + 8183*eta + 614*ipow(Si,2)*(-1 + 3*eta)) - 2304*(-9 - 11*eta + ipow(Si,2)*(-2 + 6*eta))))/6912.))/(2.*(1 + k)*n);
    Bp(0,2) = 0;
    Cp(0,2) = (C2beta*(1 + ipow(Ci,2))*(-2 + 5*ipow(e,2) - (23*ipow(e,4))/8. + (65*ipow(e,6))/144. - (85*ipow(e,8))/1152.) + C2beta*x*(2.6666666666666665 + (83*ipow(e,4)*(1 - 3*eta))/6. - 8*eta + (40*ipow(e,2)*(-1 + 3*eta))/3. + (110*ipow(e,6)*(-1 + 3*eta))/27. - (455*ipow(e,8)*(-1 + 3*eta))/864. + ((1 + ipow(Ci,2))*(-1152*ipow(e,2)*(-315 + 151*eta + 4*ipow(Si,2)*(-1 + 3*eta)) - 56*ipow(e,6)*(-909 + 925*eta + 80*ipow(Si,2)*(-1 + 3*eta)) + 144*ipow(e,4)*(-1101 + 1309*eta + 86*ipow(Si,2)*(-1 + 3*eta)) + ipow(e,8)*(13176 + 8183*eta + 614*ipow(Si,2)*(-1 + 3*eta)) - 2304*(-9 - 11*eta + ipow(Si,2)*(-2 + 6*eta))))/6912.))/(2.*(1 + k)*n);
    Dp(0,2) = 0;
    Ax(0,2) = (C2beta*Ci*(2304*(6 + x*(-13 + eta + 4*ipow(Si,2)*(-1 + 3*eta))) - 1152*ipow(e,2)*(30 + x*(275 - 31*eta + 16*ipow(Si,2)*(-1 + 3*eta))) + 144*ipow(e,4)*(138 + x*(769 - 313*eta + 80*ipow(Si,2)*(-1 + 3*eta))) + ipow(e,8)*(510 + x*(-14996 - 2723*eta + 296*ipow(Si,2)*(-1 + 3*eta))) - 8*ipow(e,6)*(390 + x*(4603 - 1195*eta + 320*ipow(Si,2)*(-1 + 3*eta)))))/(6912.*(1 + k)*n);
    Bx(0,2) = 0;
    Cx(0,2) = (Ci*S2beta*(2304*(6 + x*(-13 + eta + 4*ipow(Si,2)*(-1 + 3*eta))) - 1152*ipow(e,2)*(30 + x*(275 - 31*eta + 16*ipow(Si,2)*(-1 + 3*eta))) + 144*ipow(e,4)*(138 + x*(769 - 313*eta + 80*ipow(Si,2)*(-1 + 3*eta))) + ipow(e,8)*(510 + x*(-14996 - 2723*eta + 296*ipow(Si,2)*(-1 + 3*eta))) - 8*ipow(e,6)*(390 + x*(4603 - 1195*eta + 320*ipow(Si,2)*(-1 + 3*eta)))))/(6912.*(1 + k)*n);
    Dx(0,2) = 0;
    Ap(0,3) = (-3*(1 + ipow(Ci,2))*(16384 - 98304*ipow(e,2) + 151296*ipow(e,4) - 89344*ipow(e,6) + 23997*ipow(e,8))*S3beta*Si*sqrt(x)*delta)/(65536.*(1 + k)*n);
    Bp(0,3) = 0;
    Cp(0,3) = (3*C3beta*(1 + ipow(Ci,2))*(16384 - 98304*ipow(e,2) + 151296*ipow(e,4) - 89344*ipow(e,6) + 23997*ipow(e,8))*Si*sqrt(x)*delta)/(65536.*(1 + k)*n);
    Dp(0,3) = 0;
    Ax(0,3) = (-3*C3beta*Ci*(16384 - 98304*ipow(e,2) + 151296*ipow(e,4) - 89344*ipow(e,6) + 23997*ipow(e,8))*Si*sqrt(x)*delta)/(32768.*(1 + k)*n);
    Bx(0,3) = 0;
    Cx(0,3) = (-3*Ci*(16384 - 98304*ipow(e,2) + 151296*ipow(e,4) - 89344*ipow(e,6) + 23997*ipow(e,8))*S3beta*Si*sqrt(x)*delta)/(32768.*(1 + k)*n);
    Dx(0,3) = 0;
    Ap(0,4) = -((1 + ipow(Ci,2))*(2304 - 25344*ipow(e,2) + 72864*ipow(e,4) - 88576*ipow(e,6) + 53663*ipow(e,8))*S4beta*ipow(Si,2)*x*(-1 + 3*eta))/(3456.*(1 + k)*n);
    Bp(0,4) = 0;
    Cp(0,4) = (C4beta*(1 + ipow(Ci,2))*(2304 - 25344*ipow(e,2) + 72864*ipow(e,4) - 88576*ipow(e,6) + 53663*ipow(e,8))*ipow(Si,2)*x*(-1 + 3*eta))/(3456.*(1 + k)*n);
    Dp(0,4) = 0;
    Ax(0,4) = -(C4beta*Ci*(2304 - 25344*ipow(e,2) + 72864*ipow(e,4) - 88576*ipow(e,6) + 53663*ipow(e,8))*ipow(Si,2)*x*(-1 + 3*eta))/(1728.*(1 + k)*n);
    Bx(0,4) = 0;
    Cx(0,4) = -(Ci*(2304 - 25344*ipow(e,2) + 72864*ipow(e,4) - 88576*ipow(e,6) + 53663*ipow(e,8))*S4beta*ipow(Si,2)*x*(-1 + 3*eta))/(1728.*(1 + k)*n);
    Dx(0,4) = 0;
    Ap(1,0) = 0;
    Bp(1,0) = (e*(-4*(-1 + ipow(Ci,2))*x*(9216*(-27 + eta) - 2256*ipow(e,4)*(45 + eta) + 1152*ipow(e,2)*(-33 + 23*eta) + ipow(e,6)*(-97209 + 71*eta)) - 3*ipow(Si,2)*(9216*(-8 + (1 + ipow(Ci,2))*x*(-1 + 3*eta)) + 1152*ipow(e,2)*(8 + 9*(1 + ipow(Ci,2))*x*(-1 + 3*eta)) - 48*ipow(e,4)*(8 + 19*(1 + ipow(Ci,2))*x*(-1 + 3*eta)) + ipow(e,6)*(8 + 29*(1 + ipow(Ci,2))*x*(-1 + 3*eta)))))/(221184.*n);
    Cp(1,0) = 0;
    Dp(1,0) = 0;
    Ax(1,0) = -(Ci*e*(-9216 + 5760*ipow(e,2) + 528*ipow(e,4) + 457*ipow(e,6))*ipow(Si,2)*x*(-1 + 3*eta))/(18432.*n);
    Bx(1,0) = 0;
    Cx(1,0) = 0;
    Dx(1,0) = 0;
    Ap(1,1) = (e*(-9*(-96 + 128*ipow(e,2) - 38*ipow(e,4) + 5*ipow(e,6)) + ipow(Ci,2)*(-288 + 576*ipow(e,2) - 210*ipow(e,4) + 23*ipow(e,6)))*S1beta*Si*sqrt(x)*delta)/(288.*(2 + k)*n);
    Bp(1,1) = (C1beta*e*(ipow(Ci,2)*(288 - 576*ipow(e,2) + 210*ipow(e,4) - 23*ipow(e,6)) + 9*(-96 + 128*ipow(e,2) - 38*ipow(e,4) + 5*ipow(e,6)))*Si*sqrt(x)*delta)/(288.*(2 + k)*n);
    Cp(1,1) = (C1beta*e*(ipow(Ci,2)*(288 - 576*ipow(e,2) + 210*ipow(e,4) - 23*ipow(e,6)) + 9*(-96 + 128*ipow(e,2) - 38*ipow(e,4) + 5*ipow(e,6)))*Si*sqrt(x)*delta)/(288.*(2 + k)*n);
    Dp(1,1) = (e*(ipow(Ci,2)*(288 - 576*ipow(e,2) + 210*ipow(e,4) - 23*ipow(e,6)) + 9*(-96 + 128*ipow(e,2) - 38*ipow(e,4) + 5*ipow(e,6)))*S1beta*Si*sqrt(x)*delta)/(288.*(2 + k)*n);
    Ax(1,1) = (C1beta*Ci*e*(288 - 288*ipow(e,2) + 66*ipow(e,4) - 11*ipow(e,6))*Si*sqrt(x)*delta)/(144.*(2 + k)*n);
    Bx(1,1) = (Ci*e*(288 - 288*ipow(e,2) + 66*ipow(e,4) - 11*ipow(e,6))*S1beta*Si*sqrt(x)*delta)/(144.*(2 + k)*n);
    Cx(1,1) = (Ci*e*(288 - 288*ipow(e,2) + 66*ipow(e,4) - 11*ipow(e,6))*S1beta*Si*sqrt(x)*delta)/(144.*(2 + k)*n);
    Dx(1,1) = (C1beta*Ci*e*(-288 + 288*ipow(e,2) - 66*ipow(e,4) + 11*ipow(e,6))*Si*sqrt(x)*delta)/(144.*(2 + k)*n);
    Ap(1,2) = -((-((1 + ipow(Ci,2))*(-6*e + (23*ipow(e,3))/2. - (721*ipow(e,5))/96. + (1645*ipow(e,7))/768.)*S2beta) - S2beta*x*(e*(14.666666666666666 - 44*eta) + (203*ipow(e,3)*(-1 + 3*eta))/4. - (472*ipow(e,5)*(-1 + 3*eta))/9. + (296965*ipow(e,7)*(-1 + 3*eta))/13824. - ((1 + ipow(Ci,2))*(9216*e*(174 - 144*eta + ipow(Si,2)*(-1 + 3*eta)) - 48*ipow(e,5)*(-96306 + 113908*eta + 8827*ipow(Si,2)*(-1 + 3*eta)) + 1152*ipow(e,3)*(249*ipow(Si,2)*(-1 + 3*eta) + 266*(-24 + 17*eta)) + ipow(e,7)*(194177*ipow(Si,2)*(-1 + 3*eta) + 6*(-293792 + 371135*eta))))/55296.) + 2*(1 + k)*((1 + ipow(Ci,2))*(-3*e + (79*ipow(e,3))/8. - (1447*ipow(e,5))/192. + (6353*ipow(e,7))/3072.)*S2beta + S2beta*x*(e*(12.333333333333334 - 37*eta) + (377*ipow(e,3)*(-1 + 3*eta))/8. - (29813*ipow(e,5)*(-1 + 3*eta))/576. + (593611*ipow(e,7)*(-1 + 3*eta))/27648. - ((1 + ipow(Ci,2))*(1152*ipow(e,3)*(-5847 + 4109*eta + 237*ipow(Si,2)*(-1 + 3*eta)) - 48*ipow(e,5)*(-99147 + 113081*eta + 8669*ipow(Si,2)*(-1 + 3*eta)) - 9216*e*(ipow(Si,2)*(-1 + 3*eta) + 3*(-5 + 39*eta)) + ipow(e,7)*(195109*ipow(Si,2)*(-1 + 3*eta) + 3*(-529649 + 733091*eta))))/55296.)))/((1 + 2*k)*(3 + 2*k)*n));
    Bp(1,2) = -((C2beta*(1 + ipow(Ci,2))*(-3*e + (79*ipow(e,3))/8. - (1447*ipow(e,5))/192. + (6353*ipow(e,7))/3072.) + C2beta*x*(e*(12.333333333333334 - 37*eta) + (377*ipow(e,3)*(-1 + 3*eta))/8. - (29813*ipow(e,5)*(-1 + 3*eta))/576. + (593611*ipow(e,7)*(-1 + 3*eta))/27648. - ((1 + ipow(Ci,2))*(1152*ipow(e,3)*(-5847 + 4109*eta + 237*ipow(Si,2)*(-1 + 3*eta)) - 48*ipow(e,5)*(-99147 + 113081*eta + 8669*ipow(Si,2)*(-1 + 3*eta)) - 9216*e*(ipow(Si,2)*(-1 + 3*eta) + 3*(-5 + 39*eta)) + ipow(e,7)*(195109*ipow(Si,2)*(-1 + 3*eta) + 3*(-529649 + 733091*eta))))/55296.) + 2*(1 + k)*(C2beta*(1 + ipow(Ci,2))*(6*e - (23*ipow(e,3))/2. + (721*ipow(e,5))/96. - (1645*ipow(e,7))/768.) + C2beta*x*((44*e*(-1 + 3*eta))/3. - (203*ipow(e,3)*(-1 + 3*eta))/4. + (472*ipow(e,5)*(-1 + 3*eta))/9. - (296965*ipow(e,7)*(-1 + 3*eta))/13824. + ((1 + ipow(Ci,2))*(9216*e*(174 - 144*eta + ipow(Si,2)*(-1 + 3*eta)) - 48*ipow(e,5)*(-96306 + 113908*eta + 8827*ipow(Si,2)*(-1 + 3*eta)) + 1152*ipow(e,3)*(249*ipow(Si,2)*(-1 + 3*eta) + 266*(-24 + 17*eta)) + ipow(e,7)*(194177*ipow(Si,2)*(-1 + 3*eta) + 6*(-293792 + 371135*eta))))/55296.)))/((1 + 2*k)*(3 + 2*k)*n));
    Cp(1,2) = -(((C2beta*(1 + ipow(Ci,2))*e*(-4608 + 8832*ipow(e,2) - 5768*ipow(e,4) + 1645*ipow(e,6)))/768. - C2beta*x*((44*e*(-1 + 3*eta))/3. - (203*ipow(e,3)*(-1 + 3*eta))/4. + (472*ipow(e,5)*(-1 + 3*eta))/9. - (296965*ipow(e,7)*(-1 + 3*eta))/13824. + ((1 + ipow(Ci,2))*(9216*e*(174 - 144*eta + ipow(Si,2)*(-1 + 3*eta)) - 48*ipow(e,5)*(-96306 + 113908*eta + 8827*ipow(Si,2)*(-1 + 3*eta)) + 1152*ipow(e,3)*(249*ipow(Si,2)*(-1 + 3*eta) + 266*(-24 + 17*eta)) + ipow(e,7)*(194177*ipow(Si,2)*(-1 + 3*eta) + 6*(-293792 + 371135*eta))))/55296.) - 2*(1 + k)*(C2beta*(1 + ipow(Ci,2))*(-3*e + (79*ipow(e,3))/8. - (1447*ipow(e,5))/192. + (6353*ipow(e,7))/3072.) + C2beta*x*(e*(12.333333333333334 - 37*eta) + (377*ipow(e,3)*(-1 + 3*eta))/8. - (29813*ipow(e,5)*(-1 + 3*eta))/576. + (593611*ipow(e,7)*(-1 + 3*eta))/27648. - ((1 + ipow(Ci,2))*(1152*ipow(e,3)*(-5847 + 4109*eta + 237*ipow(Si,2)*(-1 + 3*eta)) - 48*ipow(e,5)*(-99147 + 113081*eta + 8669*ipow(Si,2)*(-1 + 3*eta)) - 9216*e*(ipow(Si,2)*(-1 + 3*eta) + 3*(-5 + 39*eta)) + ipow(e,7)*(195109*ipow(Si,2)*(-1 + 3*eta) + 3*(-529649 + 733091*eta))))/55296.)))/((1 + 2*k)*(3 + 2*k)*n));
    Dp(1,2) = -(((1 + ipow(Ci,2))*(-3*e + (79*ipow(e,3))/8. - (1447*ipow(e,5))/192. + (6353*ipow(e,7))/3072.)*S2beta + S2beta*x*(e*(12.333333333333334 - 37*eta) + (377*ipow(e,3)*(-1 + 3*eta))/8. - (29813*ipow(e,5)*(-1 + 3*eta))/576. + (593611*ipow(e,7)*(-1 + 3*eta))/27648. - ((1 + ipow(Ci,2))*(1152*ipow(e,3)*(-5847 + 4109*eta + 237*ipow(Si,2)*(-1 + 3*eta)) - 48*ipow(e,5)*(-99147 + 113081*eta + 8669*ipow(Si,2)*(-1 + 3*eta)) - 9216*e*(ipow(Si,2)*(-1 + 3*eta) + 3*(-5 + 39*eta)) + ipow(e,7)*(195109*ipow(Si,2)*(-1 + 3*eta) + 3*(-529649 + 733091*eta))))/55296.) - 2*(1 + k)*((1 + ipow(Ci,2))*(-6*e + (23*ipow(e,3))/2. - (721*ipow(e,5))/96. + (1645*ipow(e,7))/768.)*S2beta + S2beta*x*(e*(14.666666666666666 - 44*eta) + (203*ipow(e,3)*(-1 + 3*eta))/4. - (472*ipow(e,5)*(-1 + 3*eta))/9. + (296965*ipow(e,7)*(-1 + 3*eta))/13824. - ((1 + ipow(Ci,2))*(9216*e*(174 - 144*eta + ipow(Si,2)*(-1 + 3*eta)) - 48*ipow(e,5)*(-96306 + 113908*eta + 8827*ipow(Si,2)*(-1 + 3*eta)) + 1152*ipow(e,3)*(249*ipow(Si,2)*(-1 + 3*eta) + 266*(-24 + 17*eta)) + ipow(e,7)*(194177*ipow(Si,2)*(-1 + 3*eta) + 6*(-293792 + 371135*eta))))/55296.)))/((1 + 2*k)*(3 + 2*k)*n));
    Ax(1,2) = -(C2beta*Ci*e*(ipow(e,6)*(110268 + 228708*k + 15*x*(54790 - 26124*eta + 6707*ipow(Si,2)*(-1 + 3*eta)) + k*x*(1990672 - 836880*eta + 203393*ipow(Si,2)*(-1 + 3*eta))) - 9216*(k*(36 + 35*ipow(Si,2)*x*(-1 + 3*eta) - 4*x*(11 + 3*eta)) + 6*x*(-29 + ipow(Si,2)*(-2 + 6*eta))) + 1152*ipow(e,2)*(396 + 3*x*(1422 - 188*eta + 99*ipow(Si,2)*(-1 + 3*eta)) + k*(948 + x*(9432 - 1432*eta + 657*ipow(Si,2)*(-1 + 3*eta)))) - 48*ipow(e,4)*(8712 + 6*x*(12095 - 4000*eta + 1033*ipow(Si,2)*(-1 + 3*eta)) + k*(17364 + x*(138668 - 47284*eta + 12475*ipow(Si,2)*(-1 + 3*eta))))))/(27648.*(1 + 2*k)*(3 + 2*k)*n);
    Bx(1,2) = (Ci*e*S2beta*(9216*(3*(36 + x*(188 - 12*eta + 19*ipow(Si,2)*(-1 + 3*eta))) + 4*k*(36 + x*(130 - 12*eta + 23*ipow(Si,2)*(-1 + 3*eta)))) - 1152*ipow(e,2)*(8*k*(276 + x*(2583 - 434*eta + 180*ipow(Si,2)*(-1 + 3*eta))) + 3*(420 + x*(3744 - 680*eta + 261*ipow(Si,2)*(-1 + 3*eta)))) + 48*ipow(e,4)*(3*(5748 + x*(41908 - 15284*eta + 4211*ipow(Si,2)*(-1 + 3*eta))) + 4*k*(8652 + x*(66098 - 23284*eta + 6277*ipow(Si,2)*(-1 + 3*eta)))) - ipow(e,6)*(8*k*(59220 + x*(584411 - 222510*eta + 51394*ipow(Si,2)*(-1 + 3*eta))) + 3*(81684 + x*(894872 - 314400*eta + 69253*ipow(Si,2)*(-1 + 3*eta))))))/(55296.*(1 + 2*k)*(3 + 2*k)*n);
    Cx(1,2) = -(Ci*e*S2beta*(ipow(e,6)*(110268 + 228708*k + 15*x*(54790 - 26124*eta + 6707*ipow(Si,2)*(-1 + 3*eta)) + k*x*(1990672 - 836880*eta + 203393*ipow(Si,2)*(-1 + 3*eta))) - 9216*(k*(36 + 35*ipow(Si,2)*x*(-1 + 3*eta) - 4*x*(11 + 3*eta)) + 6*x*(-29 + ipow(Si,2)*(-2 + 6*eta))) + 1152*ipow(e,2)*(396 + 3*x*(1422 - 188*eta + 99*ipow(Si,2)*(-1 + 3*eta)) + k*(948 + x*(9432 - 1432*eta + 657*ipow(Si,2)*(-1 + 3*eta)))) - 48*ipow(e,4)*(8712 + 6*x*(12095 - 4000*eta + 1033*ipow(Si,2)*(-1 + 3*eta)) + k*(17364 + x*(138668 - 47284*eta + 12475*ipow(Si,2)*(-1 + 3*eta))))))/(27648.*(1 + 2*k)*(3 + 2*k)*n);
    Dx(1,2) = (C2beta*Ci*e*(-9216*(3*(36 + x*(188 - 12*eta + 19*ipow(Si,2)*(-1 + 3*eta))) + 4*k*(36 + x*(130 - 12*eta + 23*ipow(Si,2)*(-1 + 3*eta)))) + 1152*ipow(e,2)*(8*k*(276 + x*(2583 - 434*eta + 180*ipow(Si,2)*(-1 + 3*eta))) + 3*(420 + x*(3744 - 680*eta + 261*ipow(Si,2)*(-1 + 3*eta)))) - 48*ipow(e,4)*(3*(5748 + x*(41908 - 15284*eta + 4211*ipow(Si,2)*(-1 + 3*eta))) + 4*k*(8652 + x*(66098 - 23284*eta + 6277*ipow(Si,2)*(-1 + 3*eta)))) + ipow(e,6)*(8*k*(59220 + x*(584411 - 222510*eta + 51394*ipow(Si,2)*(-1 + 3*eta))) + 3*(81684 + x*(894872 - 314400*eta + 69253*ipow(Si,2)*(-1 + 3*eta))))))/(55296.*(1 + 2*k)*(3 + 2*k)*n);
    Ap(1,3) = ((1 + ipow(Ci,2))*e*(-144*(4 + 15*k) + 216*ipow(e,2)*(36 + 65*k) - 9*ipow(e,4)*(1700 + 2683*k) + 4*ipow(e,6)*(2839 + 4308*k))*S3beta*Si*sqrt(x)*delta)/(144.*(2 + 3*k)*(4 + 3*k)*n);
    Bp(1,3) = -(C3beta*(1 + ipow(Ci,2))*e*(-144*(28 + 33*k) + 72*ipow(e,2)*(196 + 261*k) + 4*ipow(e,6)*(2971 + 4407*k) - 3*ipow(e,4)*(6164 + 8847*k))*Si*sqrt(x)*delta)/(144.*(2 + 3*k)*(4 + 3*k)*n);
    Cp(1,3) = -(C3beta*(1 + ipow(Ci,2))*e*(-144*(4 + 15*k) + 216*ipow(e,2)*(36 + 65*k) - 9*ipow(e,4)*(1700 + 2683*k) + 4*ipow(e,6)*(2839 + 4308*k))*Si*sqrt(x)*delta)/(144.*(2 + 3*k)*(4 + 3*k)*n);
    Dp(1,3) = -((1 + ipow(Ci,2))*e*(-144*(28 + 33*k) + 72*ipow(e,2)*(196 + 261*k) + 4*ipow(e,6)*(2971 + 4407*k) - 3*ipow(e,4)*(6164 + 8847*k))*S3beta*Si*sqrt(x)*delta)/(144.*(2 + 3*k)*(4 + 3*k)*n);
    Ax(1,3) = (C3beta*Ci*e*(-144*(4 + 15*k) + 216*ipow(e,2)*(36 + 65*k) - 9*ipow(e,4)*(1700 + 2683*k) + 4*ipow(e,6)*(2839 + 4308*k))*Si*sqrt(x)*delta)/(72.*(2 + 3*k)*(4 + 3*k)*n);
    Bx(1,3) = (Ci*e*(-144*(28 + 33*k) + 72*ipow(e,2)*(196 + 261*k) + 4*ipow(e,6)*(2971 + 4407*k) - 3*ipow(e,4)*(6164 + 8847*k))*S3beta*Si*sqrt(x)*delta)/(72.*(2 + 3*k)*(4 + 3*k)*n);
    Cx(1,3) = (Ci*e*(-144*(4 + 15*k) + 216*ipow(e,2)*(36 + 65*k) - 9*ipow(e,4)*(1700 + 2683*k) + 4*ipow(e,6)*(2839 + 4308*k))*S3beta*Si*sqrt(x)*delta)/(72.*(2 + 3*k)*(4 + 3*k)*n);
    Dx(1,3) = (C3beta*Ci*e*(144*(28 + 33*k) - 72*ipow(e,2)*(196 + 261*k) - 4*ipow(e,6)*(2971 + 4407*k) + 3*ipow(e,4)*(6164 + 8847*k))*Si*sqrt(x)*delta)/(72.*(2 + 3*k)*(4 + 3*k)*n);
    Ap(1,4) = ((1 + ipow(Ci,2))*e*(-9216*(165 + 382*k) + 10368*ipow(e,2)*(2235 + 3538*k) - 48*ipow(e,4)*(1586445 + 2266406*k) + ipow(e,6)*(104934855 + 143766634*k))*S4beta*ipow(Si,2)*x*(-1 + 3*eta))/(110592.*(3 + 4*k)*(5 + 4*k)*n);
    Bp(1,4) = -(C4beta*(1 + ipow(Ci,2))*e*(-9216*(1545 + 1736*k) + 10368*ipow(e,2)*(8655 + 10424*k) - 48*ipow(e,4)*(4306485 + 5439688*k) + ipow(e,6)*(238770915 + 310654232*k))*ipow(Si,2)*x*(-1 + 3*eta))/(221184.*(3 + 4*k)*(5 + 4*k)*n);
    Cp(1,4) = -(C4beta*(1 + ipow(Ci,2))*e*(-9216*(165 + 382*k) + 10368*ipow(e,2)*(2235 + 3538*k) - 48*ipow(e,4)*(1586445 + 2266406*k) + ipow(e,6)*(104934855 + 143766634*k))*ipow(Si,2)*x*(-1 + 3*eta))/(110592.*(3 + 4*k)*(5 + 4*k)*n);
    Dp(1,4) = -((1 + ipow(Ci,2))*e*(-9216*(1545 + 1736*k) + 10368*ipow(e,2)*(8655 + 10424*k) - 48*ipow(e,4)*(4306485 + 5439688*k) + ipow(e,6)*(238770915 + 310654232*k))*S4beta*ipow(Si,2)*x*(-1 + 3*eta))/(221184.*(3 + 4*k)*(5 + 4*k)*n);
    Ax(1,4) = (C4beta*Ci*e*(-9216*(165 + 382*k) + 10368*ipow(e,2)*(2235 + 3538*k) - 48*ipow(e,4)*(1586445 + 2266406*k) + ipow(e,6)*(104934855 + 143766634*k))*ipow(Si,2)*x*(-1 + 3*eta))/(55296.*(3 + 4*k)*(5 + 4*k)*n);
    Bx(1,4) = (Ci*e*(-9216*(1545 + 1736*k) + 10368*ipow(e,2)*(8655 + 10424*k) - 48*ipow(e,4)*(4306485 + 5439688*k) + ipow(e,6)*(238770915 + 310654232*k))*S4beta*ipow(Si,2)*x*(-1 + 3*eta))/(110592.*(3 + 4*k)*(5 + 4*k)*n);
    Cx(1,4) = (Ci*e*(-9216*(165 + 382*k) + 10368*ipow(e,2)*(2235 + 3538*k) - 48*ipow(e,4)*(1586445 + 2266406*k) + ipow(e,6)*(104934855 + 143766634*k))*S4beta*ipow(Si,2)*x*(-1 + 3*eta))/(55296.*(3 + 4*k)*(5 + 4*k)*n);
    Dx(1,4) = -(C4beta*Ci*e*(-9216*(1545 + 1736*k) + 10368*ipow(e,2)*(8655 + 10424*k) - 48*ipow(e,4)*(4306485 + 5439688*k) + ipow(e,6)*(238770915 + 310654232*k))*ipow(Si,2)*x*(-1 + 3*eta))/(110592.*(3 + 4*k)*(5 + 4*k)*n);
    Ap(2,0) = 0;
    Bp(2,0) = (ipow(e,2)*(-2*(-1 + ipow(Ci,2))*x*(600*ipow(e,2)*(9 + 7*eta) - 360*(45 + 11*eta) - 15*ipow(e,4)*(309 + 59*eta) + ipow(e,6)*(-2907 + 83*eta)) - 3*ipow(Si,2)*(-720*(2 + (1 + ipow(Ci,2))*x*(-1 + 3*eta)) - 60*ipow(e,4)*(1 + 3*(1 + ipow(Ci,2))*x*(-1 + 3*eta)) + 120*ipow(e,2)*(4 + 7*(1 + ipow(Ci,2))*x*(-1 + 3*eta)) + ipow(e,6)*(4 + 17*(1 + ipow(Ci,2))*x*(-1 + 3*eta)))))/(8640.*n);
    Cp(2,0) = 0;
    Dp(2,0) = 0;
    Ax(2,0) = -(Ci*ipow(e,2)*(-180 + 150*ipow(e,2) - 15*ipow(e,4) + 8*ipow(e,6))*ipow(Si,2)*x*(-1 + 3*eta))/(360.*n);
    Bx(2,0) = 0;
    Cx(2,0) = 0;
    Dx(2,0) = 0;
    Ap(2,1) = (ipow(e,2)*(ipow(e,6)*(71595 - 74658*k) + 103680*(-4 + 5*k) - 5760*ipow(e,2)*(-132 + 133*k) + 30*ipow(e,4)*(-11631 + 11522*k) + ipow(Ci,2)*(-11520*(-12 + 23*k) + 1920*ipow(e,2)*(-228 + 221*k) - 30*ipow(e,4)*(-6909 + 7174*k) + ipow(e,6)*(-49137 + 44614*k)))*S1beta*Si*sqrt(x)*delta)/(92160.*(-1 + k)*(3 + k)*n);
    Bp(2,1) = (C1beta*ipow(e,2)*(ipow(Ci,2)*(ipow(e,6)*(84705 - 93751*k) + 11520*(-57 + 35*k) - 1920*ipow(e,2)*(-435 + 449*k) + 30*ipow(e,4)*(-14613 + 14083*k)) + 3*(ipow(e,4)*(229350 - 231530*k) - 34560*(-11 + 9*k) + 1920*ipow(e,2)*(-267 + 265*k) + ipow(e,6)*(-50793 + 48751*k)))*Si*sqrt(x)*delta)/(184320.*(-1 + k)*(3 + k)*n);
    Cp(2,1) = (C1beta*ipow(e,2)*(ipow(e,4)*(348930 - 345660*k) - 103680*(-4 + 5*k) + 5760*ipow(e,2)*(-132 + 133*k) + ipow(e,6)*(-71595 + 74658*k) + ipow(Ci,2)*(ipow(e,6)*(49137 - 44614*k) + 11520*(-12 + 23*k) - 1920*ipow(e,2)*(-228 + 221*k) + 30*ipow(e,4)*(-6909 + 7174*k)))*Si*sqrt(x)*delta)/(92160.*(-1 + k)*(3 + k)*n);
    Dp(2,1) = (ipow(e,2)*(ipow(Ci,2)*(ipow(e,6)*(84705 - 93751*k) + 11520*(-57 + 35*k) - 1920*ipow(e,2)*(-435 + 449*k) + 30*ipow(e,4)*(-14613 + 14083*k)) + 3*(ipow(e,4)*(229350 - 231530*k) - 34560*(-11 + 9*k) + 1920*ipow(e,2)*(-267 + 265*k) + ipow(e,6)*(-50793 + 48751*k)))*S1beta*Si*sqrt(x)*delta)/(184320.*(-1 + k)*(3 + k)*n);
    Ax(2,1) = (C1beta*Ci*ipow(e,2)*(ipow(e,6)*(11229 - 15022*k) + 11520*(-12 + 11*k) - 1920*ipow(e,2)*(-84 + 89*k) + 30*ipow(e,4)*(-2361 + 2174*k))*Si*sqrt(x)*delta)/(46080.*(-1 + k)*(3 + k)*n);
    Bx(2,1) = (Ci*ipow(e,2)*(ipow(e,6)*(33837 - 26251*k) + 11520*(-21 + 23*k) - 1920*ipow(e,2)*(-183 + 173*k) + 30*ipow(e,4)*(-4161 + 4535*k))*S1beta*Si*sqrt(x)*delta)/(92160.*(-1 + k)*(3 + k)*n);
    Cx(2,1) = (Ci*ipow(e,2)*(ipow(e,6)*(11229 - 15022*k) + 11520*(-12 + 11*k) - 1920*ipow(e,2)*(-84 + 89*k) + 30*ipow(e,4)*(-2361 + 2174*k))*S1beta*Si*sqrt(x)*delta)/(46080.*(-1 + k)*(3 + k)*n);
    Dx(2,1) = (C1beta*Ci*ipow(e,2)*(-11520*(-21 + 23*k) + 1920*ipow(e,2)*(-183 + 173*k) - 30*ipow(e,4)*(-4161 + 4535*k) + ipow(e,6)*(-33837 + 26251*k))*Si*sqrt(x)*delta)/(92160.*(-1 + k)*(3 + k)*n);
    Ap(2,2) = (ipow(e,2)*S2beta*(15*ipow(e,4)*(1212 + x*(6141 + 12322*eta - 1668*ipow(Si,2)*(-1 + 3*eta))) - 720*(-12 + x*(-31 - 74*eta + 4*ipow(Si,2)*(-1 + 3*eta))) + 240*ipow(e,2)*(-90 + x*(-647 - 687*eta + 74*ipow(Si,2)*(-1 + 3*eta))) + ipow(e,6)*(-7062 + x*(-26045 - 96721*eta + 14174*ipow(Si,2)*(-1 + 3*eta))) + ipow(Ci,2)*(240*ipow(e,2)*(-90 + 74*ipow(Si,2)*x*(-1 + 3*eta) + 3*x*(-409 + 351*eta)) - 720*(-12 + x*(-87 + 94*eta + 4*ipow(Si,2)*(-1 + 3*eta))) - 15*ipow(e,4)*(-1212 + x*(-16917 + 20006*eta + 1668*ipow(Si,2)*(-1 + 3*eta))) + ipow(e,6)*(-7062 + x*(-110361 + 156227*eta + 14174*ipow(Si,2)*(-1 + 3*eta))))))/(2160.*(2 + k)*n);
    Bp(2,2) = (C2beta*ipow(e,2)*(ipow(e,6)*(7062 + x*(26045 + 96721*eta - 14174*ipow(Si,2)*(-1 + 3*eta))) + 720*(-12 + x*(-31 - 74*eta + 4*ipow(Si,2)*(-1 + 3*eta))) - 240*ipow(e,2)*(-90 + x*(-647 - 687*eta + 74*ipow(Si,2)*(-1 + 3*eta))) + 15*ipow(e,4)*(-1212 + x*(-6141 - 12322*eta + 1668*ipow(Si,2)*(-1 + 3*eta))) + ipow(Ci,2)*(-240*ipow(e,2)*(-90 + 74*ipow(Si,2)*x*(-1 + 3*eta) + 3*x*(-409 + 351*eta)) + ipow(e,6)*(7062 + x*(110361 - 156227*eta - 14174*ipow(Si,2)*(-1 + 3*eta))) + 720*(-12 + x*(-87 + 94*eta + 4*ipow(Si,2)*(-1 + 3*eta))) + 15*ipow(e,4)*(-1212 + x*(-16917 + 20006*eta + 1668*ipow(Si,2)*(-1 + 3*eta))))))/(2160.*(2 + k)*n);
    Cp(2,2) = (C2beta*ipow(e,2)*(ipow(e,6)*(7062 + x*(26045 + 96721*eta - 14174*ipow(Si,2)*(-1 + 3*eta))) + 720*(-12 + x*(-31 - 74*eta + 4*ipow(Si,2)*(-1 + 3*eta))) - 240*ipow(e,2)*(-90 + x*(-647 - 687*eta + 74*ipow(Si,2)*(-1 + 3*eta))) + 15*ipow(e,4)*(-1212 + x*(-6141 - 12322*eta + 1668*ipow(Si,2)*(-1 + 3*eta))) + ipow(Ci,2)*(-240*ipow(e,2)*(-90 + 74*ipow(Si,2)*x*(-1 + 3*eta) + 3*x*(-409 + 351*eta)) + ipow(e,6)*(7062 + x*(110361 - 156227*eta - 14174*ipow(Si,2)*(-1 + 3*eta))) + 720*(-12 + x*(-87 + 94*eta + 4*ipow(Si,2)*(-1 + 3*eta))) + 15*ipow(e,4)*(-1212 + x*(-16917 + 20006*eta + 1668*ipow(Si,2)*(-1 + 3*eta))))))/(2160.*(2 + k)*n);
    Dp(2,2) = (ipow(e,2)*S2beta*(ipow(e,6)*(7062 + x*(26045 + 96721*eta - 14174*ipow(Si,2)*(-1 + 3*eta))) + 720*(-12 + x*(-31 - 74*eta + 4*ipow(Si,2)*(-1 + 3*eta))) - 240*ipow(e,2)*(-90 + x*(-647 - 687*eta + 74*ipow(Si,2)*(-1 + 3*eta))) + 15*ipow(e,4)*(-1212 + x*(-6141 - 12322*eta + 1668*ipow(Si,2)*(-1 + 3*eta))) + ipow(Ci,2)*(-240*ipow(e,2)*(-90 + 74*ipow(Si,2)*x*(-1 + 3*eta) + 3*x*(-409 + 351*eta)) + ipow(e,6)*(7062 + x*(110361 - 156227*eta - 14174*ipow(Si,2)*(-1 + 3*eta))) + 720*(-12 + x*(-87 + 94*eta + 4*ipow(Si,2)*(-1 + 3*eta))) + 15*ipow(e,4)*(-1212 + x*(-16917 + 20006*eta + 1668*ipow(Si,2)*(-1 + 3*eta))))))/(2160.*(2 + k)*n);
    Ax(2,2) = (C2beta*Ci*ipow(e,2)*(ipow(e,6)*(-7062 + x*(-68203 + 29753*eta - 6905*ipow(Si,2)*(-1 + 3*eta))) + 720*(12 + x*(59 - 10*eta + 10*ipow(Si,2)*(-1 + 3*eta))) - 240*ipow(e,2)*(90 + x*(937 - 183*eta + 71*ipow(Si,2)*(-1 + 3*eta))) + 15*ipow(e,4)*(1212 + x*(11529 - 3842*eta + 1026*ipow(Si,2)*(-1 + 3*eta)))))/(1080.*(2 + k)*n);
    Bx(2,2) = (Ci*ipow(e,2)*S2beta*(ipow(e,6)*(-7062 + x*(-68203 + 29753*eta - 6905*ipow(Si,2)*(-1 + 3*eta))) + 720*(12 + x*(59 - 10*eta + 10*ipow(Si,2)*(-1 + 3*eta))) - 240*ipow(e,2)*(90 + x*(937 - 183*eta + 71*ipow(Si,2)*(-1 + 3*eta))) + 15*ipow(e,4)*(1212 + x*(11529 - 3842*eta + 1026*ipow(Si,2)*(-1 + 3*eta)))))/(1080.*(2 + k)*n);
    Cx(2,2) = (Ci*ipow(e,2)*S2beta*(ipow(e,6)*(-7062 + x*(-68203 + 29753*eta - 6905*ipow(Si,2)*(-1 + 3*eta))) + 720*(12 + x*(59 - 10*eta + 10*ipow(Si,2)*(-1 + 3*eta))) - 240*ipow(e,2)*(90 + x*(937 - 183*eta + 71*ipow(Si,2)*(-1 + 3*eta))) + 15*ipow(e,4)*(1212 + x*(11529 - 3842*eta + 1026*ipow(Si,2)*(-1 + 3*eta)))))/(1080.*(2 + k)*n);
    Dx(2,2) = (C2beta*Ci*ipow(e,2)*(-720*(12 + x*(59 - 10*eta + 10*ipow(Si,2)*(-1 + 3*eta))) + 240*ipow(e,2)*(90 + x*(937 - 183*eta + 71*ipow(Si,2)*(-1 + 3*eta))) - 15*ipow(e,4)*(1212 + x*(11529 - 3842*eta + 1026*ipow(Si,2)*(-1 + 3*eta))) + ipow(e,6)*(7062 + x*(68203 - 29753*eta + 6905*ipow(Si,2)*(-1 + 3*eta)))))/(1080.*(2 + k)*n);
    Ap(2,3) = ((1 + ipow(Ci,2))*ipow(e,2)*(-11520*(60 + 161*k) + 3200*ipow(e,2)*(860 + 2541*k) - 30*ipow(e,4)*(145565 + 436714*k) + ipow(e,6)*(3432145 + 10291162*k))*S3beta*Si*sqrt(x)*delta)/(30720.*(1 + 3*k)*(5 + 3*k)*n);
    Bp(2,3) = -(C3beta*(1 + ipow(Ci,2))*ipow(e,2)*(-11520*(265 + 909*k) + 9600*ipow(e,2)*(1655 + 5043*k) - 30*ipow(e,4)*(873485 + 2620341*k) + ipow(e,6)*(20566505 + 61731153*k))*Si*sqrt(x)*delta)/(184320.*(1 + 3*k)*(5 + 3*k)*n);
    Cp(2,3) = -(C3beta*(1 + ipow(Ci,2))*ipow(e,2)*(-11520*(60 + 161*k) + 3200*ipow(e,2)*(860 + 2541*k) - 30*ipow(e,4)*(145565 + 436714*k) + ipow(e,6)*(3432145 + 10291162*k))*Si*sqrt(x)*delta)/(30720.*(1 + 3*k)*(5 + 3*k)*n);
    Dp(2,3) = -((1 + ipow(Ci,2))*ipow(e,2)*(-11520*(265 + 909*k) + 9600*ipow(e,2)*(1655 + 5043*k) - 30*ipow(e,4)*(873485 + 2620341*k) + ipow(e,6)*(20566505 + 61731153*k))*S3beta*Si*sqrt(x)*delta)/(184320.*(1 + 3*k)*(5 + 3*k)*n);
    Ax(2,3) = (C3beta*Ci*ipow(e,2)*(-11520*(60 + 161*k) + 3200*ipow(e,2)*(860 + 2541*k) - 30*ipow(e,4)*(145565 + 436714*k) + ipow(e,6)*(3432145 + 10291162*k))*Si*sqrt(x)*delta)/(15360.*(1 + 3*k)*(5 + 3*k)*n);
    Bx(2,3) = (Ci*ipow(e,2)*(-11520*(265 + 909*k) + 9600*ipow(e,2)*(1655 + 5043*k) - 30*ipow(e,4)*(873485 + 2620341*k) + ipow(e,6)*(20566505 + 61731153*k))*S3beta*Si*sqrt(x)*delta)/(92160.*(1 + 3*k)*(5 + 3*k)*n);
    Cx(2,3) = (Ci*ipow(e,2)*(-11520*(60 + 161*k) + 3200*ipow(e,2)*(860 + 2541*k) - 30*ipow(e,4)*(145565 + 436714*k) + ipow(e,6)*(3432145 + 10291162*k))*S3beta*Si*sqrt(x)*delta)/(15360.*(1 + 3*k)*(5 + 3*k)*n);
    Dx(2,3) = (C3beta*Ci*ipow(e,2)*(11520*(265 + 909*k) - 9600*ipow(e,2)*(1655 + 5043*k) + 30*ipow(e,4)*(873485 + 2620341*k) - ipow(e,6)*(20566505 + 61731153*k))*Si*sqrt(x)*delta)/(92160.*(1 + 3*k)*(5 + 3*k)*n);
    Ap(2,4) = ((1 + ipow(Ci,2))*ipow(e,2)*(-1440*(285 + 514*k) - 810*ipow(e,4)*(7887 + 15698*k) + 240*ipow(e,2)*(10617 + 20686*k) + ipow(e,6)*(8314053 + 16617670*k))*S4beta*ipow(Si,2)*x*(-1 + 3*eta))/(17280.*(1 + 2*k)*(3 + 2*k)*n);
    Bp(2,4) = -(C4beta*(1 + ipow(Ci,2))*ipow(e,2)*(-1440*(201 + 458*k) - 810*ipow(e,4)*(7773 + 15622*k) + 240*ipow(e,2)*(9795 + 20138*k) + ipow(e,6)*(8298399 + 16607234*k))*ipow(Si,2)*x*(-1 + 3*eta))/(17280.*(1 + 2*k)*(3 + 2*k)*n);
    Cp(2,4) = -(C4beta*(1 + ipow(Ci,2))*ipow(e,2)*(-1440*(285 + 514*k) - 810*ipow(e,4)*(7887 + 15698*k) + 240*ipow(e,2)*(10617 + 20686*k) + ipow(e,6)*(8314053 + 16617670*k))*ipow(Si,2)*x*(-1 + 3*eta))/(17280.*(1 + 2*k)*(3 + 2*k)*n);
    Dp(2,4) = -((1 + ipow(Ci,2))*ipow(e,2)*(-1440*(201 + 458*k) - 810*ipow(e,4)*(7773 + 15622*k) + 240*ipow(e,2)*(9795 + 20138*k) + ipow(e,6)*(8298399 + 16607234*k))*S4beta*ipow(Si,2)*x*(-1 + 3*eta))/(17280.*(1 + 2*k)*(3 + 2*k)*n);
    Ax(2,4) = (C4beta*Ci*ipow(e,2)*(-1440*(285 + 514*k) - 810*ipow(e,4)*(7887 + 15698*k) + 240*ipow(e,2)*(10617 + 20686*k) + ipow(e,6)*(8314053 + 16617670*k))*ipow(Si,2)*x*(-1 + 3*eta))/(8640.*(1 + 2*k)*(3 + 2*k)*n);
    Bx(2,4) = (Ci*ipow(e,2)*(-1440*(201 + 458*k) - 810*ipow(e,4)*(7773 + 15622*k) + 240*ipow(e,2)*(9795 + 20138*k) + ipow(e,6)*(8298399 + 16607234*k))*S4beta*ipow(Si,2)*x*(-1 + 3*eta))/(8640.*(1 + 2*k)*(3 + 2*k)*n);
    Cx(2,4) = (Ci*ipow(e,2)*(-1440*(285 + 514*k) - 810*ipow(e,4)*(7887 + 15698*k) + 240*ipow(e,2)*(10617 + 20686*k) + ipow(e,6)*(8314053 + 16617670*k))*S4beta*ipow(Si,2)*x*(-1 + 3*eta))/(8640.*(1 + 2*k)*(3 + 2*k)*n);
    Dx(2,4) = -(C4beta*Ci*ipow(e,2)*(-1440*(201 + 458*k) - 810*ipow(e,4)*(7773 + 15622*k) + 240*ipow(e,2)*(9795 + 20138*k) + ipow(e,6)*(8298399 + 16607234*k))*ipow(Si,2)*x*(-1 + 3*eta))/(8640.*(1 + 2*k)*(3 + 2*k)*n);
    Ap(3,0) = 0;
    Bp(3,0) = (ipow(e,3)*(4*(-1 + ipow(Ci,2))*x*(640*(63 + 23*eta) - 120*ipow(e,2)*(233 + 141*eta) + 3*ipow(e,4)*(4765 + 1917*eta)) + 3*ipow(Si,2)*(640*(8 + 9*(1 + ipow(Ci,2))*x*(-1 + 3*eta)) - 360*ipow(e,2)*(8 + 19*(1 + ipow(Ci,2))*x*(-1 + 3*eta)) + 81*ipow(e,4)*(8 + 29*(1 + ipow(Ci,2))*x*(-1 + 3*eta)))))/(40960.*n);
    Cp(3,0) = 0;
    Dp(3,0) = 0;
    Ax(3,0) = (9*Ci*ipow(e,3)*(640 - 680*ipow(e,2) + 181*ipow(e,4))*ipow(Si,2)*x*(-1 + 3*eta))/(10240.*n);
    Bx(3,0) = 0;
    Cx(3,0) = 0;
    Dx(3,0) = 0;
    Ap(3,1) = (ipow(e,3)*(160*(-16 + 9*k) - 30*ipow(e,2)*(-164 + 83*k) + 3*ipow(e,4)*(-996 + 497*k) + ipow(Ci,2)*(ipow(e,4)*(1860 - 937*k) + 160*(8 - 5*k) + 10*ipow(e,2)*(-292 + 147*k)))*S1beta*Si*sqrt(x)*delta)/(160.*(-2 + k)*(4 + k)*n);
    Bp(3,1) = (C1beta*ipow(e,3)*(8960 + ipow(e,4)*(8940 - 4479*k) - 4000*k + 30*ipow(e,2)*(-500 + 247*k) + ipow(Ci,2)*(ipow(e,2)*(8840 - 4390*k) + 160*(-32 + 13*k) + ipow(e,4)*(-5636 + 2797*k)))*Si*sqrt(x)*delta)/(480.*(-2 + k)*(4 + k)*n);
    Cp(3,1) = (C1beta*ipow(e,3)*(160*(16 - 9*k) + 30*ipow(e,2)*(-164 + 83*k) - 3*ipow(e,4)*(-996 + 497*k) + ipow(Ci,2)*(ipow(e,2)*(2920 - 1470*k) + 160*(-8 + 5*k) + ipow(e,4)*(-1860 + 937*k)))*Si*sqrt(x)*delta)/(160.*(-2 + k)*(4 + k)*n);
    Dp(3,1) = (ipow(e,3)*(8960 + ipow(e,4)*(8940 - 4479*k) - 4000*k + 30*ipow(e,2)*(-500 + 247*k) + ipow(Ci,2)*(ipow(e,2)*(8840 - 4390*k) + 160*(-32 + 13*k) + ipow(e,4)*(-5636 + 2797*k)))*S1beta*Si*sqrt(x)*delta)/(480.*(-2 + k)*(4 + k)*n);
    Ax(3,1) = (C1beta*Ci*ipow(e,3)*(ipow(e,2)*(1000 - 510*k) + 320*(-2 + k) + ipow(e,4)*(-564 + 277*k))*Si*sqrt(x)*delta)/(80.*(-2 + k)*(4 + k)*n);
    Bx(3,1) = (Ci*ipow(e,3)*(ipow(e,2)*(3080 - 1510*k) + 960*(-2 + k) + ipow(e,4)*(-1652 + 841*k))*S1beta*Si*sqrt(x)*delta)/(240.*(-2 + k)*(4 + k)*n);
    Cx(3,1) = (Ci*ipow(e,3)*(ipow(e,2)*(1000 - 510*k) + 320*(-2 + k) + ipow(e,4)*(-564 + 277*k))*S1beta*Si*sqrt(x)*delta)/(80.*(-2 + k)*(4 + k)*n);
    Dx(3,1) = (C1beta*Ci*ipow(e,3)*(ipow(e,4)*(1652 - 841*k) - 960*(-2 + k) + 10*ipow(e,2)*(-308 + 151*k))*Si*sqrt(x)*delta)/(240.*(-2 + k)*(4 + k)*n);
    Ap(3,2) = -(S2beta*(24*(1 + ipow(Ci,2))*ipow(e,3)*(50560 - 134140*ipow(e,2) + 131777*ipow(e,4)) - x*(-7569920*ipow(e,3)*(-1 + 3*eta) + 26678080*ipow(e,5)*(-1 + 3*eta) - 34071124*ipow(e,7)*(-1 + 3*eta) + (1 + ipow(Ci,2))*(640*ipow(e,3)*(-19872 + 20306*eta + 1237*ipow(Si,2)*(-1 + 3*eta)) + ipow(e,7)*(-51668616 + 62634586*eta + 5480261*ipow(Si,2)*(-1 + 3*eta)) - 40*ipow(e,5)*(94427*ipow(Si,2)*(-1 + 3*eta) + 2*(-620847 + 601444*eta)))) + 2*(1 + k)*(-6*(1 + ipow(Ci,2))*ipow(e,3)*(65920 - 179480*ipow(e,2) + 175339*ipow(e,4)) + x*(-2543360*ipow(e,3)*(-1 + 3*eta) + 8890640*ipow(e,5)*(-1 + 3*eta) - 11357542*ipow(e,7)*(-1 + 3*eta) + (1 + ipow(Ci,2))*(640*ipow(e,3)*(-6351 + 6773*eta + 421*ipow(Si,2)*(-1 + 3*eta)) - 40*ipow(e,5)*(-417627 + 401329*eta + 31441*ipow(Si,2)*(-1 + 3*eta)) + ipow(e,7)*(-17085003 + 20864513*eta + 1827413*ipow(Si,2)*(-1 + 3*eta)))))))/(30720.*(-1 + 2*k)*(5 + 2*k)*n);
    Bp(3,2) = -(C2beta*(9*(-6*(1 + ipow(Ci,2))*ipow(e,3)*(65920 - 179480*ipow(e,2) + 175339*ipow(e,4)) + x*(-2543360*ipow(e,3)*(-1 + 3*eta) + 8890640*ipow(e,5)*(-1 + 3*eta) - 11357542*ipow(e,7)*(-1 + 3*eta) + (1 + ipow(Ci,2))*(640*ipow(e,3)*(-6351 + 6773*eta + 421*ipow(Si,2)*(-1 + 3*eta)) - 40*ipow(e,5)*(-417627 + 401329*eta + 31441*ipow(Si,2)*(-1 + 3*eta)) + ipow(e,7)*(-17085003 + 20864513*eta + 1827413*ipow(Si,2)*(-1 + 3*eta))))) + 2*(1 + k)*(24*(1 + ipow(Ci,2))*ipow(e,3)*(50560 - 134140*ipow(e,2) + 131777*ipow(e,4)) + x*(7569920*ipow(e,3)*(-1 + 3*eta) - 26678080*ipow(e,5)*(-1 + 3*eta) + 34071124*ipow(e,7)*(-1 + 3*eta) - (1 + ipow(Ci,2))*(640*ipow(e,3)*(-19872 + 20306*eta + 1237*ipow(Si,2)*(-1 + 3*eta)) + ipow(e,7)*(-51668616 + 62634586*eta + 5480261*ipow(Si,2)*(-1 + 3*eta)) - 40*ipow(e,5)*(94427*ipow(Si,2)*(-1 + 3*eta) + 2*(-620847 + 601444*eta)))))))/(92160.*(-1 + 2*k)*(5 + 2*k)*n);
    Cp(3,2) = -(C2beta*(-24*(1 + ipow(Ci,2))*ipow(e,3)*(50560 - 134140*ipow(e,2) + 131777*ipow(e,4)) + x*(-7569920*ipow(e,3)*(-1 + 3*eta) + 26678080*ipow(e,5)*(-1 + 3*eta) - 34071124*ipow(e,7)*(-1 + 3*eta) + (1 + ipow(Ci,2))*(640*ipow(e,3)*(-19872 + 20306*eta + 1237*ipow(Si,2)*(-1 + 3*eta)) + ipow(e,7)*(-51668616 + 62634586*eta + 5480261*ipow(Si,2)*(-1 + 3*eta)) - 40*ipow(e,5)*(94427*ipow(Si,2)*(-1 + 3*eta) + 2*(-620847 + 601444*eta)))) - 2*(1 + k)*(-6*(1 + ipow(Ci,2))*ipow(e,3)*(65920 - 179480*ipow(e,2) + 175339*ipow(e,4)) + x*(-2543360*ipow(e,3)*(-1 + 3*eta) + 8890640*ipow(e,5)*(-1 + 3*eta) - 11357542*ipow(e,7)*(-1 + 3*eta) + (1 + ipow(Ci,2))*(640*ipow(e,3)*(-6351 + 6773*eta + 421*ipow(Si,2)*(-1 + 3*eta)) - 40*ipow(e,5)*(-417627 + 401329*eta + 31441*ipow(Si,2)*(-1 + 3*eta)) + ipow(e,7)*(-17085003 + 20864513*eta + 1827413*ipow(Si,2)*(-1 + 3*eta)))))))/(30720.*(-1 + 2*k)*(5 + 2*k)*n);
    Dp(3,2) = -(S2beta*(9*(-6*(1 + ipow(Ci,2))*ipow(e,3)*(65920 - 179480*ipow(e,2) + 175339*ipow(e,4)) + x*(-2543360*ipow(e,3)*(-1 + 3*eta) + 8890640*ipow(e,5)*(-1 + 3*eta) - 11357542*ipow(e,7)*(-1 + 3*eta) + (1 + ipow(Ci,2))*(640*ipow(e,3)*(-6351 + 6773*eta + 421*ipow(Si,2)*(-1 + 3*eta)) - 40*ipow(e,5)*(-417627 + 401329*eta + 31441*ipow(Si,2)*(-1 + 3*eta)) + ipow(e,7)*(-17085003 + 20864513*eta + 1827413*ipow(Si,2)*(-1 + 3*eta))))) - 2*(1 + k)*(-24*(1 + ipow(Ci,2))*ipow(e,3)*(50560 - 134140*ipow(e,2) + 131777*ipow(e,4)) + x*(-7569920*ipow(e,3)*(-1 + 3*eta) + 26678080*ipow(e,5)*(-1 + 3*eta) - 34071124*ipow(e,7)*(-1 + 3*eta) + (1 + ipow(Ci,2))*(640*ipow(e,3)*(-19872 + 20306*eta + 1237*ipow(Si,2)*(-1 + 3*eta)) + ipow(e,7)*(-51668616 + 62634586*eta + 5480261*ipow(Si,2)*(-1 + 3*eta)) - 40*ipow(e,5)*(94427*ipow(Si,2)*(-1 + 3*eta) + 2*(-620847 + 601444*eta)))))))/(92160.*(-1 + 2*k)*(5 + 2*k)*n);
    Ax(3,2) = (C2beta*Ci*ipow(e,3)*(640*(-660 + x*(-5230 + 940*eta - 575*ipow(Si,2)*(-1 + 3*eta)) + k*(1236 + x*(8728 - 1624*eta + 1145*ipow(Si,2)*(-1 + 3*eta)))) - 40*ipow(e,2)*(-30*(888 + x*(9841 - 2220*eta + 802*ipow(Si,2)*(-1 + 3*eta))) + k*(53844 + x*(612988 - 135860*eta + 48251*ipow(Si,2)*(-1 + 3*eta)))) + ipow(e,4)*(-5*(211716 + x*(2364118 - 774300*eta + 202715*ipow(Si,2)*(-1 + 3*eta))) + k*(2104068 + x*(22812464 - 7656400*eta + 2023945*ipow(Si,2)*(-1 + 3*eta))))))/(15360.*(-1 + 2*k)*(5 + 2*k)*n);
    Bx(3,2) = (Ci*ipow(e,3)*S2beta*(640*(-5*(708 + x*(4544 - 872*eta + 685*ipow(Si,2)*(-1 + 3*eta))) + 8*k*(948 + x*(6979 - 1282*eta + 860*ipow(Si,2)*(-1 + 3*eta)))) - 40*ipow(e,2)*(-5*(32532 + x*(376804 - 82580*eta + 29003*ipow(Si,2)*(-1 + 3*eta))) + 4*k*(80484 + x*(908218 - 202460*eta + 72311*ipow(Si,2)*(-1 + 3*eta)))) + ipow(e,4)*(-5*(1257204 + x*(13355992 - 4559200*eta + 1213085*ipow(Si,2)*(-1 + 3*eta))) + 8*k*(1581324 + x*(17316527 - 5763950*eta + 1518760*ipow(Si,2)*(-1 + 3*eta))))))/(92160.*(-1 + 2*k)*(5 + 2*k)*n);
    Cx(3,2) = (Ci*ipow(e,3)*S2beta*(640*(-660 + x*(-5230 + 940*eta - 575*ipow(Si,2)*(-1 + 3*eta)) + k*(1236 + x*(8728 - 1624*eta + 1145*ipow(Si,2)*(-1 + 3*eta)))) - 40*ipow(e,2)*(-30*(888 + x*(9841 - 2220*eta + 802*ipow(Si,2)*(-1 + 3*eta))) + k*(53844 + x*(612988 - 135860*eta + 48251*ipow(Si,2)*(-1 + 3*eta)))) + ipow(e,4)*(-5*(211716 + x*(2364118 - 774300*eta + 202715*ipow(Si,2)*(-1 + 3*eta))) + k*(2104068 + x*(22812464 - 7656400*eta + 2023945*ipow(Si,2)*(-1 + 3*eta))))))/(15360.*(-1 + 2*k)*(5 + 2*k)*n);
    Dx(3,2) = (C2beta*Ci*ipow(e,3)*(-640*(-5*(708 + x*(4544 - 872*eta + 685*ipow(Si,2)*(-1 + 3*eta))) + 8*k*(948 + x*(6979 - 1282*eta + 860*ipow(Si,2)*(-1 + 3*eta)))) + 40*ipow(e,2)*(-5*(32532 + x*(376804 - 82580*eta + 29003*ipow(Si,2)*(-1 + 3*eta))) + 4*k*(80484 + x*(908218 - 202460*eta + 72311*ipow(Si,2)*(-1 + 3*eta)))) + ipow(e,4)*(6286020 + 5*x*(13355992 - 4559200*eta + 1213085*ipow(Si,2)*(-1 + 3*eta)) - 8*k*(1581324 + x*(17316527 - 5763950*eta + 1518760*ipow(Si,2)*(-1 + 3*eta))))))/(92160.*(-1 + 2*k)*(5 + 2*k)*n);
    Ap(3,3) = (-27*(1 + ipow(Ci,2))*ipow(e,3)*(40 - 180*ipow(e,2) + 309*ipow(e,4))*S3beta*Si*sqrt(x)*delta)/(80.*(2 + k)*n);
    Bp(3,3) = (27*C3beta*(1 + ipow(Ci,2))*ipow(e,3)*(40 - 180*ipow(e,2) + 309*ipow(e,4))*Si*sqrt(x)*delta)/(80.*(2 + k)*n);
    Cp(3,3) = (27*C3beta*(1 + ipow(Ci,2))*ipow(e,3)*(40 - 180*ipow(e,2) + 309*ipow(e,4))*Si*sqrt(x)*delta)/(80.*(2 + k)*n);
    Dp(3,3) = (27*(1 + ipow(Ci,2))*ipow(e,3)*(40 - 180*ipow(e,2) + 309*ipow(e,4))*S3beta*Si*sqrt(x)*delta)/(80.*(2 + k)*n);
    Ax(3,3) = (-27*C3beta*Ci*ipow(e,3)*(40 - 180*ipow(e,2) + 309*ipow(e,4))*Si*sqrt(x)*delta)/(40.*(2 + k)*n);
    Bx(3,3) = (-27*Ci*ipow(e,3)*(40 - 180*ipow(e,2) + 309*ipow(e,4))*S3beta*Si*sqrt(x)*delta)/(40.*(2 + k)*n);
    Cx(3,3) = (-27*Ci*ipow(e,3)*(40 - 180*ipow(e,2) + 309*ipow(e,4))*S3beta*Si*sqrt(x)*delta)/(40.*(2 + k)*n);
    Dx(3,3) = (27*C3beta*Ci*ipow(e,3)*(40 - 180*ipow(e,2) + 309*ipow(e,4))*Si*sqrt(x)*delta)/(40.*(2 + k)*n);
    Ap(3,4) = -((1 + ipow(Ci,2))*ipow(e,3)*(640*(9695 + 39154*k) - 40*ipow(e,2)*(1048005 + 4195558*k) + ipow(e,4)*(109697707 + 438791450*k))*S4beta*ipow(Si,2)*x*(-1 + 3*eta))/(61440.*(1 + 4*k)*(7 + 4*k)*n);
    Bp(3,4) = (C4beta*(1 + ipow(Ci,2))*ipow(e,3)*(640*(59479 + 235672*k) - 40*ipow(e,2)*(6300413 + 25180424*k) + ipow(e,4)*(658188419 + 2632749944*k))*ipow(Si,2)*x*(-1 + 3*eta))/(368640.*(1 + 4*k)*(7 + 4*k)*n);
    Cp(3,4) = (C4beta*(1 + ipow(Ci,2))*ipow(e,3)*(640*(9695 + 39154*k) - 40*ipow(e,2)*(1048005 + 4195558*k) + ipow(e,4)*(109697707 + 438791450*k))*ipow(Si,2)*x*(-1 + 3*eta))/(61440.*(1 + 4*k)*(7 + 4*k)*n);
    Dp(3,4) = ((1 + ipow(Ci,2))*ipow(e,3)*(640*(59479 + 235672*k) - 40*ipow(e,2)*(6300413 + 25180424*k) + ipow(e,4)*(658188419 + 2632749944*k))*S4beta*ipow(Si,2)*x*(-1 + 3*eta))/(368640.*(1 + 4*k)*(7 + 4*k)*n);
    Ax(3,4) = -(C4beta*Ci*ipow(e,3)*(640*(9695 + 39154*k) - 40*ipow(e,2)*(1048005 + 4195558*k) + ipow(e,4)*(109697707 + 438791450*k))*ipow(Si,2)*x*(-1 + 3*eta))/(30720.*(1 + 4*k)*(7 + 4*k)*n);
    Bx(3,4) = -(Ci*ipow(e,3)*(640*(59479 + 235672*k) - 40*ipow(e,2)*(6300413 + 25180424*k) + ipow(e,4)*(658188419 + 2632749944*k))*S4beta*ipow(Si,2)*x*(-1 + 3*eta))/(184320.*(1 + 4*k)*(7 + 4*k)*n);
    Cx(3,4) = -(Ci*ipow(e,3)*(640*(9695 + 39154*k) - 40*ipow(e,2)*(1048005 + 4195558*k) + ipow(e,4)*(109697707 + 438791450*k))*S4beta*ipow(Si,2)*x*(-1 + 3*eta))/(30720.*(1 + 4*k)*(7 + 4*k)*n);
    Dx(3,4) = (C4beta*Ci*ipow(e,3)*(640*(59479 + 235672*k) - 40*ipow(e,2)*(6300413 + 25180424*k) + ipow(e,4)*(658188419 + 2632749944*k))*ipow(Si,2)*x*(-1 + 3*eta))/(184320.*(1 + 4*k)*(7 + 4*k)*n);
    Ap(4,0) = 0;
    Bp(4,0) = (ipow(e,4)*(2*(-1 + ipow(Ci,2))*x*(15*(81 + 35*eta) - 12*ipow(e,2)*(102 + 59*eta) + 4*ipow(e,4)*(162 + 83*eta)) + 3*ipow(Si,2)*(-48*ipow(e,2)*(1 + 3*(1 + ipow(Ci,2))*x*(-1 + 3*eta)) + 15*(4 + 7*(1 + ipow(Ci,2))*x*(-1 + 3*eta)) + 4*ipow(e,4)*(4 + 17*(1 + ipow(Ci,2))*x*(-1 + 3*eta)))))/(540.*n);
    Cp(4,0) = 0;
    Dp(4,0) = 0;
    Ax(4,0) = (Ci*ipow(e,4)*(120 - 156*ipow(e,2) + 65*ipow(e,4))*ipow(Si,2)*x*(-1 + 3*eta))/(180.*n);
    Bx(4,0) = 0;
    Cx(4,0) = 0;
    Dx(4,0) = 0;
    Ap(4,1) = (ipow(e,4)*(ipow(Ci,2)*(ipow(e,4)*(28490385 - 9544423*k) - 960*(-16395 + 6221*k) + 192*ipow(e,2)*(-190365 + 64427*k)) + 3*(960*(-9795 + 3541*k) - 64*ipow(e,2)*(-311295 + 105241*k) + ipow(e,4)*(-14951865 + 4992127*k)))*S1beta*Si*sqrt(x)*delta)/(737280.*(-3 + k)*(5 + k)*n);
    Bp(4,1) = (C1beta*ipow(e,4)*(ipow(e,4)*(22473765 - 7478997*k) - 2880*(-5415 + 1667*k) + 192*ipow(e,2)*(-158415 + 52067*k) + ipow(Ci,2)*(960*(-9615 + 2827*k) - 192*ipow(e,2)*(-97005 + 31849*k) + 7*ipow(e,4)*(-2047785 + 679193*k)))*Si*sqrt(x)*delta)/(368640.*(-3 + k)*(5 + k)*n);
    Cp(4,1) = (C1beta*ipow(e,4)*(ipow(Ci,2)*(960*(-16395 + 6221*k) - 192*ipow(e,2)*(-190365 + 64427*k) + 7*ipow(e,4)*(-4070055 + 1363489*k)) - 3*(960*(-9795 + 3541*k) - 64*ipow(e,2)*(-311295 + 105241*k) + ipow(e,4)*(-14951865 + 4992127*k)))*Si*sqrt(x)*delta)/(737280.*(-3 + k)*(5 + k)*n);
    Dp(4,1) = (ipow(e,4)*(ipow(e,4)*(22473765 - 7478997*k) - 2880*(-5415 + 1667*k) + 192*ipow(e,2)*(-158415 + 52067*k) + ipow(Ci,2)*(960*(-9615 + 2827*k) - 192*ipow(e,2)*(-97005 + 31849*k) + 7*ipow(e,4)*(-2047785 + 679193*k)))*S1beta*Si*sqrt(x)*delta)/(368640.*(-3 + k)*(5 + k)*n);
    Ax(4,1) = (C1beta*Ci*ipow(e,4)*(960*(-6495 + 2201*k) - 192*ipow(e,2)*(-60465 + 20407*k) + ipow(e,4)*(-8182605 + 2715979*k))*Si*sqrt(x)*delta)/(368640.*(-3 + k)*(5 + k)*n);
    Bx(4,1) = (Ci*ipow(e,4)*(960*(-3315 + 1087*k) - 192*ipow(e,2)*(-30705 + 10109*k) + ipow(e,4)*(-4069635 + 1362323*k))*S1beta*Si*sqrt(x)*delta)/(184320.*(-3 + k)*(5 + k)*n);
    Cx(4,1) = (Ci*ipow(e,4)*(960*(-6495 + 2201*k) - 192*ipow(e,2)*(-60465 + 20407*k) + ipow(e,4)*(-8182605 + 2715979*k))*S1beta*Si*sqrt(x)*delta)/(368640.*(-3 + k)*(5 + k)*n);
    Dx(4,1) = (C1beta*Ci*ipow(e,4)*(ipow(e,4)*(4069635 - 1362323*k) - 960*(-3315 + 1087*k) + 192*ipow(e,2)*(-30705 + 10109*k))*Si*sqrt(x)*delta)/(184320.*(-3 + k)*(5 + k)*n);
    Ap(4,2) = -(((1 + ipow(Ci,2))*ipow(e,4)*(117360 - 337992*ipow(e,2) + 377795*ipow(e,4))*S2beta)/1440. - 4*S2beta*x*((-2909*ipow(e,4)*(-1 + 3*eta))/18. + (8869*ipow(e,6)*(-1 + 3*eta))/15. - (715607*ipow(e,8)*(-1 + 3*eta))/864. + ((1 + ipow(Ci,2))*(240*ipow(e,4)*(-38721 + 40323*eta + 2906*ipow(Si,2)*(-1 + 3*eta)) - 24*ipow(e,6)*(-1464639 + 1530977*eta + 128304*ipow(Si,2)*(-1 + 3*eta)) + ipow(e,8)*(-41803791 + 52276225*eta + 4732786*ipow(Si,2)*(-1 + 3*eta))))/34560.) + 2*(1 + k)*((1 + ipow(Ci,2))*((-161*ipow(e,4))/8. + (2821*ipow(e,6))/48. - (377449*ipow(e,8))/5760.)*S2beta + (S2beta*x*(-5612160*ipow(e,4)*(-1 + 3*eta) + 20436480*ipow(e,6)*(-1 + 3*eta) - 28622632*ipow(e,8)*(-1 + 3*eta) + (1 + ipow(Ci,2))*(240*ipow(e,4)*(-37905 + 40353*eta + 2926*ipow(Si,2)*(-1 + 3*eta)) - 24*ipow(e,6)*(-1469667 + 1531795*eta + 128304*ipow(Si,2)*(-1 + 3*eta)) + ipow(e,8)*(-41695383 + 52259459*eta + 4732550*ipow(Si,2)*(-1 + 3*eta)))))/34560.))/(4.*(-1 + k)*(3 + k)*n);
    Bp(4,2) = -((C2beta*(-6*(1 + ipow(Ci,2))*ipow(e,4)*(115920 - 338520*ipow(e,2) + 377449*ipow(e,4)) + x*(-5612160*ipow(e,4)*(-1 + 3*eta) + 20436480*ipow(e,6)*(-1 + 3*eta) - 28622632*ipow(e,8)*(-1 + 3*eta) + (1 + ipow(Ci,2))*(240*ipow(e,4)*(-37905 + 40353*eta + 2926*ipow(Si,2)*(-1 + 3*eta)) - 24*ipow(e,6)*(-1469667 + 1531795*eta + 128304*ipow(Si,2)*(-1 + 3*eta)) + ipow(e,8)*(-41695383 + 52259459*eta + 4732550*ipow(Si,2)*(-1 + 3*eta))))))/8640. + 2*(1 + k)*((C2beta*(1 + ipow(Ci,2))*ipow(e,4)*(117360 - 337992*ipow(e,2) + 377795*ipow(e,4)))/5760. + C2beta*x*((2909*ipow(e,4)*(-1 + 3*eta))/18. - (8869*ipow(e,6)*(-1 + 3*eta))/15. + (715607*ipow(e,8)*(-1 + 3*eta))/864. - ((1 + ipow(Ci,2))*(240*ipow(e,4)*(-38721 + 40323*eta + 2906*ipow(Si,2)*(-1 + 3*eta)) - 24*ipow(e,6)*(-1464639 + 1530977*eta + 128304*ipow(Si,2)*(-1 + 3*eta)) + ipow(e,8)*(-41803791 + 52276225*eta + 4732786*ipow(Si,2)*(-1 + 3*eta))))/34560.)))/(4.*(-1 + k)*(3 + k)*n);
    Cp(4,2) = -(-2*(1 + k)*(C2beta*(1 + ipow(Ci,2))*((-161*ipow(e,4))/8. + (2821*ipow(e,6))/48. - (377449*ipow(e,8))/5760.) + (C2beta*x*(-5612160*ipow(e,4)*(-1 + 3*eta) + 20436480*ipow(e,6)*(-1 + 3*eta) - 28622632*ipow(e,8)*(-1 + 3*eta) + (1 + ipow(Ci,2))*(240*ipow(e,4)*(-37905 + 40353*eta + 2926*ipow(Si,2)*(-1 + 3*eta)) - 24*ipow(e,6)*(-1469667 + 1531795*eta + 128304*ipow(Si,2)*(-1 + 3*eta)) + ipow(e,8)*(-41695383 + 52259459*eta + 4732550*ipow(Si,2)*(-1 + 3*eta)))))/34560.) - 4*((C2beta*(1 + ipow(Ci,2))*ipow(e,4)*(117360 - 337992*ipow(e,2) + 377795*ipow(e,4)))/5760. + C2beta*x*((2909*ipow(e,4)*(-1 + 3*eta))/18. - (8869*ipow(e,6)*(-1 + 3*eta))/15. + (715607*ipow(e,8)*(-1 + 3*eta))/864. - ((1 + ipow(Ci,2))*(240*ipow(e,4)*(-38721 + 40323*eta + 2906*ipow(Si,2)*(-1 + 3*eta)) - 24*ipow(e,6)*(-1464639 + 1530977*eta + 128304*ipow(Si,2)*(-1 + 3*eta)) + ipow(e,8)*(-41803791 + 52276225*eta + 4732786*ipow(Si,2)*(-1 + 3*eta))))/34560.)))/(4.*(-1 + k)*(3 + k)*n);
    Dp(4,2) = -((S2beta*(-6*(1 + ipow(Ci,2))*ipow(e,4)*(115920 - 338520*ipow(e,2) + 377449*ipow(e,4)) + x*(-5612160*ipow(e,4)*(-1 + 3*eta) + 20436480*ipow(e,6)*(-1 + 3*eta) - 28622632*ipow(e,8)*(-1 + 3*eta) + (1 + ipow(Ci,2))*(240*ipow(e,4)*(-37905 + 40353*eta + 2926*ipow(Si,2)*(-1 + 3*eta)) - 24*ipow(e,6)*(-1469667 + 1531795*eta + 128304*ipow(Si,2)*(-1 + 3*eta)) + ipow(e,8)*(-41695383 + 52259459*eta + 4732550*ipow(Si,2)*(-1 + 3*eta))))))/8640. - 2*(1 + k)*(-((1 + ipow(Ci,2))*ipow(e,4)*(117360 - 337992*ipow(e,2) + 377795*ipow(e,4))*S2beta)/5760. + S2beta*x*((-2909*ipow(e,4)*(-1 + 3*eta))/18. + (8869*ipow(e,6)*(-1 + 3*eta))/15. - (715607*ipow(e,8)*(-1 + 3*eta))/864. + ((1 + ipow(Ci,2))*(240*ipow(e,4)*(-38721 + 40323*eta + 2906*ipow(Si,2)*(-1 + 3*eta)) - 24*ipow(e,6)*(-1464639 + 1530977*eta + 128304*ipow(Si,2)*(-1 + 3*eta)) + ipow(e,8)*(-41803791 + 52276225*eta + 4732786*ipow(Si,2)*(-1 + 3*eta))))/34560.)))/(4.*(-1 + k)*(3 + k)*n);
    Ax(4,2) = (C2beta*Ci*ipow(e,4)*(-24*ipow(e,2)*(-84366 + 84630*k - 3*x*(344649 - 84389*eta + 28176*ipow(Si,2)*(-1 + 3*eta)) + k*x*(1043907 - 254515*eta + 84576*ipow(Si,2)*(-1 + 3*eta))) + 240*(-3*(990 + x*(9319 - 1851*eta + 968*ipow(Si,2)*(-1 + 3*eta))) + k*(2898 + x*(26213 - 5277*eta + 2920*ipow(Si,2)*(-1 + 3*eta)))) + ipow(e,4)*(-3*(756282 + x*(9199745 - 3118033*eta + 807820*ipow(Si,2)*(-1 + 3*eta))) + k*(2264694 + x*(27384067 - 9325511*eta + 2423108*ipow(Si,2)*(-1 + 3*eta))))))/(34560.*(-1 + k)*(3 + k)*n);
    Bx(4,2) = (Ci*ipow(e,4)*S2beta*(-24*ipow(e,2)*(-84762 + 84498*k - 3*x*(349629 - 85063*eta + 28200*ipow(Si,2)*(-1 + 3*eta)) + k*x*(1038927 - 253841*eta + 84552*ipow(Si,2)*(-1 + 3*eta))) + 240*(-3*(954 + x*(8447 - 1713*eta + 976*ipow(Si,2)*(-1 + 3*eta))) + k*(2934 + x*(27085 - 5415*eta + 2912*ipow(Si,2)*(-1 + 3*eta)))) + ipow(e,4)*(-3*(754206 + x*(9092161 - 3103739*eta + 807644*ipow(Si,2)*(-1 + 3*eta))) + k*(2266770 + x*(27491651 - 9339805*eta + 2423284*ipow(Si,2)*(-1 + 3*eta))))))/(34560.*(-1 + k)*(3 + k)*n);
    Cx(4,2) = (Ci*ipow(e,4)*S2beta*(-24*ipow(e,2)*(-84366 + 84630*k - 3*x*(344649 - 84389*eta + 28176*ipow(Si,2)*(-1 + 3*eta)) + k*x*(1043907 - 254515*eta + 84576*ipow(Si,2)*(-1 + 3*eta))) + 240*(-3*(990 + x*(9319 - 1851*eta + 968*ipow(Si,2)*(-1 + 3*eta))) + k*(2898 + x*(26213 - 5277*eta + 2920*ipow(Si,2)*(-1 + 3*eta)))) + ipow(e,4)*(-3*(756282 + x*(9199745 - 3118033*eta + 807820*ipow(Si,2)*(-1 + 3*eta))) + k*(2264694 + x*(27384067 - 9325511*eta + 2423108*ipow(Si,2)*(-1 + 3*eta))))))/(34560.*(-1 + k)*(3 + k)*n);
    Dx(4,2) = (C2beta*Ci*ipow(e,4)*(24*ipow(e,2)*(-84762 + 84498*k - 3*x*(349629 - 85063*eta + 28200*ipow(Si,2)*(-1 + 3*eta)) + k*x*(1038927 - 253841*eta + 84552*ipow(Si,2)*(-1 + 3*eta))) + ipow(e,4)*(2262618 - 2266770*k + k*x*(-27491651 + 9339805*eta - 2423284*ipow(Si,2)*(-1 + 3*eta)) + 3*x*(9092161 - 3103739*eta + 807644*ipow(Si,2)*(-1 + 3*eta))) - 240*(-3*(954 + x*(8447 - 1713*eta + 976*ipow(Si,2)*(-1 + 3*eta))) + k*(2934 + x*(27085 - 5415*eta + 2912*ipow(Si,2)*(-1 + 3*eta))))))/(34560.*(-1 + k)*(3 + k)*n);
    Ap(4,3) = -((1 + ipow(Ci,2))*ipow(e,4)*(320*(-59087 + 176361*k) - 192*ipow(e,2)*(-450779 + 1353053*k) + ipow(e,4)*(-158280983 + 474762481*k))*S3beta*Si*sqrt(x)*delta)/(245760.*(-1 + 3*k)*(7 + 3*k)*n);
    Bp(4,3) = (C3beta*(1 + ipow(Ci,2))*ipow(e,4)*(960*(-29281 + 88293*k) - 192*ipow(e,2)*(-676795 + 2029311*k) + ipow(e,4)*(-237351065 + 712173897*k))*Si*sqrt(x)*delta)/(368640.*(-1 + 3*k)*(7 + 3*k)*n);
    Cp(4,3) = (C3beta*(1 + ipow(Ci,2))*ipow(e,4)*(320*(-59087 + 176361*k) - 192*ipow(e,2)*(-450779 + 1353053*k) + ipow(e,4)*(-158280983 + 474762481*k))*Si*sqrt(x)*delta)/(245760.*(-1 + 3*k)*(7 + 3*k)*n);
    Dp(4,3) = ((1 + ipow(Ci,2))*ipow(e,4)*(960*(-29281 + 88293*k) - 192*ipow(e,2)*(-676795 + 2029311*k) + ipow(e,4)*(-237351065 + 712173897*k))*S3beta*Si*sqrt(x)*delta)/(368640.*(-1 + 3*k)*(7 + 3*k)*n);
    Ax(4,3) = (C3beta*Ci*ipow(e,4)*(ipow(e,4)*(158280983 - 474762481*k) + 320*(59087 - 176361*k) + 192*ipow(e,2)*(-450779 + 1353053*k))*Si*sqrt(x)*delta)/(122880.*(-1 + 3*k)*(7 + 3*k)*n);
    Bx(4,3) = (Ci*ipow(e,4)*(ipow(e,4)*(237351065 - 712173897*k) - 960*(-29281 + 88293*k) + 192*ipow(e,2)*(-676795 + 2029311*k))*S3beta*Si*sqrt(x)*delta)/(184320.*(-1 + 3*k)*(7 + 3*k)*n);
    Cx(4,3) = (Ci*ipow(e,4)*(ipow(e,4)*(158280983 - 474762481*k) + 320*(59087 - 176361*k) + 192*ipow(e,2)*(-450779 + 1353053*k))*S3beta*Si*sqrt(x)*delta)/(122880.*(-1 + 3*k)*(7 + 3*k)*n);
    Dx(4,3) = (C3beta*Ci*ipow(e,4)*(960*(-29281 + 88293*k) - 192*ipow(e,2)*(-676795 + 2029311*k) + ipow(e,4)*(-237351065 + 712173897*k))*Si*sqrt(x)*delta)/(184320.*(-1 + 3*k)*(7 + 3*k)*n);
    Ap(4,4) = (-64*(1 + ipow(Ci,2))*ipow(e,4)*(120 - 792*ipow(e,2) + 2119*ipow(e,4))*S4beta*ipow(Si,2)*x*(-1 + 3*eta))/(135.*(2 + k)*n);
    Bp(4,4) = (64*C4beta*(1 + ipow(Ci,2))*ipow(e,4)*(120 - 792*ipow(e,2) + 2119*ipow(e,4))*ipow(Si,2)*x*(-1 + 3*eta))/(135.*(2 + k)*n);
    Cp(4,4) = (64*C4beta*(1 + ipow(Ci,2))*ipow(e,4)*(120 - 792*ipow(e,2) + 2119*ipow(e,4))*ipow(Si,2)*x*(-1 + 3*eta))/(135.*(2 + k)*n);
    Dp(4,4) = (64*(1 + ipow(Ci,2))*ipow(e,4)*(120 - 792*ipow(e,2) + 2119*ipow(e,4))*S4beta*ipow(Si,2)*x*(-1 + 3*eta))/(135.*(2 + k)*n);
    Ax(4,4) = (-128*C4beta*Ci*ipow(e,4)*(120 - 792*ipow(e,2) + 2119*ipow(e,4))*ipow(Si,2)*x*(-1 + 3*eta))/(135.*(2 + k)*n);
    Bx(4,4) = (-128*Ci*ipow(e,4)*(120 - 792*ipow(e,2) + 2119*ipow(e,4))*S4beta*ipow(Si,2)*x*(-1 + 3*eta))/(135.*(2 + k)*n);
    Cx(4,4) = (-128*Ci*ipow(e,4)*(120 - 792*ipow(e,2) + 2119*ipow(e,4))*S4beta*ipow(Si,2)*x*(-1 + 3*eta))/(135.*(2 + k)*n);
    Dx(4,4) = (128*C4beta*Ci*ipow(e,4)*(120 - 792*ipow(e,2) + 2119*ipow(e,4))*ipow(Si,2)*x*(-1 + 3*eta))/(135.*(2 + k)*n);
    Ap(5,0) = 0;
    Bp(5,0) = (-125*ipow(e,5)*(4*(-1 + ipow(Ci,2))*x*(-24*(99 + 47*eta) + ipow(e,2)*(3087 + 1775*eta)) + 3*ipow(Si,2)*(-24*(8 + 19*(1 + ipow(Ci,2))*x*(-1 + 3*eta)) + 25*ipow(e,2)*(8 + 29*(1 + ipow(Ci,2))*x*(-1 + 3*eta)))))/(221184.*n);
    Cp(5,0) = 0;
    Dp(5,0) = 0;
    Ax(5,0) = (-625*Ci*ipow(e,5)*(-24 + 37*ipow(e,2))*ipow(Si,2)*x*(-1 + 3*eta))/(18432.*n);
    Bx(5,0) = 0;
    Cx(5,0) = 0;
    Dx(5,0) = 0;
    Ap(5,1) = (ipow(e,5)*(-9*(2496 - 664*k + ipow(e,2)*(-5844 + 1483*k)) + ipow(Ci,2)*(13248 - 3624*k + ipow(e,2)*(-32892 + 8369*k)))*S1beta*Si*sqrt(x)*delta)/(288.*(-4 + k)*(6 + k)*n);
    Bp(5,1) = (C1beta*ipow(e,5)*(ipow(Ci,2)*(ipow(e,2)*(167964 - 41261*k) + 24*(-3072 + 703*k)) + 9*(40*(336 - 79*k) + ipow(e,2)*(-29748 + 7327*k)))*Si*sqrt(x)*delta)/(1440.*(-4 + k)*(6 + k)*n);
    Cp(5,1) = (C1beta*ipow(e,5)*(ipow(Ci,2)*(ipow(e,2)*(32892 - 8369*k) + 24*(-552 + 151*k)) + 9*(2496 - 664*k + ipow(e,2)*(-5844 + 1483*k)))*Si*sqrt(x)*delta)/(288.*(-4 + k)*(6 + k)*n);
    Dp(5,1) = (ipow(e,5)*(ipow(Ci,2)*(ipow(e,2)*(167964 - 41261*k) + 24*(-3072 + 703*k)) + 9*(40*(336 - 79*k) + ipow(e,2)*(-29748 + 7327*k)))*S1beta*Si*sqrt(x)*delta)/(1440.*(-4 + k)*(6 + k)*n);
    Ax(5,1) = (C1beta*Ci*ipow(e,5)*(ipow(e,2)*(9852 - 2489*k) + 24*(-192 + 49*k))*Si*sqrt(x)*delta)/(144.*(-4 + k)*(6 + k)*n);
    Bx(5,1) = (Ci*ipow(e,5)*(ipow(e,2)*(49884 - 12341*k) + 24*(-984 + 241*k))*S1beta*Si*sqrt(x)*delta)/(720.*(-4 + k)*(6 + k)*n);
    Cx(5,1) = (Ci*ipow(e,5)*(ipow(e,2)*(9852 - 2489*k) + 24*(-192 + 49*k))*S1beta*Si*sqrt(x)*delta)/(144.*(-4 + k)*(6 + k)*n);
    Dx(5,1) = (C1beta*Ci*ipow(e,5)*(23616 - 5784*k + ipow(e,2)*(-49884 + 12341*k))*Si*sqrt(x)*delta)/(720.*(-4 + k)*(6 + k)*n);
    Ap(5,2) = (ipow(e,5)*S2beta*(ipow(e,2)*(2*k*(-5295906 + 9859861*ipow(Si,2)*x*(-1 + 3*eta) - 25*x*(1567357 + 3022665*eta)) + 21*(755748 + x*(5565778 + 10794912*eta - 1408135*ipow(Si,2)*(-1 + 3*eta)))) - 24*(21*(10176 + x*(63212 + 121702*eta - 13397*ipow(Si,2)*(-1 + 3*eta))) + 2*k*(-70314 + x*(-413729 - 861805*eta + 94265*ipow(Si,2)*(-1 + 3*eta)))) + ipow(Ci,2)*(ipow(e,2)*(-21*(-755748 + 1408135*ipow(Si,2)*x*(-1 + 3*eta) + 6*x*(-2420287 + 2678820*eta)) + 2*k*(-5295906 + 9859861*ipow(Si,2)*x*(-1 + 3*eta) + 3*x*(-33965801 + 37524603*eta))) - 24*(21*(10176 + x*(160320 - 169622*eta - 13397*ipow(Si,2)*(-1 + 3*eta))) + 2*k*(-70314 + x*(-1097211 + 1188641*eta + 94265*ipow(Si,2)*(-1 + 3*eta)))))))/(55296.*(-3 + 2*k)*(7 + 2*k)*n);
    Bp(5,2) = -(((5*C2beta*(18*(1 + ipow(Ci,2))*ipow(e,5)*(-93752 + 294217*ipow(e,2)) + x*(-16403568*ipow(e,5)*(-1 + 3*eta) + 62713478*ipow(e,7)*(-1 + 3*eta) - (1 + ipow(Ci,2))*(-24*ipow(e,5)*(-1097211 + 1188641*eta + 94265*ipow(Si,2)*(-1 + 3*eta)) + ipow(e,7)*(9859861*ipow(Si,2)*(-1 + 3*eta) + 3*(-33965801 + 37524603*eta))))))/55296. + 2*(1 + k)*(C2beta*(1 + ipow(Ci,2))*((29527*ipow(e,5))/960. - (73507*ipow(e,7))/768.) + (C2beta*x*(81749568*ipow(e,5)*(-1 + 3*eta) - 313501780*ipow(e,7)*(-1 + 3*eta) + (1 + ipow(Ci,2))*(-24*ipow(e,5)*(-5561142 + 5939344*eta + 469867*ipow(Si,2)*(-1 + 3*eta)) + ipow(e,7)*(49290557*ipow(Si,2)*(-1 + 3*eta) + 6*(-84791828 + 93779823*eta)))))/276480.))/((-3 + 2*k)*(7 + 2*k)*n));
    Cp(5,2) = (C2beta*ipow(e,5)*(504*(10176 + x*(63212 + 121702*eta - 13397*ipow(Si,2)*(-1 + 3*eta))) + 48*k*(-70314 + x*(-413729 - 861805*eta + 94265*ipow(Si,2)*(-1 + 3*eta))) + ipow(e,2)*(21*(-755748 + 1408135*ipow(Si,2)*x*(-1 + 3*eta) - 2*x*(2782889 + 5397456*eta)) + 2*k*(5295906 + x*(39183925 + 75566625*eta - 9859861*ipow(Si,2)*(-1 + 3*eta)))) + ipow(Ci,2)*(504*(10176 + x*(160320 - 169622*eta - 13397*ipow(Si,2)*(-1 + 3*eta))) + 48*k*(-70314 + x*(-1097211 + 1188641*eta + 94265*ipow(Si,2)*(-1 + 3*eta))) + ipow(e,2)*(21*(-755748 + 1408135*ipow(Si,2)*x*(-1 + 3*eta) + 6*x*(-2420287 + 2678820*eta)) + k*(10591812 + x*(203794806 - 225147618*eta - 19719722*ipow(Si,2)*(-1 + 3*eta)))))))/(55296.*(-3 + 2*k)*(7 + 2*k)*n);
    Dp(5,2) = -(S2beta*(25*(18*(1 + ipow(Ci,2))*ipow(e,5)*(-93752 + 294217*ipow(e,2)) + x*(-16403568*ipow(e,5)*(-1 + 3*eta) + 62713478*ipow(e,7)*(-1 + 3*eta) - (1 + ipow(Ci,2))*(-24*ipow(e,5)*(-1097211 + 1188641*eta + 94265*ipow(Si,2)*(-1 + 3*eta)) + ipow(e,7)*(9859861*ipow(Si,2)*(-1 + 3*eta) + 3*(-33965801 + 37524603*eta))))) - 2*(1 + k)*(72*(1 + ipow(Ci,2))*ipow(e,5)*(-118108 + 367535*ipow(e,2)) + x*(-81749568*ipow(e,5)*(-1 + 3*eta) + 313501780*ipow(e,7)*(-1 + 3*eta) - (1 + ipow(Ci,2))*(-24*ipow(e,5)*(-5561142 + 5939344*eta + 469867*ipow(Si,2)*(-1 + 3*eta)) + ipow(e,7)*(49290557*ipow(Si,2)*(-1 + 3*eta) + 6*(-84791828 + 93779823*eta)))))))/(276480.*(-3 + 2*k)*(7 + 2*k)*n);
    Ax(5,2) = (C2beta*Ci*ipow(e,5)*(24*(-42*(5088 + x*(55883 - 11980*eta + 5440*ipow(Si,2)*(-1 + 3*eta))) + k*(140628 + x*(1510940 - 326836*eta + 153211*ipow(Si,2)*(-1 + 3*eta)))) + ipow(e,2)*(k*(-10591812 - 11637017*ipow(Si,2)*x*(-1 + 3*eta) + 16*x*(-8817583 + 2312949*eta)) + 21*(755748 + x*(10043750 - 2639004*eta + 830851*ipow(Si,2)*(-1 + 3*eta))))))/(27648.*(-3 + 2*k)*(7 + 2*k)*n);
    Bx(5,2) = (Ci*ipow(e,5)*S2beta*(96*k*(354324 + 19*x*(203054 - 43684*eta + 20089*ipow(Si,2)*(-1 + 3*eta))) - 504*(99924 + x*(1063876 - 230996*eta + 109691*ipow(Si,2)*(-1 + 3*eta))) + ipow(e,2)*(21*(7568820 + x*(100906328 - 26451168*eta + 8313613*ipow(Si,2)*(-1 + 3*eta))) - 8*k*(13231260 + x*(176000039 - 46213134*eta + 14542444*ipow(Si,2)*(-1 + 3*eta))))))/(276480.*(-3 + 2*k)*(7 + 2*k)*n);
    Cx(5,2) = (Ci*ipow(e,5)*S2beta*(24*(-42*(5088 + x*(55883 - 11980*eta + 5440*ipow(Si,2)*(-1 + 3*eta))) + k*(140628 + x*(1510940 - 326836*eta + 153211*ipow(Si,2)*(-1 + 3*eta)))) + ipow(e,2)*(k*(-10591812 - 11637017*ipow(Si,2)*x*(-1 + 3*eta) + 16*x*(-8817583 + 2312949*eta)) + 21*(755748 + x*(10043750 - 2639004*eta + 830851*ipow(Si,2)*(-1 + 3*eta))))))/(27648.*(-3 + 2*k)*(7 + 2*k)*n);
    Dx(5,2) = (C2beta*Ci*ipow(e,5)*(-96*k*(354324 + 19*x*(203054 - 43684*eta + 20089*ipow(Si,2)*(-1 + 3*eta))) + 504*(99924 + x*(1063876 - 230996*eta + 109691*ipow(Si,2)*(-1 + 3*eta))) + ipow(e,2)*(-21*(7568820 + x*(100906328 - 26451168*eta + 8313613*ipow(Si,2)*(-1 + 3*eta))) + 8*k*(13231260 + x*(176000039 - 46213134*eta + 14542444*ipow(Si,2)*(-1 + 3*eta))))))/(276480.*(-3 + 2*k)*(7 + 2*k)*n);
    Ap(5,3) = ((1 + ipow(Ci,2))*ipow(e,5)*(39384 - 58959*k + 16*ipow(e,2)*(-11672 + 17511*k))*S3beta*Si*sqrt(x)*delta)/(144.*(-2 + 3*k)*(8 + 3*k)*n);
    Bp(5,3) = -(C3beta*(1 + ipow(Ci,2))*ipow(e,5)*(196296 - 295029*k + 16*ipow(e,2)*(-58376 + 87549*k))*Si*sqrt(x)*delta)/(720.*(-2 + 3*k)*(8 + 3*k)*n);
    Cp(5,3) = -(C3beta*(1 + ipow(Ci,2))*ipow(e,5)*(39384 - 58959*k + 16*ipow(e,2)*(-11672 + 17511*k))*Si*sqrt(x)*delta)/(144.*(-2 + 3*k)*(8 + 3*k)*n);
    Dp(5,3) = -((1 + ipow(Ci,2))*ipow(e,5)*(196296 - 295029*k + 16*ipow(e,2)*(-58376 + 87549*k))*S3beta*Si*sqrt(x)*delta)/(720.*(-2 + 3*k)*(8 + 3*k)*n);
    Ax(5,3) = (C3beta*Ci*ipow(e,5)*(39384 - 58959*k + 16*ipow(e,2)*(-11672 + 17511*k))*Si*sqrt(x)*delta)/(72.*(-2 + 3*k)*(8 + 3*k)*n);
    Bx(5,3) = (Ci*ipow(e,5)*(196296 - 295029*k + 16*ipow(e,2)*(-58376 + 87549*k))*S3beta*Si*sqrt(x)*delta)/(360.*(-2 + 3*k)*(8 + 3*k)*n);
    Cx(5,3) = (Ci*ipow(e,5)*(39384 - 58959*k + 16*ipow(e,2)*(-11672 + 17511*k))*S3beta*Si*sqrt(x)*delta)/(72.*(-2 + 3*k)*(8 + 3*k)*n);
    Dx(5,3) = -(C3beta*Ci*ipow(e,5)*(196296 - 295029*k + 16*ipow(e,2)*(-58376 + 87549*k))*Si*sqrt(x)*delta)/(360.*(-2 + 3*k)*(8 + 3*k)*n);
    Ap(5,4) = ((1 + ipow(Ci,2))*ipow(e,5)*(51667848 - 206619024*k + ipow(e,2)*(-342215325 + 1368888442*k))*S4beta*ipow(Si,2)*x*(-1 + 3*eta))/(110592.*(-1 + 4*k)*(9 + 4*k)*n);
    Bp(5,4) = -(C4beta*(1 + ipow(Ci,2))*ipow(e,5)*(516442824 - 2066294976*k + ipow(e,2)*(-3422275389 + 13688830136*k))*ipow(Si,2)*x*(-1 + 3*eta))/(1.10592e6*(-1 + 4*k)*(9 + 4*k)*n);
    Cp(5,4) = -(C4beta*(1 + ipow(Ci,2))*ipow(e,5)*(51667848 - 206619024*k + ipow(e,2)*(-342215325 + 1368888442*k))*ipow(Si,2)*x*(-1 + 3*eta))/(110592.*(-1 + 4*k)*(9 + 4*k)*n);
    Dp(5,4) = -((1 + ipow(Ci,2))*ipow(e,5)*(516442824 - 2066294976*k + ipow(e,2)*(-3422275389 + 13688830136*k))*S4beta*ipow(Si,2)*x*(-1 + 3*eta))/(1.10592e6*(-1 + 4*k)*(9 + 4*k)*n);
    Ax(5,4) = (C4beta*Ci*ipow(e,5)*(51667848 - 206619024*k + ipow(e,2)*(-342215325 + 1368888442*k))*ipow(Si,2)*x*(-1 + 3*eta))/(55296.*(-1 + 4*k)*(9 + 4*k)*n);
    Bx(5,4) = (Ci*ipow(e,5)*(516442824 - 2066294976*k + ipow(e,2)*(-3422275389 + 13688830136*k))*S4beta*ipow(Si,2)*x*(-1 + 3*eta))/(552960.*(-1 + 4*k)*(9 + 4*k)*n);
    Cx(5,4) = (Ci*ipow(e,5)*(51667848 - 206619024*k + ipow(e,2)*(-342215325 + 1368888442*k))*S4beta*ipow(Si,2)*x*(-1 + 3*eta))/(55296.*(-1 + 4*k)*(9 + 4*k)*n);
    Dx(5,4) = -(C4beta*Ci*ipow(e,5)*(516442824 - 2066294976*k + ipow(e,2)*(-3422275389 + 13688830136*k))*ipow(Si,2)*x*(-1 + 3*eta))/(552960.*(-1 + 4*k)*(9 + 4*k)*n);
    Ap(6,0) = 0;
    Bp(6,0) = (-9*ipow(e,6)*(2*(-1 + ipow(Ci,2))*x*(-7*(117 + 59*eta) + 3*ipow(e,2)*(431 + 249*eta)) + 3*ipow(Si,2)*(-28*(1 + 3*(1 + ipow(Ci,2))*x*(-1 + 3*eta)) + 9*ipow(e,2)*(4 + 17*(1 + ipow(Ci,2))*x*(-1 + 3*eta)))))/(2240.*n);
    Cp(6,0) = 0;
    Dp(6,0) = 0;
    Ax(6,0) = (-81*Ci*ipow(e,6)*(-14 + 25*ipow(e,2))*ipow(Si,2)*x*(-1 + 3*eta))/(1120.*n);
    Bx(6,0) = 0;
    Cx(6,0) = 0;
    Dx(6,0) = 0;
    Ap(6,1) = (ipow(e,6)*(-93933490 + ipow(e,2)*(241825395 - 49101954*k) + 19757948*k + ipow(Ci,2)*(57333430 - 12280436*k + ipow(e,2)*(-153808585 + 31322342*k)))*S1beta*Si*sqrt(x)*delta)/(645120.*(-5 + k)*(7 + k)*n);
    Bp(6,1) = (C1beta*ipow(e,6)*(66399410 - 12632382*k + ipow(Ci,2)*(-41386870 + ipow(e,2)*(104719265 - 20570103*k) + 7734874*k) + 3*ipow(e,2)*(-54694185 + 10775087*k))*Si*sqrt(x)*delta)/(430080.*(-5 + k)*(7 + k)*n);
    Cp(6,1) = (C1beta*ipow(e,6)*(93933490 - 19757948*k + 3*ipow(e,2)*(-80608465 + 16367318*k) + ipow(Ci,2)*(ipow(e,2)*(153808585 - 31322342*k) + 14*(-4095245 + 877174*k)))*Si*sqrt(x)*delta)/(645120.*(-5 + k)*(7 + k)*n);
    Dp(6,1) = (ipow(e,6)*(66399410 - 12632382*k + ipow(Ci,2)*(-41386870 + ipow(e,2)*(104719265 - 20570103*k) + 7734874*k) + 3*ipow(e,2)*(-54694185 + 10775087*k))*S1beta*Si*sqrt(x)*delta)/(430080.*(-5 + k)*(7 + k)*n);
    Ax(6,1) = (C1beta*Ci*ipow(e,6)*(ipow(e,2)*(44008405 - 8889806*k) + 42*(-435715 + 89018*k))*Si*sqrt(x)*delta)/(322560.*(-5 + k)*(7 + k)*n);
    Bx(6,1) = (Ci*ipow(e,6)*(ipow(e,2)*(29681645 - 5877579*k) + 14*(-893305 + 174911*k))*S1beta*Si*sqrt(x)*delta)/(215040.*(-5 + k)*(7 + k)*n);
    Cx(6,1) = (Ci*ipow(e,6)*(ipow(e,2)*(44008405 - 8889806*k) + 42*(-435715 + 89018*k))*S1beta*Si*sqrt(x)*delta)/(322560.*(-5 + k)*(7 + k)*n);
    Dx(6,1) = (C1beta*Ci*ipow(e,6)*(12506270 - 2448754*k + ipow(e,2)*(-29681645 + 5877579*k))*Si*sqrt(x)*delta)/(215040.*(-5 + k)*(7 + k)*n);
    Ap(6,2) = -(S2beta*(-9*(6*(1 + ipow(Ci,2))*ipow(e,6)*(-19166 + 64165*ipow(e,2)) + x*(-1298248*ipow(e,6)*(-1 + 3*eta) + 5236780*ipow(e,8)*(-1 + 3*eta) - (1 + ipow(Ci,2))*(-7*ipow(e,6)*(-296049 + 324862*eta + 27284*ipow(Si,2)*(-1 + 3*eta)) + ipow(e,8)*(-8123013 + 9388295*eta + 849158*ipow(Si,2)*(-1 + 3*eta))))) + (1 + k)*(6*(1 + ipow(Ci,2))*ipow(e,6)*(-57190 + 192529*ipow(e,2)) + x*(-3904040*ipow(e,6)*(-1 + 3*eta) + 15714172*ipow(e,8)*(-1 + 3*eta) - (1 + ipow(Ci,2))*(-7*ipow(e,6)*(-879669 + 975110*eta + 81988*ipow(Si,2)*(-1 + 3*eta)) + ipow(e,8)*(-24380433 + 28171499*eta + 2547950*ipow(Si,2)*(-1 + 3*eta)))))))/(15120.*(-2 + k)*(4 + k)*n);
    Bp(6,2) = (C2beta*ipow(e,6)*(ipow(e,2)*(-770184 + 384990*k + k*x*(2886233 + 6322045*eta - 849158*ipow(Si,2)*(-1 + 3*eta)) + 12*x*(-481669 - 1054081*eta + 141566*ipow(Si,2)*(-1 + 3*eta))) + 7*(k*(-16428 + 27284*ipow(Si,2)*x*(-1 + 3*eta) - 5*x*(22117 + 46306*eta)) + 4*(8148 + x*(52841 + 116630*eta - 13676*ipow(Si,2)*(-1 + 3*eta)))) + ipow(Ci,2)*(ipow(e,2)*(-770184 + 384990*k + k*x*(8123013 - 9388295*eta - 849158*ipow(Si,2)*(-1 + 3*eta)) + 12*x*(-1354785 + 1565267*eta + 141566*ipow(Si,2)*(-1 + 3*eta))) + 7*(32592 + x*(583620 - 650248*eta - 54704*ipow(Si,2)*(-1 + 3*eta)) + k*(-16428 + x*(-296049 + 324862*eta + 27284*ipow(Si,2)*(-1 + 3*eta)))))))/(5040.*(-2 + k)*(4 + k)*n);
    Cp(6,2) = -(C2beta*(54*(1 + ipow(Ci,2))*ipow(e,6)*(-19166 + 64165*ipow(e,2)) - 9*x*(1298248*ipow(e,6)*(-1 + 3*eta) - 5236780*ipow(e,8)*(-1 + 3*eta) + (1 + ipow(Ci,2))*(-7*ipow(e,6)*(-296049 + 324862*eta + 27284*ipow(Si,2)*(-1 + 3*eta)) + ipow(e,8)*(-8123013 + 9388295*eta + 849158*ipow(Si,2)*(-1 + 3*eta)))) - (1 + k)*(6*(1 + ipow(Ci,2))*ipow(e,6)*(-57190 + 192529*ipow(e,2)) + x*(-3904040*ipow(e,6)*(-1 + 3*eta) + 15714172*ipow(e,8)*(-1 + 3*eta) - (1 + ipow(Ci,2))*(-7*ipow(e,6)*(-879669 + 975110*eta + 81988*ipow(Si,2)*(-1 + 3*eta)) + ipow(e,8)*(-24380433 + 28171499*eta + 2547950*ipow(Si,2)*(-1 + 3*eta)))))))/(15120.*(-2 + k)*(4 + k)*n);
    Dp(6,2) = (ipow(e,6)*S2beta*(ipow(e,2)*(-770184 + 384990*k + k*x*(2886233 + 6322045*eta - 849158*ipow(Si,2)*(-1 + 3*eta)) + 12*x*(-481669 - 1054081*eta + 141566*ipow(Si,2)*(-1 + 3*eta))) + 7*(k*(-16428 + 27284*ipow(Si,2)*x*(-1 + 3*eta) - 5*x*(22117 + 46306*eta)) + 4*(8148 + x*(52841 + 116630*eta - 13676*ipow(Si,2)*(-1 + 3*eta)))) + ipow(Ci,2)*(ipow(e,2)*(-770184 + 384990*k + k*x*(8123013 - 9388295*eta - 849158*ipow(Si,2)*(-1 + 3*eta)) + 12*x*(-1354785 + 1565267*eta + 141566*ipow(Si,2)*(-1 + 3*eta))) + 7*(32592 + x*(583620 - 650248*eta - 54704*ipow(Si,2)*(-1 + 3*eta)) + k*(-16428 + x*(-296049 + 324862*eta + 27284*ipow(Si,2)*(-1 + 3*eta)))))))/(5040.*(-2 + k)*(4 + k)*n);
    Ax(6,2) = (C2beta*Ci*ipow(e,6)*(7*(-98832 + 49020*k - 4*x*(307261 - 70366*eta + 28574*ipow(Si,2)*(-1 + 3*eta)) + k*x*(600809 - 138530*eta + 57442*ipow(Si,2)*(-1 + 3*eta))) + ipow(e,2)*(2309736 - 1155174*k + k*x*(-16523347 + 4600241*eta - 1380593*ipow(Si,2)*(-1 + 3*eta)) + 4*x*(8254565 - 2299471*eta + 689935*ipow(Si,2)*(-1 + 3*eta)))))/(7560.*(-2 + k)*(4 + k)*n);
    Bx(6,2) = (Ci*ipow(e,6)*S2beta*(ipow(e,2)*(770184 - 384990*k + k*x*(-5504623 + 1533125*eta - 460037*ipow(Si,2)*(-1 + 3*eta)) + 12*x*(918227 - 255593*eta + 76713*ipow(Si,2)*(-1 + 3*eta))) + 7*(-4*(8148 + x*(99373 - 22966*eta + 9590*ipow(Si,2)*(-1 + 3*eta))) + k*(16428 + x*(203317 - 46666*eta + 19082*ipow(Si,2)*(-1 + 3*eta))))))/(2520.*(-2 + k)*(4 + k)*n);
    Cx(6,2) = (Ci*ipow(e,6)*S2beta*(7*(-98832 + 49020*k - 4*x*(307261 - 70366*eta + 28574*ipow(Si,2)*(-1 + 3*eta)) + k*x*(600809 - 138530*eta + 57442*ipow(Si,2)*(-1 + 3*eta))) + ipow(e,2)*(2309736 - 1155174*k + k*x*(-16523347 + 4600241*eta - 1380593*ipow(Si,2)*(-1 + 3*eta)) + 4*x*(8254565 - 2299471*eta + 689935*ipow(Si,2)*(-1 + 3*eta)))))/(7560.*(-2 + k)*(4 + k)*n);
    Dx(6,2) = (C2beta*Ci*ipow(e,6)*(28*(8148 + x*(99373 - 22966*eta + 9590*ipow(Si,2)*(-1 + 3*eta))) - 7*k*(16428 + x*(203317 - 46666*eta + 19082*ipow(Si,2)*(-1 + 3*eta))) + ipow(e,2)*(-12*(64182 + x*(918227 - 255593*eta + 76713*ipow(Si,2)*(-1 + 3*eta))) + k*(384990 + x*(5504623 - 1533125*eta + 460037*ipow(Si,2)*(-1 + 3*eta))))))/(2520.*(-2 + k)*(4 + k)*n);
    Ap(6,3) = (3*(1 + ipow(Ci,2))*ipow(e,6)*(1861482 - 1859564*k + ipow(e,2)*(-9166887 + 9167514*k))*S3beta*Si*sqrt(x)*delta)/(71680.*(-1 + k)*(3 + k)*n);
    Bp(6,3) = (-3*C3beta*(1 + ipow(Ci,2))*ipow(e,6)*(3717210 - 3721046*k + 3*ipow(e,2)*(-6111885 + 6111467*k))*Si*sqrt(x)*delta)/(143360.*(-1 + k)*(3 + k)*n);
    Cp(6,3) = (-3*C3beta*(1 + ipow(Ci,2))*ipow(e,6)*(1861482 - 1859564*k + ipow(e,2)*(-9166887 + 9167514*k))*Si*sqrt(x)*delta)/(71680.*(-1 + k)*(3 + k)*n);
    Dp(6,3) = (-3*(1 + ipow(Ci,2))*ipow(e,6)*(3717210 - 3721046*k + 3*ipow(e,2)*(-6111885 + 6111467*k))*S3beta*Si*sqrt(x)*delta)/(143360.*(-1 + k)*(3 + k)*n);
    Ax(6,3) = (3*C3beta*Ci*ipow(e,6)*(1861482 - 1859564*k + ipow(e,2)*(-9166887 + 9167514*k))*Si*sqrt(x)*delta)/(35840.*(-1 + k)*(3 + k)*n);
    Bx(6,3) = (3*Ci*ipow(e,6)*(3717210 - 3721046*k + 3*ipow(e,2)*(-6111885 + 6111467*k))*S3beta*Si*sqrt(x)*delta)/(71680.*(-1 + k)*(3 + k)*n);
    Cx(6,3) = (3*Ci*ipow(e,6)*(1861482 - 1859564*k + ipow(e,2)*(-9166887 + 9167514*k))*S3beta*Si*sqrt(x)*delta)/(35840.*(-1 + k)*(3 + k)*n);
    Dx(6,3) = (3*C3beta*Ci*ipow(e,6)*(-3717210 + 3721046*k - 3*ipow(e,2)*(-6111885 + 6111467*k))*Si*sqrt(x)*delta)/(71680.*(-1 + k)*(3 + k)*n);
    Ap(6,4) = ((1 + ipow(Ci,2))*ipow(e,6)*(54692470 - 109373012*k + ipow(e,2)*(-367185665 + 734375734*k))*S4beta*ipow(Si,2)*x*(-1 + 3*eta))/(120960.*(-1 + 2*k)*(5 + 2*k)*n);
    Bp(6,4) = -(C4beta*(1 + ipow(Ci,2))*ipow(e,6)*(18227510 - 36458996*k + ipow(e,2)*(-122396445 + 244791422*k))*ipow(Si,2)*x*(-1 + 3*eta))/(40320.*(-1 + 2*k)*(5 + 2*k)*n);
    Cp(6,4) = -(C4beta*(1 + ipow(Ci,2))*ipow(e,6)*(54692470 - 109373012*k + ipow(e,2)*(-367185665 + 734375734*k))*ipow(Si,2)*x*(-1 + 3*eta))/(120960.*(-1 + 2*k)*(5 + 2*k)*n);
    Dp(6,4) = -((1 + ipow(Ci,2))*ipow(e,6)*(18227510 - 36458996*k + ipow(e,2)*(-122396445 + 244791422*k))*S4beta*ipow(Si,2)*x*(-1 + 3*eta))/(40320.*(-1 + 2*k)*(5 + 2*k)*n);
    Ax(6,4) = (C4beta*Ci*ipow(e,6)*(54692470 - 109373012*k + ipow(e,2)*(-367185665 + 734375734*k))*ipow(Si,2)*x*(-1 + 3*eta))/(60480.*(-1 + 2*k)*(5 + 2*k)*n);
    Bx(6,4) = (Ci*ipow(e,6)*(18227510 - 36458996*k + ipow(e,2)*(-122396445 + 244791422*k))*S4beta*ipow(Si,2)*x*(-1 + 3*eta))/(20160.*(-1 + 2*k)*(5 + 2*k)*n);
    Cx(6,4) = (Ci*ipow(e,6)*(54692470 - 109373012*k + ipow(e,2)*(-367185665 + 734375734*k))*S4beta*ipow(Si,2)*x*(-1 + 3*eta))/(60480.*(-1 + 2*k)*(5 + 2*k)*n);
    Dx(6,4) = -(C4beta*Ci*ipow(e,6)*(18227510 - 36458996*k + ipow(e,2)*(-122396445 + 244791422*k))*ipow(Si,2)*x*(-1 + 3*eta))/(20160.*(-1 + 2*k)*(5 + 2*k)*n);
    Ap(7,0) = 0;
    Bp(7,0) = (16807*ipow(e,7)*(4*(-1 + ipow(Ci,2))*x*(135 + 71*eta) + 3*ipow(Si,2)*(8 + 29*(1 + ipow(Ci,2))*x*(-1 + 3*eta))))/(1.10592e6*n);
    Cp(7,0) = 0;
    Dp(7,0) = 0;
    Ax(7,0) = (117649*Ci*ipow(e,7)*ipow(Si,2)*x*(-1 + 3*eta))/(92160.*n);
    Bx(7,0) = 0;
    Cx(7,0) = 0;
    Dx(7,0) = 0;
    Ap(7,1) = (ipow(e,7)*(-370008 + ipow(Ci,2)*(231240 - 40727*k) + 64341*k)*S1beta*Si*sqrt(x)*delta)/(1440.*(-6 + k)*(8 + k)*n);
    Bp(7,1) = (C1beta*ipow(e,7)*(2718360 - 434349*k + ipow(Ci,2)*(-1723656 + 271967*k))*Si*sqrt(x)*delta)/(10080.*(-6 + k)*(8 + k)*n);
    Cp(7,1) = (C1beta*ipow(e,7)*(370008 - 64341*k + ipow(Ci,2)*(-231240 + 40727*k))*Si*sqrt(x)*delta)/(1440.*(-6 + k)*(8 + k)*n);
    Dp(7,1) = (ipow(e,7)*(2718360 - 434349*k + ipow(Ci,2)*(-1723656 + 271967*k))*S1beta*Si*sqrt(x)*delta)/(10080.*(-6 + k)*(8 + k)*n);
    Ax(7,1) = (C1beta*Ci*ipow(e,7)*(-69384 + 11807*k)*Si*sqrt(x)*delta)/(720.*(-6 + k)*(8 + k)*n);
    Bx(7,1) = (Ci*ipow(e,7)*(-497352 + 81191*k)*S1beta*Si*sqrt(x)*delta)/(5040.*(-6 + k)*(8 + k)*n);
    Cx(7,1) = (Ci*ipow(e,7)*(-69384 + 11807*k)*S1beta*Si*sqrt(x)*delta)/(720.*(-6 + k)*(8 + k)*n);
    Dx(7,1) = (C1beta*Ci*ipow(e,7)*(497352 - 81191*k)*Si*sqrt(x)*delta)/(5040.*(-6 + k)*(8 + k)*n);
    Ap(7,2) = (ipow(e,7)*S2beta*(45*(-2057076 + 4096027*ipow(Si,2)*x*(-1 + 3*eta) - 2*x*(7516757 + 16513656*eta)) + k*(36824868 + x*(261888602 + 597622866*eta - 73830986*ipow(Si,2)*(-1 + 3*eta))) + ipow(Ci,2)*(45*(-2057076 + 4096027*ipow(Si,2)*x*(-1 + 3*eta) + 42*x*(-991381 + 1113956*eta)) + k*(36824868 + x*(741986886 - 842671986*eta - 73830986*ipow(Si,2)*(-1 + 3*eta))))))/(276480.*(-5 + 2*k)*(9 + 2*k)*n);
    Bp(7,2) = (C2beta*ipow(e,7)*(643422690 - 45*x*(-100877273 - 232756809*eta + 28723439*ipow(Si,2)*(-1 + 3*eta)) + 2*k*(-129393288 + x*(-938396732 - 2083851906*eta + 258152201*ipow(Si,2)*(-1 + 3*eta))) + ipow(Ci,2)*(643422690 - 45*x*(-287717439 + 327763689*eta + 28723439*ipow(Si,2)*(-1 + 3*eta)) + 2*k*(-129393288 + 258152201*ipow(Si,2)*x*(-1 + 3*eta) + 6*x*(-435949496 + 491341471*eta)))))/(1.93536e6*(-5 + 2*k)*(9 + 2*k)*n);
    Cp(7,2) = (C2beta*ipow(e,7)*(92568420 - 36824868*k + 45*x*(15033514 + 33027312*eta - 4096027*ipow(Si,2)*(-1 + 3*eta)) + k*x*(-261888602 - 597622866*eta + 73830986*ipow(Si,2)*(-1 + 3*eta)) + ipow(Ci,2)*(-45*(-2057076 + 4096027*ipow(Si,2)*x*(-1 + 3*eta) + 42*x*(-991381 + 1113956*eta)) + k*(-36824868 + x*(-741986886 + 842671986*eta + 73830986*ipow(Si,2)*(-1 + 3*eta))))))/(276480.*(-5 + 2*k)*(9 + 2*k)*n);
    Dp(7,2) = (ipow(e,7)*S2beta*(643422690 - 45*x*(-100877273 - 232756809*eta + 28723439*ipow(Si,2)*(-1 + 3*eta)) + 2*k*(-129393288 + x*(-938396732 - 2083851906*eta + 258152201*ipow(Si,2)*(-1 + 3*eta))) + ipow(Ci,2)*(643422690 - 45*x*(-287717439 + 327763689*eta + 28723439*ipow(Si,2)*(-1 + 3*eta)) + 2*k*(-129393288 + 258152201*ipow(Si,2)*x*(-1 + 3*eta) + 6*x*(-435949496 + 491341471*eta)))))/(1.93536e6*(-5 + 2*k)*(9 + 2*k)*n);
    Ax(7,2) = (C2beta*Ci*ipow(e,7)*(-45*(2057076 + x*(28335758 - 6879420*eta + 2555095*ipow(Si,2)*(-1 + 3*eta))) + k*(36824868 + x*(501937744 - 122524560*eta + 46193585*ipow(Si,2)*(-1 + 3*eta)))))/(138240.*(-5 + 2*k)*(9 + 2*k)*n);
    Bx(7,2) = (Ci*ipow(e,7)*S2beta*(-45*(28596564 + x*(388594712 - 95006880*eta + 35973205*ipow(Si,2)*(-1 + 3*eta))) + 8*k*(64696644 + x*(888523427 - 216049230*eta + 80586430*ipow(Si,2)*(-1 + 3*eta)))))/(1.93536e6*(-5 + 2*k)*(9 + 2*k)*n);
    Cx(7,2) = (Ci*ipow(e,7)*S2beta*(-45*(2057076 + x*(28335758 - 6879420*eta + 2555095*ipow(Si,2)*(-1 + 3*eta))) + k*(36824868 + x*(501937744 - 122524560*eta + 46193585*ipow(Si,2)*(-1 + 3*eta)))))/(138240.*(-5 + 2*k)*(9 + 2*k)*n);
    Dx(7,2) = (C2beta*Ci*ipow(e,7)*(45*(28596564 + x*(388594712 - 95006880*eta + 35973205*ipow(Si,2)*(-1 + 3*eta))) - 8*k*(64696644 + x*(888523427 - 216049230*eta + 80586430*ipow(Si,2)*(-1 + 3*eta)))))/(1.93536e6*(-5 + 2*k)*(9 + 2*k)*n);
    Ap(7,3) = -((1 + ipow(Ci,2))*ipow(e,7)*(-1116560 + 836907*k)*S3beta*Si*sqrt(x)*delta)/(720.*(-4 + 3*k)*(10 + 3*k)*n);
    Bp(7,3) = (C3beta*(1 + ipow(Ci,2))*ipow(e,7)*(-7809080 + 5860401*k)*Si*sqrt(x)*delta)/(5040.*(-4 + 3*k)*(10 + 3*k)*n);
    Cp(7,3) = (C3beta*(1 + ipow(Ci,2))*ipow(e,7)*(-1116560 + 836907*k)*Si*sqrt(x)*delta)/(720.*(-4 + 3*k)*(10 + 3*k)*n);
    Dp(7,3) = ((1 + ipow(Ci,2))*ipow(e,7)*(-7809080 + 5860401*k)*S3beta*Si*sqrt(x)*delta)/(5040.*(-4 + 3*k)*(10 + 3*k)*n);
    Ax(7,3) = (C3beta*Ci*ipow(e,7)*(1116560 - 836907*k)*Si*sqrt(x)*delta)/(360.*(-4 + 3*k)*(10 + 3*k)*n);
    Bx(7,3) = (Ci*ipow(e,7)*(7809080 - 5860401*k)*S3beta*Si*sqrt(x)*delta)/(2520.*(-4 + 3*k)*(10 + 3*k)*n);
    Cx(7,3) = (Ci*ipow(e,7)*(1116560 - 836907*k)*S3beta*Si*sqrt(x)*delta)/(360.*(-4 + 3*k)*(10 + 3*k)*n);
    Dx(7,3) = (C3beta*Ci*ipow(e,7)*(-7809080 + 5860401*k)*Si*sqrt(x)*delta)/(2520.*(-4 + 3*k)*(10 + 3*k)*n);
    Ap(7,4) = -((1 + ipow(Ci,2))*ipow(e,7)*(-2779141233 + 3705298570*k)*S4beta*ipow(Si,2)*x*(-1 + 3*eta))/(552960.*(-3 + 4*k)*(11 + 4*k)*n);
    Bp(7,4) = (C4beta*(1 + ipow(Ci,2))*ipow(e,7)*(-38904296541 + 51875518424*k)*ipow(Si,2)*x*(-1 + 3*eta))/(7.74144e6*(-3 + 4*k)*(11 + 4*k)*n);
    Cp(7,4) = (C4beta*(1 + ipow(Ci,2))*ipow(e,7)*(-2779141233 + 3705298570*k)*ipow(Si,2)*x*(-1 + 3*eta))/(552960.*(-3 + 4*k)*(11 + 4*k)*n);
    Dp(7,4) = ((1 + ipow(Ci,2))*ipow(e,7)*(-38904296541 + 51875518424*k)*S4beta*ipow(Si,2)*x*(-1 + 3*eta))/(7.74144e6*(-3 + 4*k)*(11 + 4*k)*n);
    Ax(7,4) = -(C4beta*Ci*ipow(e,7)*(-2779141233 + 3705298570*k)*ipow(Si,2)*x*(-1 + 3*eta))/(276480.*(-3 + 4*k)*(11 + 4*k)*n);
    Bx(7,4) = -(Ci*ipow(e,7)*(-38904296541 + 51875518424*k)*S4beta*ipow(Si,2)*x*(-1 + 3*eta))/(3.87072e6*(-3 + 4*k)*(11 + 4*k)*n);
    Cx(7,4) = -(Ci*ipow(e,7)*(-2779141233 + 3705298570*k)*S4beta*ipow(Si,2)*x*(-1 + 3*eta))/(276480.*(-3 + 4*k)*(11 + 4*k)*n);
    Dx(7,4) = (C4beta*Ci*ipow(e,7)*(-38904296541 + 51875518424*k)*ipow(Si,2)*x*(-1 + 3*eta))/(3.87072e6*(-3 + 4*k)*(11 + 4*k)*n);
    Ap(8,0) = 0;
    Bp(8,0) = (32*ipow(e,8)*(2*(-1 + ipow(Ci,2))*x*(153 + 83*eta) + 3*ipow(Si,2)*(4 + 17*(1 + ipow(Ci,2))*x*(-1 + 3*eta))))/(945.*n);
    Cp(8,0) = 0;
    Dp(8,0) = 0;
    Ax(8,0) = (512*Ci*ipow(e,8)*ipow(Si,2)*x*(-1 + 3*eta))/(315.*n);
    Bx(8,0) = 0;
    Cx(8,0) = 0;
    Dx(8,0) = 0;
    Ap(8,1) = (ipow(e,8)*(-9001471815 + ipow(Ci,2)*(5723933733 - 856293691*k) + 1333925337*k)*S1beta*Si*sqrt(x)*delta)/(2.064384e7*(-7 + k)*(9 + k)*n);
    Bp(8,1) = (C1beta*ipow(e,8)*(2344869513 - 322981161*k + ipow(Ci,2)*(-1506955275 + 205632107*k))*Si*sqrt(x)*delta)/(5.16096e6*(-7 + k)*(9 + k)*n);
    Cp(8,1) = (C1beta*ipow(e,8)*(9001471815 - 1333925337*k + ipow(Ci,2)*(-5723933733 + 856293691*k))*Si*sqrt(x)*delta)/(2.064384e7*(-7 + k)*(9 + k)*n);
    Dp(8,1) = (ipow(e,8)*(2344869513 - 322981161*k + ipow(Ci,2)*(-1506955275 + 205632107*k))*S1beta*Si*sqrt(x)*delta)/(5.16096e6*(-7 + k)*(9 + k)*n);
    Ax(8,1) = (C1beta*Ci*ipow(e,8)*(-1638769041 + 238815823*k)*Si*sqrt(x)*delta)/(1.032192e7*(-7 + k)*(9 + k)*n);
    Bx(8,1) = (Ci*ipow(e,8)*(-418957119 + 58674527*k)*S1beta*Si*sqrt(x)*delta)/(2.58048e6*(-7 + k)*(9 + k)*n);
    Cx(8,1) = (Ci*ipow(e,8)*(-1638769041 + 238815823*k)*S1beta*Si*sqrt(x)*delta)/(1.032192e7*(-7 + k)*(9 + k)*n);
    Dx(8,1) = (C1beta*Ci*ipow(e,8)*(418957119 - 58674527*k)*Si*sqrt(x)*delta)/(2.58048e6*(-7 + k)*(9 + k)*n);
    Ap(8,2) = (ipow(e,8)*S2beta*(-70487460 + 23402508*k + k*x*(176629895 + 426405274*eta - 54704996*ipow(Si,2)*(-1 + 3*eta)) + 15*x*(-36239335 - 84948242*eta + 10931668*ipow(Si,2)*(-1 + 3*eta)) + ipow(Ci,2)*(-70487460 + 23402508*k + k*x*(520624839 - 605579558*eta - 54704996*ipow(Si,2)*(-1 + 3*eta)) + 15*x*(-104907687 + 121056814*eta + 10931668*ipow(Si,2)*(-1 + 3*eta)))))/(483840.*(-3 + k)*(5 + k)*n);
    Bp(8,2) = (C2beta*ipow(e,8)*(35068770 - 11736246*k + 15*x*(17548820 + 42682129*eta - 5471666*ipow(Si,2)*(-1 + 3*eta)) + k*x*(-90027490 - 212578613*eta + 27335002*ipow(Si,2)*(-1 + 3*eta)) + ipow(Ci,2)*(35068770 - 11736246*k + 15*x*(51964644 - 60565343*eta - 5471666*ipow(Si,2)*(-1 + 3*eta)) + k*x*(-261780018 + 302678971*eta + 27335002*ipow(Si,2)*(-1 + 3*eta)))))/(241920.*(-3 + k)*(5 + k)*n);
    Cp(8,2) = (C2beta*ipow(e,8)*(70487460 - 23402508*k + 15*x*(36239335 + 84948242*eta - 10931668*ipow(Si,2)*(-1 + 3*eta)) + k*x*(-176629895 - 426405274*eta + 54704996*ipow(Si,2)*(-1 + 3*eta)) + ipow(Ci,2)*(70487460 - 23402508*k + 15*x*(104907687 - 121056814*eta - 10931668*ipow(Si,2)*(-1 + 3*eta)) + k*x*(-520624839 + 605579558*eta + 54704996*ipow(Si,2)*(-1 + 3*eta)))))/(483840.*(-3 + k)*(5 + k)*n);
    Dp(8,2) = (ipow(e,8)*S2beta*(35068770 - 11736246*k + 15*x*(17548820 + 42682129*eta - 5471666*ipow(Si,2)*(-1 + 3*eta)) + k*x*(-90027490 - 212578613*eta + 27335002*ipow(Si,2)*(-1 + 3*eta)) + ipow(Ci,2)*(35068770 - 11736246*k + 15*x*(51964644 - 60565343*eta - 5471666*ipow(Si,2)*(-1 + 3*eta)) + k*x*(-261780018 + 302678971*eta + 27335002*ipow(Si,2)*(-1 + 3*eta)))))/(241920.*(-3 + k)*(5 + k)*n);
    Ax(8,2) = (C2beta*Ci*ipow(e,8)*(-15*(4699164 + x*(70573511 - 18054286*eta + 6235420*ipow(Si,2)*(-1 + 3*eta))) + k*(23402508 + x*(348627367 - 89587142*eta + 31293740*ipow(Si,2)*(-1 + 3*eta)))))/(241920.*(-3 + k)*(5 + k)*n);
    Bx(8,2) = (Ci*ipow(e,8)*S2beta*(-15*(2337918 + x*(34756732 - 8941607*eta + 3132290*ipow(Si,2)*(-1 + 3*eta))) + k*(11736246 + x*(175903754 - 45050179*eta + 15603130*ipow(Si,2)*(-1 + 3*eta)))))/(120960.*(-3 + k)*(5 + k)*n);
    Cx(8,2) = (Ci*ipow(e,8)*S2beta*(-15*(4699164 + x*(70573511 - 18054286*eta + 6235420*ipow(Si,2)*(-1 + 3*eta))) + k*(23402508 + x*(348627367 - 89587142*eta + 31293740*ipow(Si,2)*(-1 + 3*eta)))))/(241920.*(-3 + k)*(5 + k)*n);
    Dx(8,2) = (C2beta*Ci*ipow(e,8)*(35068770 - 11736246*k + k*x*(-175903754 + 45050179*eta - 15603130*ipow(Si,2)*(-1 + 3*eta)) + 15*x*(34756732 - 8941607*eta + 3132290*ipow(Si,2)*(-1 + 3*eta))))/(120960.*(-3 + k)*(5 + k)*n);
    Ap(8,3) = -((1 + ipow(Ci,2))*ipow(e,8)*(-21620404355 + 12967107613*k)*S3beta*Si*sqrt(x)*delta)/(6.88128e6*(-5 + 3*k)*(11 + 3*k)*n);
    Bp(8,3) = (C3beta*(1 + ipow(Ci,2))*ipow(e,8)*(-16206477485 + 9727737741*k)*Si*sqrt(x)*delta)/(5.16096e6*(-5 + 3*k)*(11 + 3*k)*n);
    Cp(8,3) = (C3beta*(1 + ipow(Ci,2))*ipow(e,8)*(-21620404355 + 12967107613*k)*Si*sqrt(x)*delta)/(6.88128e6*(-5 + 3*k)*(11 + 3*k)*n);
    Dp(8,3) = ((1 + ipow(Ci,2))*ipow(e,8)*(-16206477485 + 9727737741*k)*S3beta*Si*sqrt(x)*delta)/(5.16096e6*(-5 + 3*k)*(11 + 3*k)*n);
    Ax(8,3) = (C3beta*Ci*ipow(e,8)*(21620404355 - 12967107613*k)*Si*sqrt(x)*delta)/(3.44064e6*(-5 + 3*k)*(11 + 3*k)*n);
    Bx(8,3) = (Ci*ipow(e,8)*(16206477485 - 9727737741*k)*S3beta*Si*sqrt(x)*delta)/(2.58048e6*(-5 + 3*k)*(11 + 3*k)*n);
    Cx(8,3) = (Ci*ipow(e,8)*(21620404355 - 12967107613*k)*S3beta*Si*sqrt(x)*delta)/(3.44064e6*(-5 + 3*k)*(11 + 3*k)*n);
    Dx(8,3) = (C3beta*Ci*ipow(e,8)*(-16206477485 + 9727737741*k)*Si*sqrt(x)*delta)/(2.58048e6*(-5 + 3*k)*(11 + 3*k)*n);
    Ap(8,4) = -((1 + ipow(Ci,2))*ipow(e,8)*(-90701751 + 90698435*k)*S4beta*ipow(Si,2)*x*(-1 + 3*eta))/(120960.*(-1 + k)*(3 + k)*n);
    Bp(8,4) = (C4beta*(1 + ipow(Ci,2))*ipow(e,8)*(-90696777 + 90700093*k)*ipow(Si,2)*x*(-1 + 3*eta))/(120960.*(-1 + k)*(3 + k)*n);
    Cp(8,4) = (C4beta*(1 + ipow(Ci,2))*ipow(e,8)*(-90701751 + 90698435*k)*ipow(Si,2)*x*(-1 + 3*eta))/(120960.*(-1 + k)*(3 + k)*n);
    Dp(8,4) = ((1 + ipow(Ci,2))*ipow(e,8)*(-90696777 + 90700093*k)*S4beta*ipow(Si,2)*x*(-1 + 3*eta))/(120960.*(-1 + k)*(3 + k)*n);
    Ax(8,4) = -(C4beta*Ci*ipow(e,8)*(-90701751 + 90698435*k)*ipow(Si,2)*x*(-1 + 3*eta))/(60480.*(-1 + k)*(3 + k)*n);
    Bx(8,4) = -(Ci*ipow(e,8)*(-90696777 + 90700093*k)*S4beta*ipow(Si,2)*x*(-1 + 3*eta))/(60480.*(-1 + k)*(3 + k)*n);
    Cx(8,4) = -(Ci*ipow(e,8)*(-90701751 + 90698435*k)*S4beta*ipow(Si,2)*x*(-1 + 3*eta))/(60480.*(-1 + k)*(3 + k)*n);
    Dx(8,4) = (C4beta*Ci*ipow(e,8)*(-90696777 + 90700093*k)*ipow(Si,2)*x*(-1 + 3*eta))/(60480.*(-1 + k)*(3 + k)*n);
    /* End of auto-generated code */
    
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

