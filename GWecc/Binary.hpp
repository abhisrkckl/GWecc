#ifndef _Binary_hpp_
#define _Binary_hpp_ 1

#include <cmath>

/*
 * t0       in s
 * psi      in rad
 * i        in rad
 * n        in rad/s
 * l        in rad
 * gamma    in rad 
 */
struct BinaryState{
    double t;
    double psi, i;
    double n, e, l, gamma;
    bool merged;

    BinaryState(double _t, double _psi, double _i, double _n, double _e, double _l, double _gamma)
       : t(_t), psi(_psi), i(_i), n(_n), e(_e), l(_l), gamma(_gamma), merged(false) {}

    BinaryState()
       : t(NAN), psi(NAN), i(NAN), n(NAN), e(NAN), l(NAN), gamma(NAN), 
         merged(true) {} 
};

/*
 * M         in s (geometric units)
 */
class BinaryMass{

private:
    const double M, q;
    const double eta, delta, Mchirp;

public:
    BinaryMass(double _M, double _q) 
        : M(_M), 
          q(_q),
          eta(q/(1+q)/(1+q)),
          delta(fabs(1-q)/(1+q)),
          Mchirp(M * pow(eta, 3./5)) {}
    
    BinaryMass() = delete;
    void operator=(BinaryMass const&) = delete;
    
    const double mass() const{
        return M;
    }
    const double mass_ratio() const{
        return q;
    }
    const double symmetric_mass_ratio() const{
        return eta;
    }
    const double differential_mass_ratio() const{
        return delta;
    }
    const double chirp_mass() const{
        return Mchirp;
    }
};

/*
 * DL         in s (geometric units)
 * RA, DEC     in rad
 */
struct SkyPosition{
    double DL;
    double RA, DEC;
    double z;
    
    SkyPosition(double _DL, double _RA, double _DEC, double _z) : DL(_DL), RA(_RA), DEC(_DEC), z(_z) {} 
};

#endif
