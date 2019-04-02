#ifndef _Binary_hpp_
#define _Binary_hpp_ 1

#include <cmath>

/*
 * t0		in s
 * Omega	in rad
 * i		in rad
 * n 		in rad/s
 * l		in rad
 * gamma	in rad 
 */
struct BinaryState{
	double t;
	double Omega, i;
	double n, e, l, gamma;
};

/*
 * M 		in s (geometric units)
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
		  Mchirp(M * pow(eta, 3./5)) 	{}
	
	BinaryMass() 			  = delete;
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
 * DL 		in s (geometric units)
 * RA, DEC 	in rad
 */
struct SkyPosition{
	double DL;
	double RA, DEC;
};

#endif
