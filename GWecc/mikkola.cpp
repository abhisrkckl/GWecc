#include <math.h>
#include <stdio.h>
#include "mikkola.h"

template <typename T> 
constexpr int sign(T val) {
    return (T(0) < val) - (val < T(0));
}

/*
 * This is based on [Mikkola 1987].
 *    adsabs.harvard.edu/full/1987CeMec..40..329M
 *
*/


/*
 * Solves the elliptical Kepler equation.
 *
 * l    is mean anomaly    in rad
 * e    is eccentricity
 *
 * Returns the eccentric anomaly u which is the solution of the Kepler equation.
 *     l = u - e sin(u)
 */
double MIKKOLA(double l, const double e){

    constexpr double Pi=M_PI, TwoPi=2*M_PI;

    if(e==0 || l==0){
        return l;
    }

    if(e<0 || e>=1){
        fprintf(stderr,"ERROR: The eccentricity of an ellipse must lie within [0,1).\n");
        return NAN;
    }

    int sgn = sign(l);
    l*=sgn;            // l>0
    
    int ncycles = (int)(l/TwoPi);
    l -= TwoPi*ncycles;    // 0<=l<2*pi

    bool flag = l>Pi;
    if(flag){
        l = TwoPi-l;    // 0<=l<=pi
    }    

    double  alpha  = (1-e)/(4*e + 0.5),
            alpha3 = alpha*alpha*alpha;
    double  beta   = (l/2)/(4*e + 0.5),
            beta2  = beta*beta;
    
    double z = (beta>0) ? cbrt(beta + sqrt(alpha3 + beta2))
                        : cbrt(beta - sqrt(alpha3 + beta2));

    double  s = (z - alpha/z),
            s5= s*s*s*s*s;
    double  w = (s - (0.078*s5)/(1 + e)),
            w3= w*w*w;
    
    double E0 = (l + e*(3*w - 4*w3));
    double u = E0 ;
    
    double  esu  = e*sin(u);
    double  ecu  = e*cos(u);
    
    double  fu  = (u - esu - l),
            f1u = (1 - ecu),
            f2u = (esu),
            f3u = (ecu),
            f4u =-(esu);

    double  u1 = -fu/ f1u,
            u2 = -fu/(f1u + f2u*u1/2),
            u3 = -fu/(f1u + f2u*u2/2 + f3u*(u2*u2)/6.0),
            u4 = -fu/(f1u + f2u*u3/2 + f3u*(u3*u3)/6.0 + f4u*(u3*u3*u3)/24.0),
            xi = (E0 + u4);
    
    double sol = flag    ?(TwoPi-xi)    :xi;
    
    u = sgn*(sol + ncycles*TwoPi);
    
    return u;
}

std::vector<double> MIKKOLA(const std::vector<double> &ls, const double e){
    
    size_t length = ls.size();
    std::vector<double> us(length);

    for(unsigned int i=0; i<length; i++){
        us[i] = MIKKOLA(ls[i],e);
    }

    return us;
}

std::vector<double> MIKKOLA(const std::vector<double> &ls, const std::vector<double> &es){
    
    size_t length = ls.size();

    if(es.size() != length){
        fprintf(stderr,"ERROR: Vector size mismatch. Returning empty vector.");
        std::vector<double> us(0);
        return us;
    }
    else{
        std::vector<double> us(length);

        for(unsigned int i=0; i<length; i++){
            us[i] = MIKKOLA(ls[i],es[i]);
        }
        return us;
    }
}
