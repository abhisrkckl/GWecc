#ifndef _ipow_hpp_
#define _ipow_hpp_ 1

constexpr static double ipow(const double x, const unsigned int n){
    
    double result{1}, temp{x};
    
    for(unsigned int j=1; j<=n; j<<=1){
       if(n&j){
           result *= temp;
       }
       temp *= temp;
    }
    
    return result;
}


template<unsigned int n>
constexpr static double ipow(const double x){
    
    double result{1}, temp{x};
    
    for(unsigned int j=1; j<=n; j<<=1){
       if(n&j){
           result *= temp;
       }
       temp *= temp;
    }
    
    return result;
}

#endif
