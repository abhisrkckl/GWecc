#ifndef _NumericalWaveform_hpp_
#define _NumericalWaveform_hpp_

#include "EccentricResiduals.hpp"
#include "OrbitalEvolution.hpp"
#include <gsl/gsl_integration.h>
#include <tuple>

struct WaveformParams{
	const double Fp, Fx, DGW;
	const BinaryMass &bin_mass;
	const BinaryState &bin_init;
	const EvolveCoeffs_t &ev_coeffs;
};

class GSL_QAG_Integrator{
	
private:
	
	const size_t workspace_size;
	gsl_integration_workspace *workspace;
	
	int integ_type;
	
	double epsabs, epsrel;
	
public:
	GSL_QAG_Integrator(const size_t _workspace_size, const double _epsabs, const double _epsrel, const int key)
		: workspace_size(_workspace_size),
		  integ_type(key),
		  epsabs(_epsabs),
		  epsrel(_epsrel) {
		  
		workspace = gsl_integration_workspace_alloc(workspace_size);
	}
	
	~GSL_QAG_Integrator(){
		gsl_integration_workspace_free(workspace);
	}
	
	std::tuple<double,double> eval(gsl_function &func, const double x0, const double x1) const{
		double result, error;
		
		gsl_integration_qag(&func, x0, x1, epsabs, epsrel, workspace_size, integ_type, workspace, &result, &error);
		
		return std::make_tuple(result,error);
	}
	
	std::tuple<Signal1D,Signal1D> eval(gsl_function &func, const double x0, const Signal1D &x) const {
		
		const size_t length = x.size();
		Signal1D results(length);
		Signal1D errors(length);
		
		double x_prev = x0;
		double res_prev = 0;
		
		for(unsigned int i=0; i<length; i++){
			double result, error;
			gsl_integration_qag(&func, x_prev, x[i], epsabs, epsrel, workspace_size, integ_type, workspace, &result, &error);
			
			res_prev += result;
			results[i] = res_prev;
			errors[i] = error;
			
			x_prev = x[i];
		}
		
		return std::make_tuple(results,errors);
	}
	
	Signal1D eval_noerr(gsl_function &func, const double x0, const Signal1D &x) const {
		
		const size_t length = x.size();
		Signal1D results(length);
		
		double x_prev = x0;
		double res_prev = 0;
		
		for(unsigned int i=0; i<length; i++){
			double result, error;
			gsl_integration_qag(&func, x_prev, x[i], epsabs, epsrel, workspace_size, integ_type, workspace, &result, &error);
			
			res_prev += result;
			results[i] = res_prev;
			
			x_prev = x[i];
		}
		
		return results;
	}
};


#endif
