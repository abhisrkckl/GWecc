#ifndef _FourierWaveform_hpp_
#define _FourierWaveform_hpp_ 1

#include <Eigen/Dense>
#include "binary.hpp"
#include "FourierSize.incl"
#include <tuple>

struct Fourier2DCoeffs_t{
	Array_pq Ap, Bp, Cp, Dp,
		     Ax, Bx, Cx, Dx;
};

//Fourier2DCoeffs_t FourierWaveformCoeffs(const BinaryMass &binmass, const BinaryState &binstate, const double beta);
Fourier2DCoeffs_t FourierResidualCoeffs(const BinaryMass &binmass, const BinaryState &binstate, const double beta);

std::tuple<double,double> FourierResidual_pt(const Fourier2DCoeffs_t &res_coeffs, const double l, const double lambda);

#endif
