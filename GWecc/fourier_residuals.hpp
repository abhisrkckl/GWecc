#ifndef _FourierWaveform_hpp_
#define _FourierWaveform_hpp_ 1

#include <tuple>
#include <Eigen/Dense>
#include "binary.hpp"

constexpr size_t qmax = 4;
constexpr size_t pmax = 8;
typedef Eigen::Array<double,pmax+1,qmax+1> Array_pq;

struct Fourier2DCoeffs{
	Array_pq Ap, Bp, Cp, Dp,
		     Ax, Bx, Cx, Dx;
};

//Fourier2DCoeffs_t FourierWaveformCoeffs(const BinaryMass &binmass, const BinaryState &binstate, const double beta);
Fourier2DCoeffs fourier_residual_coeffs(const BinaryMass &binmass, const BinaryState &binstate, const double beta);

std::tuple<double,double> fourier_residual_pt(const Fourier2DCoeffs &res_coeffs, const double l, const double lambda);

#endif
