#ifndef _mikkola_h_
#define _mikkola_h_

#include <vector>

double mikkola(double l, const double e);

std::vector<double> mikkola(const std::vector<double> &ls, const double e);

std::vector<double> mikkola(const std::vector<double> &ls, const std::vector<double> &es);

#endif
