#ifndef _mikkola_h_
#define _mikkola_h_

#include <vector>

double MIKKOLA(double l, const double e);
double MIKKOLAh(double l, const double e);

std::vector<double> MIKKOLA(const std::vector<double> &ls, const double e);
std::vector<double> MIKKOLAh(const std::vector<double> &ls, const double e);

std::vector<double> MIKKOLA(const std::vector<double> &ls, const std::vector<double> &es);
std::vector<double> MIKKOLAh(const std::vector<double> &ls, const std::vector<double> &es);

#endif
