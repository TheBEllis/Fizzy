#include <cmath>
#include <unordered_map>
#include <vector>

namespace utils {
namespace energy_groups {

inline std::vector<double> ukaea709() {
  const int groups = 709;
  const double Emin = 1e-5;
  const double Emax = 1e9;

  std::vector<double> bounds(groups + 1);

  double logEmin = std::log10(Emin);
  double logEmax = std::log10(Emax);

  double dlog = (logEmax - logEmin) / groups;

  for (int i = 0; i <= groups; ++i)
    bounds[i] = std::pow(10.0, logEmin + i * dlog);

  return bounds;
}

inline std::vector<double> ukaea1102() {
  const int groups = 1102;
  const double Emin = 1e-5;
  const double Emax = 1e9;

  std::vector<double> bounds(groups + 1);

  double logEmin = std::log10(Emin);
  double logEmax = std::log10(Emax);

  double dlog = (logEmax - logEmin) / groups;

  for (int i = 0; i <= groups; ++i)
    bounds[i] = std::pow(10.0, logEmin + i * dlog);

  return bounds;
}

inline std::unordered_map<size_t, std::vector<double>> gamma_groups = {
    {24, {1.000e-11, 1.000e+4, 2.000e+4, 5.000e+4, 1.000e+5, 2.000e+5, 3.000e+5,
          4.000e+5,  6.000e+5, 8.000e+5, 1.000e+6, 1.220e+6, 1.440e+6, 1.660e+6,
          2.000e+6,  2.500e+6, 3.000e+6, 4.000e+6, 5.000e+6, 6.500e+6, 8.000e+6,
          1.000e+7,  1.200e+7, 1.400e+7, 2.000e+7}},
    {22, {0.0,   1.0e4, 1.0e5, 2.0e5, 4.0e5, 1.0e6, 1.5e6, 2.0e6,
          2.5e6, 3.0e6, 3.5e6, 4.0e6, 4.5e6, 5.0e6, 5.5e6, 6.0e6,
          6.5e6, 7.0e6, 7.5e6, 8.0e6, 1.0e7, 1.2e7, 1.4e7}}};

inline std::unordered_map<size_t, std::vector<double>> neutron_groups = {
    {1102, ukaea1102()}, {709, ukaea709()}};

}; // namespace energy_groups
}; // namespace utils
