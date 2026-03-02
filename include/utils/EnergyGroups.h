#include "fispactnucleardata.hpp"
#include <unordered_map>

namespace utils {
namespace energy_groups {

inline std::unordered_map<size_t, std::vector<double>> gamma_groups = {
    {24, {1.000e-11, 1.000e+4, 2.000e+4, 5.000e+4, 1.000e+5, 2.000e+5, 3.000e+5,
          4.000e+5,  6.000e+5, 8.000e+5, 1.000e+6, 1.220e+6, 1.440e+6, 1.660e+6,
          2.000e+6,  2.500e+6, 3.000e+6, 4.000e+6, 5.000e+6, 6.500e+6, 8.000e+6,
          1.000e+7,  1.200e+7, 1.400e+7, 2.000e+7}},
    {22, {0.0,   1.0e4, 1.0e5, 2.0e5, 4.0e5, 1.0e6, 1.5e6, 2.0e6,
          2.5e6, 3.0e6, 3.5e6, 4.0e6, 4.5e6, 5.0e6, 5.5e6, 6.0e6,
          6.5e6, 7.0e6, 7.5e6, 8.0e6, 1.0e7, 1.2e7, 1.4e7}}};

}; // namespace energy_groups
}; // namespace utils
