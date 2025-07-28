#ifndef PHOTON_SHARING_DATA_H
#define PHOTON_SHARING_DATA_H

#include "boost/interprocess/containers/map.hpp"
#include "boost/interprocess/containers/vector.hpp"
#include "boost/interprocess/managed_shared_memory.hpp"
#include "boost/interprocess/shared_memory_object.hpp"

#include <boost/interprocess/interprocess_fwd.hpp>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

namespace bi = boost::interprocess;
typedef boost::interprocess::allocator<
    double, boost::interprocess::managed_shared_memory::segment_manager>
    ShmemVecAllocator;

typedef boost::interprocess::vector<double, ShmemVecAllocator> BoostIpVector;

typedef boost::interprocess::allocator<
    std::pair<const int, BoostIpVector>,
    bi::managed_shared_memory::segment_manager>
    ShmemIntVecMapAllocator;

typedef bi::map<int, BoostIpVector, std::less<int>, ShmemIntVecMapAllocator>
    BoostIpIntVecMap;

typedef boost::interprocess::allocator<
    std::pair<const int, double>, bi::managed_shared_memory::segment_manager>
    ShmemIntDoubMapAllocator;

typedef bi::map<int, double, std::less<int>, ShmemIntDoubMapAllocator>
    BoostIpIntDoubMap;

class PhotonSharingData {
public:
  PhotonSharingData(bi::managed_shared_memory &segment,
                    std::unordered_map<int, std::vector<double>> &photon_flux,
                    std::unordered_map<int, double> &elem_strength,
                    int num_energy_bins, int num_local_elems,
                    double total_domain_strength, double local_domain_strength)
      : _num_energy_bins(num_energy_bins), _num_local_elems(num_local_elems),
        _local_domain_strength(local_domain_strength),
        _total_domain_strength(total_domain_strength), _communication(true),
        _vec_alloc(segment.get_segment_manager()),
        _int_vec_map_alloc(segment.get_segment_manager()),
        _int_doub_map_alloc(segment.get_segment_manager()),
        _photon_fluxes(_int_vec_map_alloc),
        _elem_strength(_int_doub_map_alloc) {

    // Initialise shared memory _photon_flux map with values from
    // FispactProblem's _photon_fluxes
    for (auto &pair : photon_flux) {

      _photon_fluxes.insert(std::pair<int, BoostIpVector>(
          pair.first,
          BoostIpVector(pair.second.begin(), pair.second.end(), _vec_alloc)));
    }

    // Do the same for elem_strength map
    for (auto &pair : elem_strength) {

      _elem_strength.insert(std::pair<int, double>(pair.first, pair.second));
    }
  }

  int _num_energy_bins;
  int _num_local_elems;
  double _local_domain_strength;
  double _total_domain_strength;

  bool _communication;

  // Allocators
  ShmemVecAllocator _vec_alloc;
  ShmemIntVecMapAllocator _int_vec_map_alloc;
  ShmemIntDoubMapAllocator _int_doub_map_alloc;

  BoostIpIntVecMap _photon_fluxes;
  BoostIpIntDoubMap _elem_strength;
};

#endif
