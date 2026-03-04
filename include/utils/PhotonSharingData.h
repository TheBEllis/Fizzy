#ifndef PHOTON_SHARING_DATA_H
#define PHOTON_SHARING_DATA_H

#include "boost/interprocess/containers/map.hpp"
#include "boost/interprocess/containers/vector.hpp"
#include "boost/interprocess/managed_shared_memory.hpp"
#include "boost/interprocess/shared_memory_object.hpp"

#include <boost/interprocess/interprocess_fwd.hpp>
#include <cstdint>
#include <map>
#include <string>
#include <sys/types.h>
#include <unordered_map>
#include <vector>

namespace bi = boost::interprocess;

/// Allocator for vector of doubles
typedef boost::interprocess::allocator<
    double, boost::interprocess::managed_shared_memory::segment_manager>
    ShmemVecAllocator;

typedef boost::interprocess::vector<double, ShmemVecAllocator> BoostIpVector;

typedef boost::interprocess::allocator<
    std::pair<const uint64_t, uint64_t>,
    bi::managed_shared_memory::segment_manager>
    ShmemIntIntMapAllocator;

typedef bi::map<uint64_t, uint64_t, std::less<uint64_t>,
                ShmemIntIntMapAllocator>
    BoostIpIntIntMap;

typedef boost::interprocess::allocator<
    std::pair<const int, double>, bi::managed_shared_memory::segment_manager>
    ShmemIntDoubMapAllocator;

typedef bi::map<int, double, std::less<int>, ShmemIntDoubMapAllocator>
    BoostIpIntDoubMap;

class PhotonSharingData {
public:
  // PhotonSharingData(
  //     bi::managed_shared_memory &segment, std::vector<double> &photon_flux,
  //     std::vector<double> &elem_strengths, uint64_t num_photon_bins,
  //     uint64_t num_local_elems, double total_domain_strength,
  //     double local_domain_strength, std::vector<double> &photon_bins,
  //     std::unordered_map<uint64_t, uint64_t> &local_element_idx_map)
  //     : _num_photon_bins(num_photon_bins), _num_local_elems(num_local_elems),
  //       _local_domain_strength(local_domain_strength),
  //       _total_domain_strength(total_domain_strength), _is_setup(false),
  //       _vec_alloc(segment.get_segment_manager()),
  //       _int_int_map_alloc(segment.get_segment_manager()),
  //       _int_doub_map_alloc(segment.get_segment_manager()),
  //       _photon_fluxes(_vec_alloc), _elem_strength(_vec_alloc),
  //       _photon_bins(_vec_alloc), _local_elem_idx_map(_int_int_map_alloc) {
  //
  //   // Initialise shared memory _photon_flux with values from
  //   // FispactProblem's _photon_fluxes
  //   _photon_fluxes =
  //       BoostIpVector(photon_flux.begin(), photon_flux.end(), _vec_alloc);
  //
  //   _elem_strength =
  //       BoostIpVector(elem_strengths.begin(), elem_strengths.end(),
  //       _vec_alloc);
  //
  //   _photon_bins =
  //       BoostIpVector(photon_bins.begin(), photon_bins.end(), _vec_alloc);
  //
  //   for (auto &[global_elem_id, local_elem_id] : local_element_idx_map) {
  //
  //     _local_elem_idx_map.insert(
  //         std::pair<uint64_t, uint64_t>(global_elem_id, local_elem_id));
  //   }
  // }

  PhotonSharingData(bi::managed_shared_memory &segment)
      : _is_setup(false), _vec_alloc(segment.get_segment_manager()),
        _int_int_map_alloc(segment.get_segment_manager()),
        _int_doub_map_alloc(segment.get_segment_manager()),
        _photon_fluxes(_vec_alloc), _elem_strength(_vec_alloc),
        _photon_bins(_vec_alloc), _local_elem_idx_map(_int_int_map_alloc) {}

  void setPhotonSpectra(const std::vector<double>::iterator spectra_begin,
                        const std::vector<double>::iterator spectra_end) {

    _photon_fluxes.assign(spectra_begin, spectra_end);
  }

  void
  setElementStrengths(const std::vector<double>::iterator elem_strength_begin,
                      const std::vector<double>::iterator elem_strength_end) {
    _elem_strength.assign(elem_strength_begin, elem_strength_end);
  }

  void setPhotonBins(const std::vector<double> &photon_bins) {
    _photon_bins.assign(photon_bins.begin(), photon_bins.end());
  }

  void setNumPhotonBins(const size_t &n_photon_bins) {
    _num_photon_bins = n_photon_bins;
  }

  void setNumLocalElems(const size_t &n_local_elems) {
    _num_local_elems = n_local_elems;
  }

  void setLocalDomainStrength(const double &local_domain_strength) {
    _local_domain_strength = local_domain_strength;
  }

  void setTotalDomainStrength(const double &total_domain_strength) {
    _total_domain_strength = total_domain_strength;
  }

  void setLocalElemIdMap(
      std::unordered_map<uint64_t, uint64_t> &local_element_idx_map) {

    for (auto &[global_elem_id, local_elem_id] : local_element_idx_map) {

      _local_elem_idx_map.insert(
          std::pair<uint64_t, uint64_t>(global_elem_id, local_elem_id));
    }
  }

  size_t _num_photon_bins;
  size_t _num_local_elems;
  double _local_domain_strength;
  double _total_domain_strength;

  bool _is_setup;

  /// Allocators
  ShmemVecAllocator _vec_alloc;
  ShmemIntIntMapAllocator _int_int_map_alloc;
  ShmemIntDoubMapAllocator _int_doub_map_alloc;

  /// Data structures
  BoostIpVector _photon_fluxes;
  BoostIpVector _elem_strength;
  BoostIpVector _photon_bins;

  BoostIpIntIntMap _local_elem_idx_map;
};

#endif
