#include <cassert>
#include <vector>

class PhotonSpectra {
public:
  PhotonSpectra(size_t nt, size_t ne, size_t nb)
      : _nt(nt), _ne(ne), _nb(nb), _data(_nt * _ne * _nb, 0) {}

  /// Accessors
  // Scalar data accessor
  double &at(size_t t, size_t e, size_t b) { return _data.at(idx(t, e, b)); }

  const double &at(size_t t, size_t e, size_t b) const {
    return _data.at(idx(t, e, b));
  }

  // Single element spectrum at time
  std::vector<double>::iterator spectrum_begin(size_t t, size_t e) {
    return _data.begin() + idx(t, e, 0);
  }

  std::vector<double>::iterator spectrum_end(size_t t, size_t e) {
    return spectrum_begin(t, e) + _nb;
  }

  // Entire spectrum at time
  std::vector<double>::iterator time_begin(size_t t) {
    return _data.begin() + idx(t, 0, 0);
  }
  std::vector<double>::iterator time_end(size_t t) {
    return time_begin(t) + (_ne * _nb);
  }

  // Data accessors for debugging
  std::vector<double> &data() { return _data; };
  double *data_buffer() { return _data.data(); }

private:
  size_t idx(size_t t, size_t e, size_t b) const {

    // dbg checks
    assert(t < _nt);
    assert(e < _ne);
    assert(b < _nb);

    return (t * _nb * _ne) + (e * _nb) + b;
  }

  /// Number of timesteps, number of elements, number of bins
  size_t _nt, _ne, _nb;
  std::vector<double> _data;
};
