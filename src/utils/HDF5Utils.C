#include "HDF5Utils.h"
#include "MooseError.h"
#include "mpi.h"
#include <string>

namespace hdf5_utils {

bool attribute_exists(hid_t obj_id, const char *name) {

  htri_t out = H5Aexists_by_name(obj_id, ".", name, H5P_DEFAULT);
  // If the attribute exists out will be 1
  return out > 0;
}

bool object_exists(hid_t object_id, const char *name) {
  htri_t out = H5LTpath_valid(object_id, name, true);
  if (out < 0) {
    std::string error =
        "Failed to check if object " + std::string(name) + " exists.";
    mooseError(error);
  }
  return (out > 0);
}

void get_shape(hid_t obj_id, hsize_t *dims) {

  auto type = H5Iget_type(obj_id);
  hid_t dspace;
  if (type == H5I_DATASET) {
    dspace = H5Dget_space(obj_id);
  } else if (type == H5I_ATTR) {
    dspace = H5Aget_space(obj_id);
  } else {
    throw std::runtime_error{
        "Expected dataset or attribute in call to get_shape."};
  }
  H5Sget_simple_extent_dims(dspace, dims, nullptr);
  H5Sclose(dspace);
}

std::string object_name(hid_t obj_id) {
  // Determine size and create buffer
  size_t size = 1 + H5Iget_name(obj_id, nullptr, 0);
  char *buffer = new char[size];

  // Read and return name
  H5Iget_name(obj_id, buffer, size);
  std::string str = buffer;
  delete[] buffer;
  return str;
}

void ensure_exists(hid_t obj_id, const char *name, bool attribute) {
  if (attribute) {
    if (!attribute_exists(obj_id, name)) {

      std::string error_msg = "Attribute " + std::string(name) +
                              " does not exist in object " +
                              object_name(obj_id);

      mooseError(error_msg);
    }
  } else {
    if (!object_exists(obj_id, name)) {

      std::string error_msg = "Object " + std::string(name) +
                              " does not exist in object " +
                              std::string(object_name(obj_id));

      mooseError(error_msg);
    }
  }
}

hid_t file_open(const char *filename, char mode, bool parallel,
                const MPI_Comm &comm) {
  bool create;
  unsigned int flags;
  switch (mode) {
  case 'r':
  case 'a':
    create = false;
    flags = (mode == 'r' ? H5F_ACC_RDONLY : H5F_ACC_RDWR);
    break;
  case 'w':
  case 'x':
    create = true;
    flags = (mode == 'x' ? H5F_ACC_EXCL : H5F_ACC_TRUNC);
    break;
  default:
    std::string error_msg = "Invalid file mode: " + std::to_string(mode);
    mooseError(error_msg);
  }

  hid_t plist = H5P_DEFAULT;
#ifdef H5_HAVE_PARALLEL
  if (parallel) {
    // Setup file access property list with parallel I/O access
    plist = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_all_coll_metadata_ops(plist, true);
    H5Pset_fapl_mpio(plist, comm, MPI_INFO_NULL);
  }
#endif

  // Open the file collectively
  hid_t file_id;
  if (create) {
    file_id = H5Fcreate(filename, flags, H5P_DEFAULT, plist);
  } else {
    file_id = H5Fopen(filename, flags, plist);
  }
  if (file_id < 0) {
    std::string error_msg = "Failed to open HDF5 file with mode " +
                            std::to_string(mode) + ": " + filename;
    mooseError(error_msg);
  }

#ifdef H5_HAVE_PARALLEL
  // Close the property list
  if (parallel)
    H5Pclose(plist);
#endif

  return file_id;
}

void file_close(hid_t file_id) { H5Fclose(file_id); }

hid_t open_group(hid_t group_id, const char *name) {
  ensure_exists(group_id, name);
  return H5Gopen(group_id, name, H5P_DEFAULT);
}

hid_t open_group(hid_t group_id, const std::string &name) {
  return open_group(group_id, name.c_str());
}

void close_group(hid_t group_id) {
  if (H5Gclose(group_id) < 0)
    mooseError("Failed to close group");
}

hid_t open_dataset(hid_t group_id, const char *name) {
  ensure_exists(group_id, name);
  return H5Dopen(group_id, name, H5P_DEFAULT);
}

void close_dataset(hid_t dataset_id) {
  if (H5Dclose(dataset_id) < 0)
    mooseError("Failed to close dataset");
}

bool using_mpio_device(hid_t obj_id) {

  // Determine file that this object is part of
  hid_t file_id = H5Iget_file_id(obj_id);

  // Get file access property list
  hid_t fapl_id = H5Fget_access_plist(file_id);

  // Get low-level driver identifier
  hid_t driver = H5Pget_driver(fapl_id);

  // Free resources
  H5Pclose(fapl_id);
  H5Fclose(file_id);

  return driver == H5FD_MPIO;
}

void write_dataset_lowlevel(hid_t obj_id, const char *name, int ndim,
                            const hsize_t *dims, hid_t mem_type_id,
                            void *buffer, bool parallel, bool indep,
                            hid_t file_space_id) {

  hid_t dataset = obj_id;
  if (name) {
    dataset = open_dataset(obj_id, name);
  }
  // If array is given, create a simple dataspace. Otherwise, create a
  // scalar datascape.
  hid_t mem_space = H5Screate_simple(ndim, dims, NULL);

  if (parallel) {
#ifdef H5_HAVE_PARALLEL
    // Set up collective vs independent I/O
    auto data_xfer_mode = indep ? H5FD_MPIO_INDEPENDENT : H5FD_MPIO_COLLECTIVE;

    // Create dataset transfer property list
    hid_t plist = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist, data_xfer_mode);

    // Write data
    H5Dwrite(dataset, mem_type_id, mem_space, file_space_id, plist, buffer);
    H5Pclose(plist);
#endif
  } else {
    H5Dwrite(dataset, mem_type_id, mem_space, file_space_id, H5P_DEFAULT,
             buffer);
  }

  // Free resources
  H5Sclose(mem_space);
  if (name) {
    H5Dclose(dataset);
  }
}

void read_dataset_lowlevel(hid_t obj_id, const char *name, hid_t mem_type_id,
                           hid_t mem_space_id, void *buffer, bool parallel,
                           bool indep, hid_t file_space_id) {

  hid_t dataset = obj_id;
  if (name) {
    dataset = open_dataset(obj_id, name);
  }

  if (parallel) {
#ifdef H5_HAVE_PARALLEL
    // Set up collective vs independent I/O
    auto data_xfer_mode = indep ? H5FD_MPIO_INDEPENDENT : H5FD_MPIO_COLLECTIVE;

    // Create dataset transfer property list
    hid_t plist = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist, data_xfer_mode);

    // Read data
    H5Dread(dataset, mem_type_id, mem_space_id, file_space_id, plist, buffer);
    H5Pclose(plist);
#endif
  } else {
    H5Dread(dataset, mem_type_id, mem_space_id, file_space_id, H5P_DEFAULT,
            buffer);
  }

  if (name) {
    H5Dclose(dataset);
  }
}

void read_double(hid_t obj_id, const char *name, double *buffer, bool parallel,
                 bool indep) {
  read_dataset_lowlevel(obj_id, name, H5T_NATIVE_DOUBLE, H5S_ALL, buffer,
                        parallel, indep);
}

void read_int(hid_t obj_id, const char *name, int *buffer, bool parallel,
              bool indep) {
  read_dataset_lowlevel(obj_id, name, H5T_NATIVE_INT, H5S_ALL, buffer, parallel,
                        indep);
}

void read_string(hid_t obj_id, const char *name,
                 std::vector<std::string> &buffer, int slen, bool parallel,
                 bool indep) {

  /// Currently easier to read into a char buffer and then use the predefined
  /// string length to seperate the char buffer into strings
  std::vector<char> char_buff(slen * buffer.size());

  /// Create datatype for "string"
  hid_t datatype = H5Tcopy(H5T_C_S1);
  H5Tset_size(datatype, slen);

  H5Tset_strpad(datatype, H5T_STR_NULLPAD);

  /// Read strings into char buffer
  read_dataset_lowlevel(obj_id, name, datatype, H5S_ALL, char_buff.data(),
                        parallel, indep);

  /// Split vector of chars into strings of length slen
  for (size_t i = 0; i < buffer.size(); i++) {
    std::string s(char_buff.data() + i * slen, slen);

    auto pos = s.find_last_not_of(std::string("\0 ", 2));
    if (pos != std::string::npos)
      s.erase(pos + 1);
    else
      s.clear();
    buffer[i] = std::move(s);
  }

  /// Free memory allocated to datatype
  H5Tclose(datatype);
}

void read_double_hyperslab(hid_t obj_id, const char *name, hsize_t ndim,
                           hsize_t dims[], hsize_t start[], hsize_t count[],
                           double *results, bool parallel, bool indep) {
  hid_t dataset = obj_id;
  if (name) {
    dataset = open_dataset(obj_id, name);
  }
  hid_t mem_space_id = H5Screate_simple(ndim, dims, nullptr);

  hid_t file_space_id = H5Dget_space(dataset);

  H5Sselect_hyperslab(file_space_id, H5S_SELECT_SET, start, nullptr, count,
                      nullptr);

  read_dataset_lowlevel(obj_id, name, H5T_NATIVE_DOUBLE, mem_space_id, results,
                        parallel, indep, file_space_id);

  H5Sclose(mem_space_id);
  H5Sclose(file_space_id);

  if (name) {
    H5Dclose(dataset);
  }
}

void write_double_hyperslab(hid_t obj_id, const char *name, hsize_t ndim,
                            hsize_t dims[], hsize_t start[], hsize_t count[],
                            double *results, bool parallel, bool indep) {

  hid_t dataset = obj_id;
  if (name) {
    dataset = open_dataset(obj_id, name);
  }
  hid_t file_space_id = H5Dget_space(dataset);

  H5Sselect_hyperslab(file_space_id, H5S_SELECT_SET, start, nullptr, count,
                      nullptr);

  write_dataset_lowlevel(obj_id, name, ndim, dims, H5T_NATIVE_DOUBLE, results,
                         parallel, indep, file_space_id);

  H5Sclose(file_space_id);
  if (name) {
    H5Dclose(dataset);
  }
}

} // namespace hdf5_utils
