#pragma once

#include "H5Cpp.h"
#include <algorithm>
#include <cstddef>
#include <cstring>
#include <hdf5/openmpi/H5Apublic.h>
#include <hdf5/openmpi/H5Dpublic.h>
#include <hdf5/openmpi/H5FDmpio.h>
#include <hdf5/openmpi/H5Fpublic.h>
#include <hdf5/openmpi/H5Ipublic.h>
#include <hdf5/openmpi/H5Ppublic.h>
#include <hdf5/openmpi/H5Spublic.h>
#include <hdf5/openmpi/H5Tpublic.h>
#include <hdf5/openmpi/H5public.h>
#include <hdf5/openmpi/H5version.h>
#include <sstream>
#include <string>
#include <vector>

#include "MooseError.h"

#include "hdf5.h"
#include "hdf5_hl.h"

namespace hdf5_utils {

/*
 *Check an attribute exists for a given object, using the object id and the
 *name of the attribute
 *
 *
 **/
bool attribute_exists(hid_t obj_id, const char *name);
/***/
bool object_exists(hid_t object_id, const char *name);

/*
 * Get the shape of an object
 *
 *
 **/
void get_shape(hid_t obj_id, hsize_t *dims);

std::string object_name(hid_t obj_id);
/***/
void ensure_exists(hid_t obj_id, const char *name, bool attribute = false);

/***/
hid_t file_open(const char *filename, char mode, bool parallel,
                const MPI_Comm &comm);

/***/
void file_close(hid_t file_id);

/***/

hid_t open_dataset(hid_t group_id, const char *name);

void close_dataset(hid_t dataset_id);

hid_t open_group(hid_t group_id, const char *name);

hid_t open_group(hid_t group_id, const std::string &name);

void close_group(hid_t group_id);

void read_dataset_lowlevel(hid_t obj_id, const char *name, hid_t mem_type_id,
                           hid_t mem_space_id, void *buffer, bool parallel,
                           bool indep = true, hid_t file_space_id = H5S_ALL);

void write_dataset_lowlevel(hid_t obj_id, const char *name, int ndim,
                            const hsize_t *dims, hid_t mem_type_id,
                            void *buffer, bool parallel, bool indep = true,
                            hid_t file_space_id = H5S_ALL);

void read_double(hid_t obj_id, const char *name, double *buffer, bool parallel,
                 bool indep = true);

void read_int(hid_t obj_id, const char *name, int *buffer, bool parallel,
              bool indep = true);

bool using_mpio_device(hid_t obj_id);

void read_double_hyperslab(hid_t group_id, const char *name, hsize_t rank,
                           hsize_t dims[], hsize_t start[], hsize_t count[],
                           double *results, bool parallel = true,
                           bool indep = true);

void write_double_hyperslab(hid_t obj_id, const char *name, hsize_t ndim,
                            hsize_t dims[3], hsize_t start[], hsize_t count[],
                            double *results, bool parallel = true,
                            bool indep = true);
} // namespace hdf5_utils
