###############################################################################
################### MOOSE Application Standard Makefile #######################
###############################################################################
#
# Optional Environment variables
# MOOSE_DIR        - Root directory of the MOOSE project
#
###############################################################################
# Use the MOOSE submodule if it exists and MOOSE_DIR is not set
MOOSE_SUBMODULE    := $(CURDIR)/moose
ifneq ($(wildcard $(MOOSE_SUBMODULE)/framework/Makefile),)
  MOOSE_DIR        ?= $(MOOSE_SUBMODULE)
else
  MOOSE_DIR        ?= $(shell dirname `pwd`)/moose
endif

# framework
FRAMEWORK_DIR      := $(MOOSE_DIR)/framework
include $(FRAMEWORK_DIR)/build.mk
include $(FRAMEWORK_DIR)/moose.mk

################################## MODULES ####################################
# To use certain physics included with MOOSE, set variables below to
# yes as needed.  Or set ALL_MODULES to yes to turn on everything (overrides
# other set variables).

ALL_MODULES                 := no

CHEMICAL_REACTIONS          := no
CONTACT                     := no
ELECTROMAGNETICS            := no
EXTERNAL_PETSC_SOLVER       := no
FLUID_PROPERTIES            := no
FSI                         := no
FUNCTIONAL_EXPANSION_TOOLS  := no
GEOCHEMISTRY                := no
HEAT_TRANSFER               := no
LEVEL_SET                   := no
MISC                        := no
NAVIER_STOKES               := no
OPTIMIZATION                := no
PERIDYNAMICS                := no
PHASE_FIELD                 := no
POROUS_FLOW                 := no
RAY_TRACING                 := no
REACTOR                     := no
RDG                         := no
RICHARDS                    := no
SOLID_MECHANICS             := no
STOCHASTIC_TOOLS            := no
THERMAL_HYDRAULICS          := no
XFEM                        := no

include $(MOOSE_DIR)/modules/modules.mk
###############################################################################

# dep apps
APPLICATION_DIR    := $(CURDIR)
APPLICATION_NAME   := Fizzy
BUILD_EXEC         := yes
GEN_REVISION       := no
include            $(FRAMEWORK_DIR)/app.mk

# ======================================================================================
# PETSc
# ======================================================================================
PETSC_DIR           ?= $(MOOSE_DIR)/petsc
PETSC_ARCH          ?= arch-moose
LIBMESH_DIR         ?= $(MOOSE_DIR)/libmesh/installed/
# Use compiler info discovered by PETSC
ifeq ($(PETSC_ARCH),)
	include $(PETSC_DIR)/$(PETSC_ARCH)/lib/petsc/conf/petscvariables
else
	include $(PETSC_DIR)/lib/petsc/conf/petscvariables
endif

# libmesh_CXX, etc, were defined in build.mk
export CXX := $(libmesh_CXX)
export CC  := $(libmesh_CC)
export FC  := $(libmesh_F90)
export FFLAGS := $(libmesh_FFLAGS)
export CFLAGS := $(libmesh_CFLAGS)
export CXXFLAGS := $(libmesh_CXXFLAGS)
export CPPFLAGS := $(libmesh_CPPFLAGS)
export LDFLAGS := $(libmesh_LDFLAGS)
export LIBS := $(libmesh_LIBS)

FISPACT_DIR ?= ${MOOSE_DIR}/../FISPACT/ubuntu/20.10
FISPACT_DIR ?= ${MOOSE_DIR}
FISPACT_INCLUDES ?= -I ${FISPACT_DIR}/include/c -I ${FISPACT_DIR}/include/cpp
FISPACT_LIB_DIR ?= ${FISPACT_DIR}/lib/linux
# FISPACT_LIB ?= 

PUGIXML_DIR ?= $(MOOSE_DIR)/../pugixml/
PUGIXML_INCLUDES ?= -I $(PUGIXML_DIR)/src/
PUGIXML_LIB_DIR ?= ${PUGIXML_DIR}/build/

HDF5_DIR ?= /usr/lib/x86_64-linux-gnu/hdf5/openmpi/lib/

ADDITIONAL_LIBS := -L$(FISPACT_LIB_DIR) -ljsonfortran -lmonitor -lfispact -lfispactapi -L$(PUGIXML_LIB_DIR) -lpugixml -L$(HDF5_DIR) -lhdf5
ADDITIONAL_LIBS += $(CC_LINKER_SLFLAG)$(FISPACT_LIB_DIR)

ADDITIONAL_CPPFLAGS += $(FISPACT_INCLUDES) ${PUGIXML_INCLUDES}

FIZZY_EXTERNAL_FLAGS = ${ADDITIONAL_LIBS} ${ADDITIONAL_CPPFLAGS}

$(app_LIB): EXTERNAL_FLAGS := $(FIZZY_EXTERNAL_FLAGS)
$(app_test_LIB): EXTERNAL_FLAGS := $(FIZZY_EXTERNAL_FLAGS)
$(app_EXEC): EXTERNAL_FLAGS := $(FIZZY_EXTERNAL_FLAGS)
