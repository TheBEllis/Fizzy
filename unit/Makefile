###############################################################################
################### MOOSE Application Standard Makefile #######################
###############################################################################
#
# Required Environment variables (one of the following)
# PACKAGES_DIR  - Location of the MOOSE redistributable package
#
# Optional Environment variables
# MOOSE_DIR     - Root directory of the MOOSE project
# FRAMEWORK_DIR - Location of the MOOSE framework
#
###############################################################################
# Use the MOOSE submodule if it exists and MOOSE_DIR is not set
MOOSE_SUBMODULE    := $(CURDIR)/../moose
ifneq ($(wildcard $(MOOSE_SUBMODULE)/framework/Makefile),)
  MOOSE_DIR        ?= $(MOOSE_SUBMODULE)
else
  MOOSE_DIR        ?= $(shell dirname `pwd`)/../moose
endif
FRAMEWORK_DIR      ?= $(MOOSE_DIR)/framework
###############################################################################

# framework
include $(FRAMEWORK_DIR)/build.mk
include $(FRAMEWORK_DIR)/moose.mk

################################## MODULES ####################################
# set desired physics modules equal to 'yes' to enable them
CHEMICAL_REACTIONS        := no
CONTACT                   := no
FLUID_PROPERTIES          := no
FSI                       := no
HEAT_TRANSFER             := no
MISC                      := no
NAVIER_STOKES             := no
PHASE_FIELD               := no
RDG                       := no
RICHARDS                  := no
STOCHASTIC_TOOLS          := no
TENSOR_MECHANICS          := no
XFEM                      := no
POROUS_FLOW               := no
LEVEL_SET                 := no
include           $(MOOSE_DIR)/modules/modules.mk
###############################################################################

# Extra stuff for GTEST
ADDITIONAL_INCLUDES := -I$(FRAMEWORK_DIR)/contrib/gtest
ADDITIONAL_LIBS     := $(FRAMEWORK_DIR)/contrib/gtest/libgtest.la

PUGIXML_DIR ?= $(MOOSE_DIR)/../pugixml/
PUGIXML_INCLUDES ?= -I $(PUGIXML_DIR)/src/
PUGIXML_LIB_DIR ?= ${PUGIXML_DIR}/build/

HDF5_DIR ?= /usr/lib/x86_64-linux-gnu/hdf5/openmpi/lib/

ADDITIONAL_LIBS += -L$(PUGIXML_LIB_DIR) -lpugixml -L$(HDF5_DIR) -lhdf5 
ADDITIONAL_CPPFLAGS += $(FISPACT_INCLUDES) ${PUGIXML_INCLUDES} -DFIZZY_UNIT_TEST
# dep apps
CURRENT_DIR        := $(shell pwd)
APPLICATION_DIR    := $(CURRENT_DIR)/..
APPLICATION_NAME   := Fizzy

FISPACT_STUB_DIR   := $(CURRENT_DIR)/stubs
FISPACT_STUB_LIB   := $(FISPACT_STUB_DIR)/libfispact_stub.a

$(FISPACT_STUB_LIB): $(FISPACT_STUB_DIR)/fispact_stub.C
	$(libmesh_CXX) $(libmesh_CXXFLAGS) -c $< -o $(FISPACT_STUB_DIR)/fispact_stub.o
	ar rcs $@ $(FISPACT_STUB_DIR)/fispact_stub.o

ADDITIONAL_LIBS += $(FISPACT_STUB_LIB)

all: $(FISPACT_STUB_LIB)

include            $(FRAMEWORK_DIR)/app.mk

APPLICATION_DIR    := $(CURRENT_DIR)
APPLICATION_NAME   := Fizzy-unit
BUILD_EXEC         := yes

DEP_APPS    ?= $(shell $(FRAMEWORK_DIR)/scripts/find_dep_apps.py $(APPLICATION_NAME))
include $(FRAMEWORK_DIR)/app.mk

# Find all the Fizzy unit test source files and include their dependencies.
fizzy_unit_srcfiles := $(shell find $(CURRENT_DIR)/src -name "*.C")
fizzy_unit_deps := $(patsubst %.C, %.$(obj-suffix).d, $(fizzy_unit_srcfiles))
-include $(fizzy_unit_deps)

###############################################################################
# Additional special case targets should be added here
