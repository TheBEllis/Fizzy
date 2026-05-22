# Fizzy tutorial 1: Just FISPACT 

This tutorial will focus on getting a base level input file prepared for a FISPACT solve using Fizzy.

There are +6 necessary components+ to a Fizzy input file,

1. [#mesh-section]
2. [#problem-section]
3. [#schedule]
4. [#input-flux]
5. [#materials]
6. [#nd]
7. [#exec]

## The Mesh id=mesh-section

!listing /tutorials/tutorial_1/input.i
         block=Mesh
         id=mesh-block
         caption=The +Mesh+ block

A Fizzy input file defines the mesh just like a standard [MOOSE](https://mooseframework.inl.gov/syntax/Mesh/) input file. There are no special requirements for a mesh using exclusively Fizzy, however restrictions do exist when coupling with [Cardinal], these are discussed in latter [tutorials](tutorials/tutorial_2/index.md).

## The Problem id=problem-section

!listing /tutorials/tutorial_1/input.i
         block=Problem
         id=problem-block
         caption=User defined +irradiation schedule+"


The problem block represents the first significant departure from a standard MOOSE input. Here, a problem of type [FispactProblem](source/problems/FispactProblem.md) must be used. The required parameters are listed in [Table 1](#fp-params).

!table id=fp-params caption=Necessary input parameters for a FispactProblem
| Parameter | Description |
| - | - |
| +molar_mass_data+ | Path to an HDF5 file containing molar mass data for relevant nuclides. |
| +fispact_schedule_uo+ | The name of the fispact schedule UserObject the user whishes to use in this simulation. |
| +fispact_nuclear_data_uo+ | The name of the fispact nuclear data UserObject the user whishes to use in this simulation. |
| +fispact_input_flux_uo+ | The name of the UserObject defining the input flux. |
| +output_inventory_time+ | If this is a Steady simulation, this parameter determines which FISPACT-II inventory step is used for output. |


## Irradiation Schedule id=schedule

!listing /tutorials/tutorial_1/input.i
         block=UserObjects/Schedule
         id=schedule-block
         caption=The +irradiation schedule+ block

This block defines the irradiation schedule of your simulation. 

| Parameter | Description |
| - | - |
| +times+ |  The times in seconds used for the irradiation schedule. |
| +flux_amplitude+ | Energy integrated projectile flux values in particles per second per cm^2 |

## Input Flux id=input-flux

!listing /tutorials/tutorial_1/input.i
         block=UserObjects/InputFlux
         id=input-flux-block
         caption=The +input-flux+ block

!alert warning 
Currently only OpenMC statepoint files are supported as inputs.

This block defines the input flux spectra for your FISPACT-II simulation.

## FISPACT material definitions id=materials

!listing /tutorials/tutorial_1/input.i
         block=UserObjects/steel
         id=materials-block
         caption=A FISPACT material definition

This block defines the materials in the simulation. Each block within the mesh must have +one and only one+ FispactMaterial assigned to it. Users can define materials by explicitly stating the concentrations of each isotope, or by defining elements, and letting FISPACT automatically decide the abundancies of each isotope.  

## Nuclear Data id=nd

!listing /tutorials/tutorial_1/input.i
         block=UserObjects/endf_nuclear_data
         id=nd-block
         caption=The user defined +Nuclear Data+ 

Here the user sets the necessary paths to their chosen nuclear data set. The [FispactNuclearData](source/userobjects/FispactNuclearDataPaths.md) UserObject will check whether the paths handed to it are valid, but currently it cannot check the validity of the nuclear data itself.

## Executioner id=exec

!listing /tutorials/tutorial_1/input.i
         block=Executioner
         id=executioner-block
         caption=The +Executioner+ block

For [Steady](source/executioners/Steady.md) simulations, the executioner block looks just like a standard MOOSE input file.

!content pagination use_title=True
                    next=tutorials/tutorial_2/index.md
