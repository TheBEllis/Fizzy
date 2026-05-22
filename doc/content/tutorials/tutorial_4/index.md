# Fizzy tutorial 4: Distributed Sampling with Cardinal 

In this tutorial, we will set up a MOOSE multiapp using Cardinal, to perform a photon transport solve. The source of the transport solve will be defined by the photon emission spectra calculated by FISPACT-II. The source will also be distributed over MPI ranks, so that different MPI ranks are responsible for sampling from subsections of the domain.

## Prerequisites

1. An MPI compatible MOOSE build
2. A valid Cardinal installation w/ OpenMC & DAGMC
3. A build of FizzyCompiledSource 

The next section will explain how to install the latter. For instructions on how to install Cardinal, please visit Cardinal's own [installation instructions](https://cardinal.cels.anl.gov/start.html).

## Background

!media media/FizzyFlow.png
       id=fizzyflow
       style=width:50%;margin-left:auto;margin-right:auto;
       prefix=Figure
       caption=Evolution of Nb93 atom count over time.

Figure 1 shows the flow of data used when performing distributed sampling with Fizzy. To perform distributed sampling, a custom OpenMC source is used, defined in FizzyCompiledSource. FizzyCompiledSource derived from OpenMC's [CompiledSource](https://docs.openmc.org/en/stable/pythonapi/generated/openmc.CompiledSource.html) functionality, which allows a source term to be defined by a compiled library. We utilise the CompiledSource functionality in order to have granular control over the weighting of sampled particles. This is necessary to perform distributed sampling, as the source terms for different MPI ranks may be of different strengths.

Once Fizzy has completed its solve, the calculated photon emission spectra need to be communicated directly to the CompiledSource. There is not currently a way to do this using MOOSE's Transfers system. Therefore, we utilise an interprocess memory segment to store the photon emission spectra calculated by Fizzy, which we can then read into the CompiledSource at run time.

## Installing FizzyCompiledSource

FizzyCompiledSource is included as a submodule of Fizzy. To build it, first initialise the submodules in the repository.

```language=bash
git submodule git submodule update --init
```

Then change directory into FizzyCompiledSource from the main Fizzy directory, and create a build folder.

```language=bash
cd Fizzy
cd FizzyCompiledSource
mkdir build && cd build
```

Run CMake from the build folder, specifying the directory of your Cardinal build.

```language=bash
cmake -DCARDINAL_DIR=/path/to/cardinal ../
make -j
```

The built library can be found in the build folder at libCompiledSource.so.

## OpenMC inputs id=openmc_inputs

Within the OpenMC Python API, users can set parameters for their OpenMC run. This tutorial won't go over the fundamentals of setting up an OpenMC solve, or installing the OpenMC python module, but for more information, the OpenMC team have extensively documented the Python API, as well as all other parts of [OpenMC](https://docs.openmc.org/en/stable/pythonapi/index.html).

To make use of distributed sampling, users must set the library they compiled in the previous step as their desired source in OpenMC. Below is an example OpenMC python file that shows users how to set the source, and that can be used a base for their own simulations.

!listing /tutorials/tutorial_4/openmc_photon.py 

## Fizzy Inputs 

To prepare the Fizzy input for communicating with Cardinal/OpenMC, the com_photon_flux parameter must be set to True. This prompts Fizzy to put the necessary photon emission rate data into an interprocess memory segment.

!listing /tutorials/tutorial_4/input.i
         block=Problem
         id=newblocks

The other requisite addition to the Fizzy input file in the MultiApps block, within which we define a Cardinal multi-app that runs after Fizzy's solve is completed.

!listing /tutorials/tutorial_4/input.i
         block=MultiApps 
         id=multiapps

The next additions are not required for the simulation to run, but might still be useful. Here we define a variable transfer from Cardinal to OpenMC. The tallies photon flux on each element will be transferred to the Auxiliary Variable photon_flux in Fizzy. 

!listing /tutorials/tutorial_4/input.i
         block=AuxVariables Transfers
         id=newblocks

## Cardinal Inputs

Listing 5 shows the example Cardinal input file. This input defines the mesh we wish to tally on, quantities we wish to tally, and the OpenMC xml files we wish to use as input.


!alert warning
The mesh used in Fizzy and Cardinal must be the same for distributed sampling to work.


!listing /tutorials/tutorial_4/photons_input.i
         id=cardinal

!alert warning 
To make use of distributed sampling, OpenMC must be aware of the mesh used in the Fizzy calculation. For OpenMC to be aware of the mesh, at least one mesh tally must exist on the mesh.


