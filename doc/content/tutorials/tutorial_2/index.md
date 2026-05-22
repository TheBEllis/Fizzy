# Fizzy tutorial 2: Adding AuxVariables, AuxKernels, and PostProcessors 

In this tutorial, the input file from [Tutorial 1](tutorials/tutorial_1/index.md) will be developed to add AuxVariables, AuxKernels and PostProcessors to help extract data from the calculated nuclear inventories.

1. [#auxvars]
2. [#postprocessors]
3. [#outputs]

## Adding AuxVars and AuxKernels id=auxvars

!listing /tutorials/tutorial_2/input.i
         block=AuxVariables 
         id=aux-block

Firstly, we add AuxVariables into the input file. Anyone who has used MOOSE before will be familiar with this syntax. We just define a block called `AuxVariables`, which lets MOOSE know we wish to define AuxVariables here. Then we define sub-blocks that represent our actual variables. Variables used in Fizzy AuxKernels must be CONSTANT order MONOMIALS.


!listing /tutorials/tutorial_2/input.i
         block=AuxKernels
         id=auxkern-block

Secondly, we add AuxKernel blocks, to which we assign AuxVariables. Here we add two slightly different AuxKernels, a [FispactNuclideKernel](source/auxkernels/FispactNuclideKernel.md), and a [FispactElementKernel](source/auxkernels/FispactElementKernel.md). 

The FispactNuclideKernel is used for extracting metrics regarding particular nuclides from the nuclear inventory. In this example, we extract the number of Niobium-93 atoms. Because we are using AuxVariables, this value will be evaluated for each element.

The FispactElementKernel is used for getting quantities representitive of the entire nuclear inventory, not just individual nuclides. Here we are modelling the total number of atoms in the entire inventory. 

To get a list of possible values of the `metric`, parameter, see the individual documentation pages for both [FispactNuclideKernel](source/auxkernels/FispactNuclideKernel.md), and [FispactElementKernel](source/auxkernels/FispactElementKernel.md). 


## Adding PostProcessors id=postprocessors

!listing /tutorials/tutorial_2/input.i
         block=Postprocessors
         id=post-block

To calculate scalar quantities over the entire domain, or subsections of the domain, users can use PostProcessors.

This example includes two postprocessors. Firstly, there is an [ElementPhotonEmission](ElementPhotonEmission.md) PostProcessor. This will return the total gamma emission rate summed over the elements defined by the element_ids parameter. If the user wishes to obtain the total gamma emission rate over a subdomain, they can use the [SubdomainPhotonEmission](SubdomainPhotonEmission.md) PostProcessor.

A [NuclideMetric](NuclideMetric.md) is included to track domain wide metrics from particular nuclides. In this example, the nuclide_contribution postprocessor outputs the total slab dose rate contributed by Cr51.

!alert warning
When providing nuclide names to AuxKernels or PostProcessors, make sure they are formatted correctly, using the correct case.

## Output id=output

!listing /tutorials/tutorial_2/input.i
         block=Outputs
         id=outputs-block
         caption=The +Outputas+ block

Now that some variables exist on the mesh, we may want to visualise the results.To do this we will need to add an exodus output block, or in this case, use the syntactical shortcut for one. 

!content pagination use_title=True
                    previous=tutorials/tutorial_1/index.md
                    next=tutorials/tutorial_3/index.md
