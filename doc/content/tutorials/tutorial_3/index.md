# Fizzy tutorial 3: Transient Simulations

In this tutorial, the input file from [Tutorial 2](tutorials/tutorial_2/index.md) will be developed into a transient simulation, allowing users to visualise how tracked inventory metrics develop over time.

1. [#times]
2. [#exec]

## Getting Times from the FispactSchedule

!listing /tutorials/tutorial_3/input.i
         block=Times
         id=times-block

To run a transient simulation, we need to get our timesteps. FISPACT is different from standard MOOSE, in that the timestepping system is not dynamic. The user sets the times they would like to be solved for in the irradiation schedule, and then FISPACT solves for those times. We need a way of turning the times defined in the schedule into something that can be passed into MOOSE's Transient executioner. For this we use the [FispactScheduleTimes](FispactScheduleTimes.md) object. The only parameter we pass to this is the name of our previously defined FispactSchedule UserObject.

## Adding PostProcessors id=postprocessors

!listing /tutorials/tutorial_3/input.i
         block=Executioner
         id=exec-block

To perform a transient solve firstly we must set our executioner type to Transient. Then we must use use a [TimeSequenceFromTimes](TimeSequenceFromTimes.md) timestepper, and pass in our previously defined FizzyTimes.


## Output id=output

Now that the simulation is transient, we can visualise the number of Nb93 atoms over time.

!media media/Nb93.mp4
       id=Nb93
       style=width:50%;margin-left:auto;margin-right:auto;
       prefix=Figure
       caption=Evolution of Nb93 atom count over time.
       loop=true

!content pagination use_title=True
                    previous=tutorials/tutorial_2/index.md
                    next=tutorials/tutorial_4/index.md
