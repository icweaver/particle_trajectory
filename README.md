Code for calculating particle trajectories in generally binary system. Also includes artificial energy loss implementation to simulate disk formation.  

Should be able to start by running **make**. This compiles sho.f90 and main.f90 into an executable called **main** which then runs to generate the output. The data tables are saved in the data folder which is then accesed by the 3-N.sh script -- which is also automatically run after typing 'make'. 

Plotting can also be done manually in the **data** folder.

#### TODO:
* clean up directories
* generalize binary parameters w/ user input
