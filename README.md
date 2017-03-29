Mkay, you should be able run this by just typing make. That compiles sho.f90 and main.f90 for you into an executable called 'main' which it then runs to generate all of that good ol' data. The data tables are saved in the data folder which is then accesed by the 3-N.sh script which is also automatically run after typing 'make'. If you want to plot it yourself though, you can just hop out of the script and go into the 'data folder'.

I've got a bunch of data files that in retrospect I could have consolidated better, but oh well. You can pretty easily figure out what each data file contains by looking at the plotting scripts, which are written in bash/gnuplot. iPython notebook included as well. 

TODO:
write a less terrible README
clean up directories

