# PyGAVel
Module for estimating 1D Velocity Model using Genetic Algorithm.

Summary:

Using this module you can calculate 1D velocity model. The method is based on Genetic Algorithm (GA) for finding the best model from the solution space. There are two objective functions (Hypoellipse, Hypo71) that you can choose according to your dataset and so on. 

Requirements:

- Pyevolve (https://github.com/perone/Pyevolve).
- Hypoellipse (https://pubs.usgs.gov/of/1999/ofr-99-0023/).
- Hypo71 (http://jclahr.com/science/software/hypo71/). This is optional. in a case you want to use Hypo71 as an objective function, you need it.
- Numpy.
- Scipy.
- Matplotlib.

Usage [objective function = Hypoellipse]:

- Copy the module "PyGAVel.py" in a working directory.
- Copy "hypoell-loc.sh", "default.cfg", "hypoel.pha", "hypoel.prm", "hypoel.sta" and "fwd_problem.sh" files into working directory.
- The "hymain" is a executable of Hypoeelipse, has been compiled for Linux64x. You can use it in case you have problem with compiling the source code.
- Copy "par.dat" into working directory and make changes in parameters if needed.
- Run the module.

Hints:

- I suggest you to create a directory and put the "hypoell-loc.sh", "hymain" and "fwd_problem.sh" into it. Add this directory's path to the PATH variable.
- Set "GENSIZE" and "POPSIZE" in "par.dat" not less than 150 to get more reliable results.
- take a look at "Pyevolve" documentation using the following address, to get more familiar with GA algorithm and how it works.
