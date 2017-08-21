# PyGAVel
Module for estimating 1D Velocity Model using Genetic Algorithm.

Summary:

Using this module you can calculate 1D velocity model. The method is based on Genetic Algorithm (GA) for finding the best model from the solution space. There are two objective function (Hypoellipse, Hypo71) that you can choose according to your dataset and so on. 

Requirements:

- Pyevolve (https://github.com/perone/Pyevolve).
- Hypoellipse (https://pubs.usgs.gov/of/1999/ofr-99-0023/).
- Hypo71 (http://jclahr.com/science/software/hypo71/). This is optional. in a case you want to use Hypo71 as an objective function, you needt it.
- Numpy.
- Scipy.
- Matplotlib.

Usage [objective function = Hypoellipse]:

- Copy the module "PyGAVel.py" in a working directory.
- Copy "hypoell-loc.sh", "default.cfg", "hypoel.pha", "hypoel.prm", "hypoel.sta" and "fwd_problem.sh" files into working directory.
- Copy "par.dat" into working directory and make changes in parameters if needed.
- Run the module.
