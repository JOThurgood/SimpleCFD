# Riemann Solver Readme

## ers_euler_1d.f90 

[See the wiki for more info](https://github.com/JOThurgood/SimpleCFD/wiki/ers_euler_1d.f90)

This is an *exact* Riemann solver for the 1D Euler equations.

It solves for user-set initial conditions, and at a user-set time  ouputs spatially sampled data as dat files for easy postprocessing. It also automatically produces plots of u(x), p(x), rho(x) and en(x) (internal energy) at the requested sampling by calling a python script.  

User should edit the following subroutines near the top of the script:

* subroutine initial_conditions
  * Give the solver your left and right states
  * i.e. pressure (pl and pr), density (rhol and rhor), and velocity (ul and ur)
  * Also set gamma.
* subroutine control
  * This controls the sampling of your output 
  * set t to the time at which you want the solution
  * set nx - number of sampling points. 1000 doesnt take long at all. 
  * set x0 - this is the location of the interface at t=0 (the "diagphram")
  * set x_min and x_max - guess... 

You can also if you like 

* Comment out calls to overwrite user initial conditions with the setup tests 1 through 5 (after "call control" near botom), to see the output for e.g. Sod's shock tube, 1,2,3 problem (double rarefaction), and Woodward and Colella blast waves.   
* Change the precision from double to single by changing the parameter "num" near the start of the script. 

Uses the algorithm described in [Toro's textbook](https://www.springer.com/gb/book/9783540252023). 

Automatically does quantitative tests against solutions for five test problems and reports if OK. 

Automatically checks user initial conditions for pressure positivity / evolution of vacuum.
