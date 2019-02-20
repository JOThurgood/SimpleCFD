# CCAPS - Cell centred approximate projection solver

Solves 2D Navier-Stokes equations in one of the following regimes:
* incompressible, uniform density inviscid flow
* incompressible, uniform density flow with viscousity 
* incompressible, variable density inviscid flow with gravity

using a cell-centred approximate projection approach, with optional gradient limiting etc.

See the wiki for full details of usage and test data: https://github.com/JOThurgood/SimpleCFD/wiki/CCAPS-:-Cell-Centred-Approximate-Projection-Solver 

## Compiling

Currently uses cmake to produce a makefile with the proper linking etc.

* Make a directory to keep cmakes junk in
    * `mkdir build` 
    * `cd build`
* Use cmake to produce a makefile (dont have to repeat this often, just if change compiler flags or add files to source), then run the makefile (do this every time you edit the source).
    * `cmake ..`
    * `make`
* Move back up to the main directory and execute CCAPS
    * `cd ..`
    * `./bin/CCAPS`
* There is a script to cleanup the output and the compiled binary
    * `./clean.sh`
