# CCAPS_vardens - Cell-Centred Approximate Projection Solver

Version of CCAPS for incompressible, variable density flow.

It is the same solver as the main CCAPS, except the code strips away all other run modes / solvers / options (e.g. homogenous density simplification mode, homogeneous + viscosity).

This is done as a preliminary step in anticipation of writing a version of CCAPS for background subtracted models, (it was becoming to much of an if-nested nightmare to have multiple models / solvers within the one code, so given the purpose of "SimpleCFD" it seemed better to simply produce different versions of CCAPS for each regime, at least for now)

See the wiki of the main code for full details of usage and test data: https://github.com/JOThurgood/SimpleCFD/wiki/CCAPS-:-Cell-Centred-Approximate-Projection-Solver 

