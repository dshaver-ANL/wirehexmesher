# wirehexmesher
Tool to create wire-wrap meshes for Nek5000

Requires matlab, python, and a fortran compiler.
Defaults to gfortran.
For ifort, edit `wire2nek/compile_wire2nek`

Input mesh parameters in `matlab/wire_mesher.m`.

Defaults to a single span with periodic BCs. 
For multiple/partial spans edit `nreps` and `partial` variables in `wire2nek/wire2nek.f`.
For inlet/outlet BCs set `ifperiodic` to false in `wire2nek/wire2nek.f`.

Boundary conditions are assigned as:

- pin/wire = 'W  ','I  ' BCID = 1
- outer hex can = 'W  ','I  ', BCID = 2
- inlet = 'v  ','t  ', BCID = 3 (or 'P  ', 'P  ', BCID = 0 for periodic)
- outlet = 'O  ','I  ', BCID = 4 (or 'P  ', 'P  ', BCID = 0 for periodic)

Originally created by Paul Fischer, Elia Merzari and Landon Brockmeyer
