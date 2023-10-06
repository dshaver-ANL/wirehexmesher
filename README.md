# wirehexmesher
Tool to create wire-wrap meshes for Nek5000

Requires matlab, python, and a fortran compiler.
Defaults to `gfortran`.
For `ifort`, edit `wire2nek/compile_wire2nek`

Input mesh parameters in `matlab/wire_mesher.m`.

Defaults to a single span with periodic BCs. 
For multiple/partial spans edit `nreps` and `partial` variables in `wire2nek/wire2nek.f`.
For inlet/outlet BCs set `ifperiodic` to false in `wire2nek/wire2nek.f`.

Boundary conditions are assigned as:

Boundary | velocity | temperature | Boundary ID
---|---|---|---
pin/wire | 'W  ' | 'I  '| 1
outer hex can | 'W  ' | 'I  '| 2
inlet | 'v  ' ('P  ')| 't  ' ('P  ') | 3 (0)
outlet | 'O  ' ('P  ') | 'I  ' ('P  ') | 4 (0)

- Values in parantheses are use when `ifperiodic` is set to true.

Originally created by Paul Fischer, Elia Merzari and Landon Brockmeyer
