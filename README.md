# wirehexmesher
Tool to create wire-wrap meshes for Nek5000

Requires matlab, python, and a fortran compiler.

Input mesh parameters in `matlab/wire_mesher.m`.

Edit `wire2nek/compile_wire2nek` to use `wire2nek_multi.f` for partial and multiple spans.

Originally created by Paul Fischer, Elia Merzari and Landon Brockmeyer
