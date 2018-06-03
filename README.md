# nbody_serial
A serial solver for the N-Body gravitational force problem in 3D.

This code was written as a starting point for Problem Set 6 in
MPCS 51087 at the University of Chicago.

## Compilation
Compilation settings can be configured at the top of the included makefile.
Settings include turning on/off MPI and OpenMP.

$> make

## Runtime Settings

- The total number of bodies to simulate can be set at runtime using the first
command line argument of the program.

- The total number of timesteps to run can be set at runtime using the second
command line argument.

- The number of threads to run per MPI rank can be set using the third
command line argument

## Running
example:

$> ./nbody 1000 10 4

will simulate 1000 bodies for 10 timesteps, using 4 OpenMP threads per rank.
