# nbody
A solver for the N-Body gravitational force problem in 3D.

This code was written for Problem Set 6 in
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
Usage:

$> ./nbody-[MODE].exe [args]

	$> -n --particles: number of bodies

	$> -i --iterations: number of timesteps

	$> -t --threads: number of threads

example:

$> ./nbody-serial.exe --particles=1000 --iterations=1000

will simulate 1000 bodies for 1000 timesteps, serially.

$> mpirun -n 4 ./nbody-parallel.exe --threads=16 --particles=1000 --iterations=1000

will simulate 1000 bodies for 1000 timesteps, using 16 OpenMP threads per rank.

## Plot
Plotter receives binary data and transforms it into mp4 format file. Default data type is double. If one has a file in float format or different endian, check get_data function and select corresponding solution.

Usage:

$> python plotter.py -f [Datafile] -s [Output].mp4 -n [Particles] -i [Iterations]

example:

$> python plotter.py -f data-nbody-parallel.dat -s vedio-nbody-parallel.mp4 -n 1000 -i 1000

## Correctness
The program generates binary data so one cannot use diff or cmp command in Unix to check whether two files are identical or not.

Use diff.py to check correctness. If nothing printed out, then two files are identical with regard to double.

$> python diff data1 data2
