#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <assert.h>
#include <stddef.h>
#include <fcntl.h>
#include <time.h>
#ifdef MPI
#include <mpi.h>
#endif
#ifdef OPENMP
#include <omp.h>
#endif

#ifndef PI
#define PI (3.14159265359);
#endif

typedef struct {
	// Location r_i = (x,y,z)
	double x;
	double y;
	double z;
	// Velocity v_i = (vx, vy, vz)
	double vx;
	double vy;
	double vz;
	// Mass
	double mass;
} Body;

// serial.c
void run_serial_problem(int nBodies, double dt, int nIters, char * fname,
	void (*func)(Body *, int));
void compute_forces(Body * bodies, double dt, int nBodies);
void randomizeBodies_default(Body * bodies, int nBodies);
void randomizeBodies_spiral(Body * bodies, int nBodies);
void randomizeBodies_clusters(Body * bodies, int nBodies);

// utils.c
double get_time(void);
void print_inputs(long nBodies, double dt, int nIters, int nthreads );
double uniform_rand(double mu, double sigma);

// parallel.c
#ifdef MPI
void run_parallel_problem(long nBodies, double dt, long nIters, char * fname);
void compute_forces_multi_set(Body * local, Body * remote, double dt, long nlocal, long nremote);
void parallel_randomizeBodies(Body * bodies, long nBodies, int mype, int nprocs);
void distributed_write_timestep(Body * local_bodies, long nBodies,
	long timestep, long nIters, int nprocs, int mype, MPI_File * fh);
long jobs_manage(long nBodies, int mype, int nprocs, long * iptr);
#endif
