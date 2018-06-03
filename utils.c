#include "nbody_header.h"

double get_time(void) {

	#ifdef MPI
	return MPI_Wtime();
	#endif

	#ifdef OPENMP
	return omp_get_wtime();
	#endif

	time_t time;
	time = clock();

	return (double) time / (double) CLOCKS_PER_SEC;
}

void print_inputs(long nBodies, double dt, int nIters, int nthreads ) {
	
	int mype = 0;
	int nprocs = 1;

	#ifdef MPI
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &mype);
	#endif

	if( mype == 0 )
	{
		printf("INPUT PARAMETERS:\n");
		printf("N Bodies =                     %ld\n", nBodies);
		printf("Timestep dt =                  %.3le\n", dt);
		printf("Number of Timesteps =          %d\n", nIters);
		#ifdef OPENMP
		printf("Number of Threads per Rank =   %d\n", nthreads);
		#endif
		#ifdef MPI
		printf("Number of MPI Ranks =          %d\n", nprocs);
		#endif
		printf("BEGINNING N-BODY SIMLUATION\n");
	}
}

double uniform_rand(double mu, double sigma) {

	/**
	* Generate random number normally distributed
	* Don't use this function for debug propose as it calls rand()
	*/

	double U1, U2, W, mult;
	static double X1, X2;
	static int call = 0;
 
	if (call == 1) {
		call = !call;
		return (mu + sigma * (double) X2);
    }
 
	do {
		U1 = -1 + ((double) rand () / RAND_MAX) * 2;
		U2 = -1 + ((double) rand () / RAND_MAX) * 2;
        W = powf (U1, 2) + pow (U2, 2);
	} while (W >= 1 || W == 0);
 
	mult = sqrtf ((-2 * logf (W)) / W);
	X1 = U1 * mult;
	X2 = U2 * mult;
 
	call = !call;
 
	return (mu + sigma * (double) X1);
}
