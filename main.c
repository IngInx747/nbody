#include "nbody_header.h"

long nBodies;
int nthreads;
int nIters;
double dt;

void setArgs(int argc, char** argv) {

	/* getopt_long stores the option index here. */
	int option_index = 0;
	int ch;

	static struct option long_options[] = {
		//{"abc", 0|no_argument|required_argument|optional_argument, flag, 'a'},
		{"particles", required_argument, 0, 'n'},
		{"iterations", optional_argument, 0, 'i'},
		{"threads", optional_argument, 0, 't'},
		{0, 0, 0, 0}
	};

    /* Detect the end of the options. */
	while ( (ch = getopt_long(argc, argv, "n:i:t:", long_options, &option_index)) != -1 ) {
		switch (ch) {
			case 'n':
				nBodies = atoi(optarg);
				break;
			case 'i':
				nIters = atoi(optarg);
				break;
			case 't':
				nthreads = atoi(optarg);
				break;
			case '?':
				printf("Unknown option\n");
				break;
			case 0:
				break;
		}
	}
}

int main(int argc, char* argv[]) {

	// Input Parameters
	nBodies = 1000;
	nthreads = 1;
	nIters = 1000;
	dt = 0.2;

	setArgs(argc, argv);

	// Initialize RNG
	srand(42);

	// Set number of OMP threads if necessary
	#ifdef OPENMP
	omp_set_num_threads(nthreads);
	#endif

	// Initialize MPI
	#ifdef MPI
	MPI_Init(&argc, &argv);
	#endif
	
	// Print Inputs
	print_inputs(nBodies, dt, nIters, nthreads);

	// Run Problem
	#ifdef MPI
	run_parallel_problem(nBodies, dt, nIters, "data-nbody-parallel.dat");
	#else
	run_serial_problem(nBodies, dt, nIters, "data-nbody-serial.dat", randomizeBodies_spiral);
	#endif

	// Finalize MPI
	#ifdef MPI
	MPI_Finalize();
	#endif

	return 0;
}
