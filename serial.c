#include "nbody_header.h"

void run_serial_problem(int nBodies, double dt, int nIters, char * fname,
	void (*randomizeBodies)(Body *, int)) {

	// Open File and Write Header Info
	int fd = open(fname, O_RDWR|O_CREAT|O_TRUNC, S_IRUSR|S_IWUSR|S_IRGRP|S_IROTH);

	// Allocate Bodies
	Body * bodies  = (Body *) calloc( nBodies, sizeof(Body) );
	randomizeBodies(bodies, nBodies);

	double start = get_time();

	// Loop over timesteps
	for (int iter = 0; iter < nIters; iter++) {

		if ((iter*100) % nIters == 0) write(STDOUT_FILENO, ".", 2);

		// Output body positions to file
		for (int i = 0; i < nBodies; i++)
			write(fd, bodies + i, 3 * sizeof(double));

		// Compute new forces & velocities for all particles
		compute_forces(bodies, dt, nBodies);

		// Update positions of all particles
		for (int i = 0 ; i < nBodies; i++)
		{
			bodies[i].x += bodies[i].vx*dt;
			bodies[i].y += bodies[i].vy*dt;
			bodies[i].z += bodies[i].vz*dt;
		}

	}

	// Close data file
	close(fd);

	double stop = get_time();

	double runtime = stop-start;
	double time_per_iter = runtime / nIters;
	long interactions = nBodies * nBodies;
	double interactions_per_sec = interactions / time_per_iter;

	printf("SIMULATION COMPLETE\n");
	printf("Runtime [s]:              %.3le\n", runtime);
	printf("Runtime per Timestep [s]: %.3le\n", time_per_iter);
	printf("interactions:             %d\n", nIters);
	printf("Interactions per sec:     %.3le\n", interactions_per_sec);
	fflush(stdout);

	free(bodies);
}

// Computes the forces between all bodies and updates
// their velocities accordingly
void compute_forces(Body * bodies, double dt, int nBodies)
{
	double G = 6.67259e-3;
	double softening = 1.0e-5;

	// For each particle in the set
	for (int i = 0; i < nBodies; i++)
	{ 
		double Fx = 0.0;
		double Fy = 0.0;
		double Fz = 0.0;

		// Compute force from all other particles in the set
		for (int j = 0; j < nBodies; j++)
		{
			// F_ij = G * [ (m_i * m_j) / distance^3 ] * (location_j - location_i) 

			// First, compute the "location_j - location_i" values for each dimension
			double dx = bodies[j].x - bodies[i].x;
			double dy = bodies[j].y - bodies[i].y;
			double dz = bodies[j].z - bodies[i].z;

			// Then, compute the distance^3 value
			// We will also include a "softening" term to prevent near infinite forces
			// for particles that come very close to each other (helps with stability)

			// distance = sqrt( dx^2 + dx^2 + dz^2 )
			double distance = sqrt(dx*dx + dy*dy + dz*dz + softening);
			double distance_cubed = distance * distance * distance;

			// Now compute G * m_2 * 1/distance^3 term, as we will be using this
			// term once for each dimension
			// NOTE: we do not include m_1 here, as when we compute the change in velocity
			// of particle 1 later, we would be dividing this out again, so just leave it out
			double m_j = bodies[j].mass;
			double mGd = G * m_j / distance_cubed;
			Fx += mGd * dx;
			Fy += mGd * dy;
			Fz += mGd * dz;
		}

		// With the total forces on particle "i" known from this batch, we can then update its velocity
		// v = (F * t) / m_i
		// NOTE: as discussed above, we have left out m_1 from previous velocity computation,
		// so this is left out here as well
		bodies[i].vx += dt*Fx;
		bodies[i].vy += dt*Fy;
		bodies[i].vz += dt*Fz;
	}
}

// Randomizes all bodies to the following default criteria
// Locations (uniform random between -1.0 < r < 1.0 )
// Velocities (uniform random between -1.0e-3 < r < 1.0e3 )
// Masses (all equal at 1.0 / nBodies)
// You should make this more exotic
void randomizeBodies_default(Body * bodies, int nBodies) {

	// velocity scaling term
	double vm = 1.0e-3;

	for (int i = 0; i < nBodies; i++) {
		// Initialize position between -1.0 and 1.0
		bodies[i].x = 2.0 * (rand() / (double)RAND_MAX) - 1.0;
		bodies[i].y = 2.0 * (rand() / (double)RAND_MAX) - 1.0;
		bodies[i].z = 2.0 * (rand() / (double)RAND_MAX) - 1.0;

		// Intialize velocities
		bodies[i].vx = 2.0*vm * (rand() / (double)RAND_MAX) - vm;
		bodies[i].vy = 2.0*vm * (rand() / (double)RAND_MAX) - vm;
		bodies[i].vz = 2.0*vm * (rand() / (double)RAND_MAX) - vm;

		// Initialize masses so that total mass of system is constant
		// regardless of how many bodies are simulated
		bodies[i].mass = 1.0 / nBodies;
	}
}

void randomizeBodies_spiral(Body * bodies, int nBodies) {

    double pi_double = 3.14159265359;
    // velocity scaling term
	double vm = 0.5e-1;

    for (int i = 0; i < nBodies; i++) {

    	double radius = 2.0 * uniform_rand(0.0, 1.0);
    	radius = (radius>0)?radius:(-radius);
        double theta = pi_double * (2.0 * (rand() / (double)RAND_MAX) - 1.0);
        double tmpx = radius * sinf(theta);
        double tmpy = radius * cosf(theta);

        // Initialize position between -1.0 and 1.0
		bodies[i].x = tmpx;
		bodies[i].y = tmpy;
		bodies[i].z = 0.01 * (2.0 * (rand() / (double)RAND_MAX) - 1.0);

		// Intialize velocities
		bodies[i].vx = vm * (tmpy/radius  + 0.01 * (2.0 * (rand() / (double)RAND_MAX) - 1.0));
		bodies[i].vy = vm * (-tmpx/radius + 0.01 * (2.0 * (rand() / (double)RAND_MAX) - 1.0));
		bodies[i].vz = vm * (0.01 * (2.0 * (rand() / (double)RAND_MAX) - 1.0f));

		// Initialize masses so that total mass of system is constant
		// regardless of how many bodies are simulated
		bodies[i].mass = 1.0 / nBodies;
    }
}

void randomizeBodies_clusters(Body * bodies, int nBodies) {

	// velocity scaling term
	double vm = 1.0e-3;

	double deviation = 2.0;

	for (int i = 0; i < nBodies; i++) {
		int tmpx = (i%2)?(deviation):(-deviation);
		// Initialize position between -1.0 and 1.0
		bodies[i].x = 2.0 * (rand() / (double)RAND_MAX) - 1.0 + tmpx;
		bodies[i].y = 2.0 * (rand() / (double)RAND_MAX) - 1.0;
		bodies[i].z = 2.0 * (rand() / (double)RAND_MAX) - 1.0;

		// Intialize velocities
		bodies[i].vx = vm * (2.0*(rand() / (double)RAND_MAX) - 1.0);
		bodies[i].vy = vm * (2.0*(rand() / (double)RAND_MAX) - 1.0);
		bodies[i].vz = vm * (2.0*(rand() / (double)RAND_MAX) - 1.0);

		// Initialize masses so that total mass of system is constant
		// regardless of how many bodies are simulated
		bodies[i].mass = 1.0 / nBodies;
	}
}
