#include <stdio.h>
#include <stdlib.h>

#define RND0_1 ((double) random() / ((long long)1<<31))
#define G 6.67408e-11
#define EPSLON 0.01

typedef struct Particle_t {
   double x; // x position
   double y; // y position
   double vx; // velocity x
   double vy; // velocity y
   double m; // mass
} particle_t;

void init_particles(long seed, long ncside, long long n_part, particle_t *par) {
    long long i;

    srandom(seed);

    for(i = 0; i < n_part; i++) {
        par[i].x = RND0_1;
        par[i].y = RND0_1;
        par[i].vx = RND0_1 / ncside / 10.0;
        par[i].vy = RND0_1 / ncside / 10.0;

        par[i].m = RND0_1 * ncside / (G * 1e6 * n_part);
    }
}

// determine the center of mass of each cell
void massCenter_each_cell() {

}

// compute the gravitational force applied to each particle
void gforce_each_part() {

}

// calculate the new velocity and then the new position of each particle
void newVelPos_each_part() {

}

int main(int argc, char *argv[]) {
    particle_t *par;
    int t;

    //input
    if (argc == 5) {
        int rand_seed = strtol(argv[1], NULL, 10);
        int grid_size = strtol(argv[2], NULL, 10);
        int n_part = strtol(argv[3], NULL, 10);
        int time_steps = strtol(argv[4], NULL, 10);

        par = (particle_t *)malloc(n_part * sizeof(particle_t));

        init_particles(rand_seed, grid_size, n_part, par);

        // for each time-step
        for (t = 0; t < time_steps; t++) {
            massCenter_each_cell();
            gforce_each_part();
            newVelPos_each_part();
        }


    } else { 
        printf("Wrong number of arguments!\n");
    }

    return 0;
}
