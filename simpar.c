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

typedef struct Cell_t {
   double x; // x position
   double y; // y position
   double m; // mass
} cell_t;

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
void massCenter_each_cell(int npar, int ncell, particle_t *par, cell_t *cell) {

    for(int i = 0; i < npar; i++)
    {
        //int n=1;
        int n= (int)floor(par[i].x * ncell) + ((int)floor(par[i].y * ncell))* ncell;

        cell[n].x += par[i].x*par[i].m;
        cell[n].y += par[i].y*par[i].m;
        cell[n].m += par[i].m;
    }
    for(int i = 0; i < ncell * ncell; i++)
    {
        cell[i].x /= cell[i].m;
        cell[i].y /= cell[i].m;
    }

    /*
    for(int i = 0; i < npar; i++)
    {
        //int n=1;
        int n= (int)floor(par[i].x * ncell) + ((int)floor(par[i].y * ncell))* ncell;
        cell[n].m += par[i].m;
        cell[n].x += (par[i].x - cell[n].x)*(par[i].m/ cell[n].m);
        cell[n].y+= (par[i].y - cell[n].y)*(par[i].m/ cell[n].m);
        
    }
    */


}

// compute the gravitational force applied to each particle
void gforce_each_part(particle_t *par) {

}

// calculate the new velocity and then the new position of each particle
void newVelPos_each_part(particle_t *par) {

}

int main(int argc, char *argv[]) {
    particle_t *par;
    int t;

    //input
    if (argc == 5) {
        long rand_seed = strtol(argv[1], NULL, 10);
        long grid_size = strtol(argv[2], NULL, 10);
        long long n_part = strtol(argv[3], NULL, 10);
        int time_steps = strtol(argv[4], NULL, 10);

        par = (particle_t *)malloc(n_part * sizeof(particle_t));

        init_particles(rand_seed, grid_size, n_part, par);

        // for each time-step
        for (t = 0; t < time_steps; t++) {
            massCenter_each_cell(par);
            gforce_each_part(par);
            newVelPos_each_part(par);
        }


    } else { 
        printf("Wrong number of arguments!\n");
    }

    return 0;
}
