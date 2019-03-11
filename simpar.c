#include <stdio.h>
#include <stdlib.h>

#include "init_particles.c"

typedef struct Particle_t {
   double x; // x position
   double y; // y position
   double vx; // velocity x
   double vy; // velocity y
   double m; //mass
} particle_t;



//determine the center of mass of each cell
void center_of_mass() {

}

//compute the gravitational force applied to each particle
void gravitational_force() {

}

// calculate the new velocity and then the new position of each particle
void update() {

}

int main(int argc, char *argv[]) {

    particle_t *par;

    //input
    if (argc == 5) {
        int rand_seed = strtol(argv[1], NULL, 10);
        int grid_size = strtol(argv[2], NULL, 10);
        int n_part = strtol(argv[3], NULL, 10);
        int time_steps = strtol(argv[4], NULL, 10);

        par = (particle_t *)malloc(n_part);

        init_particles(rand_seed, grid_size, n_part, par);
        //for each time-step
        center_of_mass();
        gravitational_force();
        update();   


    } else { 
        printf("Wrong number of arguments!\n");
    }

    return 0;
}
