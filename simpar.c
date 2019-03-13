#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define RND0_1 ((double) random() / ((long long)1<<31))
#define G 6.67408e-11
#define EPSLON 0.01

typedef struct Particle_t {
   double x; /* x position */
   double y; /* y position */
   double vx; /* velocity x */
   double vy; /* velocity y */
   double m; /* mass */
   int c; /* cell */
   double gforcex;
   double gforcey;
} particle_t;

typedef struct Cell_t {
   double x; /* x position */
   double y; /* y position */
   double m; /* mass */
} cell_t;

void init_particles(long seed, long ncside, long long n_part, particle_t *par) {
    long long i;

    srandom(seed);

    for (i = 0; i < n_part; i++) {
        par[i].x = RND0_1;
        par[i].y = RND0_1;
        par[i].vx = RND0_1 / ncside / 10.0;
        par[i].vy = RND0_1 / ncside / 10.0;
        par[i].m = RND0_1 * ncside / (G * 1e6 * n_part);

        par[i].c = (int)floor(par[i].x * ncside) + ((int)floor(par[i].y * ncside)) * ncside;
    }
}

/* determine the center of mass of each cell */
void massCenter_each_cell(int npar, int ncell, particle_t *par, cell_t *cell) {
    int n, i;
    
    for (i = 0; i < npar; i++) {
        n = (int)floor(par[i].x * ncell) + ((int)floor(par[i].y * ncell)) * ncell;

        cell[n].m += par[i].m;
        cell[n].x += (par[i].x - cell[n].x) * (par[i].m / cell[n].m);
        cell[n].y += (par[i].y - cell[n].y) * (par[i].m / cell[n].m);       
    }
}

void init_cell(cell_t *cell, long grid_size) {
    int i;

    for (i = 0; i < grid_size * grid_size; i++) {
        cell[i].x = 0;
        cell[i].y = 0;
        cell[i].m = 0;
    }
}

/* compute the gravitational force applied to each particle */
void gforce_each_part(int npar, int ncell, particle_t *par, cell_t *cell) {
    double x, y, f, d;
    int nn, c, nx, ny, vx, vy, i, n;

    for (i = 0; i < npar; i++) {
        par[i].gforcex = 0;
        par[i].gforcey = 0;
        c = par[i].c;

        for (n = 0; n < 9; n++) {
            nx = c % ncell;
            ny = floor(c / ncell);

            vx = (n % 3 - 1);
            vy = (floor(n / 3) - 1);

            nn = (nx + vx + ncell) % ncell + ((ny + vy + ncell) % ncell) * ncell;

            if (vx + nx < 0) {
                x = cell[nn].x - par[i].x - 1;
            } else if (vx + nx > ncell) {
                x = cell[nn].x - par[i].x + 1;
            } else {
                x = cell[nn].x - par[i].x;
            }
            
            if (vy + ny < 0) {
                y = cell[nn].y - par[i].y - 1;
            } else if (vy + ny > ncell) {
                y = cell[nn].y - par[i].y + 1;
            } else {
                y = cell[nn].y - par[i].y;
            }

            d = sqrt(pow(x, 2) + pow(y, 2));

            if (d < EPSLON) {
                f = 0;
            } else {
                f = G * (par[i].m * cell[nn].m) / pow(d, 2);   

                par[i].gforcex += x / (d / f);
                par[i].gforcey += y / (d / f);
            }
        }
    }
}

/* calculate the new velocity and then the new position of each particle */
void newVelPos_each_part(int npar, int ncell, particle_t *par) {
    double ax, ay;
    int i;

    for(i = 0; i < npar; i++) {
        ax = par[i].gforcex / par[i].m;
        ay = par[i].gforcey / par[i].m;

        par[i].vx += ax;
        par[i].vy += ay;

        par[i].x += par[i].vx + ax / 2;
        par[i].x = fmod(par[i].x, 1);

        par[i].y += par[i].vy + ay / 2;
        par[i].y = fmod(par[i].y, 1);

        par[i].c = (int)floor(par[i].x * ncell) + ((int)floor(par[i].y * ncell)) * ncell;
    }
}

void total_center_of_mass(particle_t *par, long npar) {
    double x = 0, y = 0, m = 0;
    int i;

    for (i = 0; i < npar; i++) {
        m += par[i].m;
        x += par[i].x * par[i].m;
        y += par[i].y * par[i].m;
    }
    x /= m;
    y /= m;
    printf("%.2f %.2f\n", x, y);
}


int main(int argc, char *argv[]) {
    particle_t *par;
    cell_t *cell;
    int t;

    /* input */
    if (argc == 5) {
        long rand_seed = strtol(argv[1], NULL, 10);
        long grid_size = strtol(argv[2], NULL, 10);
        long long n_part = strtol(argv[3], NULL, 10);
        int time_steps = strtol(argv[4], NULL, 10);

        par = (particle_t *)malloc(n_part * sizeof(particle_t));

        init_particles(rand_seed, grid_size, n_part, par);

        cell = (cell_t *)malloc(grid_size * grid_size * sizeof(cell_t));
        init_cell(cell, grid_size);

        for (t = 0; t < time_steps; t++) {
            massCenter_each_cell(n_part, grid_size, par, cell);
            gforce_each_part(n_part, grid_size, par, cell);
            newVelPos_each_part(n_part, grid_size, par);
            init_cell(cell, grid_size);
        }
        printf("%.2f %.2f\n", par[0].x, par[0].y);

        total_center_of_mass(par, n_part);
    } else { 
        printf("Wrong number of arguments!\n");
    }
    return 0;
}
