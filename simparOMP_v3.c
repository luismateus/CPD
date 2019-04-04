#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <omp.h>
#include <string.h>

#define RND0_1 ((double) random() / ((long long)1<<31))
#define G 6.67408e-11
#define EPSLON 0.0005

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
	    //assert(par[i].c >= 0 && par[i].c < ncside*ncside);
    }
}

void init_cell(cell_t *cell, long cell_n, int num_threads, cell_t *matrix) {
    int i,j, matrix_pos;

    #pragma omp parallel for private(j)
        for (i = 0; i < cell_n; i++) {
            cell[i].x = 0;
            cell[i].y = 0;
            cell[i].m = 0;
            //assert(i >= 0 && i < grid_size*grid_size);
            for (j=0; j< num_threads; j++){
                matrix_pos = i + cell_n * j;
                matrix[matrix_pos].x = 0;
                matrix[matrix_pos].y = 0;
                matrix[matrix_pos].m = 0;
                //assert(i + grid_size * grid_size * j < grid_size * grid_size * num_threads);

            }
        }
}

/* determine the center of mass of each cell */
void massCenter_each_cell(int npar, int cell_n, particle_t *par, cell_t *cell, int num_threads,cell_t *matrix) {
    int j, n, i, thread_num, matrix_pos;


    #pragma omp parallel for private(n)
        for (i = 0; i < npar; i++) {
            n = par[i].c;
            thread_num = omp_get_thread_num();
            matrix_pos = n + cell_n * thread_num;
            if(!matrix[matrix_pos].m){
                matrix[matrix_pos].m = par[i].m;
                matrix[matrix_pos].x = par[i].x*par[i].m;
                matrix[matrix_pos].y = par[i].y*par[i].m;
                //assert(n + ncell * ncell * thread_num < ncell * ncell * num_threads);
            }
            else{
                matrix[matrix_pos].m += par[i].m;
                matrix[matrix_pos].x += par[i].x*par[i].m;
                matrix[matrix_pos].y += par[i].y*par[i].m;
                //assert(n + ncell * ncell * thread_num >= 0);

            }
        }
    
    
    for(j=0; j<num_threads;j++){
        #pragma omp parallel for
            for (n = 0; n < cell_n; n++) {
                matrix_pos = n + cell_n * j;
                cell[n].m+=matrix[matrix_pos].m;
                cell[n].x+=matrix[matrix_pos].x;
                cell[n].y+=matrix[matrix_pos].y;
                //assert(n + ncell * ncell * j < ncell * ncell * num_threads);
                }
        }
    #pragma omp parallel for
        for (n = 0; n < cell_n; n++) {
            if(cell[n].m!=0){
                cell[n].x = cell[n].x/cell[n].m;
                cell[n].y = cell[n].y/cell[n].m;
                //assert(n >= 0 && n < ncell * ncell);
                
            }
        }
}

/* compute the gravitational force applied to each particle */
void gforce_each_part(int npar, int grid_size, particle_t *par, cell_t *cell) {
    double x, y, f, d;
    int nn, c, nx, ny, vx, vy, i, n;

    #pragma omp parallel for private(nn, c, nx, ny, vx, vy, n, x, y, f, d)
        for (i = 0; i < npar; i++) {
            par[i].gforcex = 0;
            par[i].gforcey = 0;
            c = par[i].c;
        

            for (n = 0; n < 9; n++) {
                nx = c % grid_size;
                ny = floor(c / grid_size);

                vx = (n % 3 - 1);
                vy = (floor(n / 3) - 1);

                nn = (nx + vx + grid_size) % grid_size + ((ny + vy + grid_size) % grid_size) * grid_size;

                if (vx + nx < 0) {
                    x = cell[nn].x - par[i].x - 1;
                } else if (vx + nx > grid_size) {
                    x = cell[nn].x - par[i].x + 1;
                } else {
                    x = cell[nn].x - par[i].x;
                }
                
                if (vy + ny < 0) {
                    y = cell[nn].y - par[i].y - 1;  
                } else if (vy + ny > grid_size) {
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
void newVelPos_each_part(int npar, int grid_size, particle_t *par) {
    double ax, ay;
    int i;

    #pragma omp parallel for private(ax, ay)
        for(i = 0; i < npar; i++) {
            ax = par[i].gforcex / par[i].m;
            ay = par[i].gforcey / par[i].m;

            par[i].vx += ax;
            par[i].vy += ay;

            par[i].x += par[i].vx + ax / 2;
            par[i].x = fmod(par[i].x + 1 , 1);

            par[i].y += par[i].vy + ay / 2;
            par[i].y = fmod(par[i].y + 1 , 1);

            par[i].c = (int)floor(par[i].x * grid_size) + ((int)floor(par[i].y * grid_size)) * grid_size;
            //assert(par[i].c >= 0);
        }
}

void total_center_of_mass(particle_t *par, long npar) {
    
    double x = 0, y = 0, m = 0;
    int i;

    #pragma parallel for private(x,y,m)
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
    cell_t *cell, *matrix;
    int t;
    int num_threads = 8; //change this
    omp_set_num_threads(num_threads);

    /* input */
    if (argc == 5) {
        long rand_seed = strtol(argv[1], NULL, 10);
        long grid_size = strtol(argv[2], NULL, 10);
        long long n_part = strtol(argv[3], NULL, 10);
        int time_steps = strtol(argv[4], NULL, 10);

        int cell_n = grid_size * grid_size;
        par = (particle_t *)malloc(n_part * sizeof(particle_t));
        init_particles(rand_seed, grid_size, n_part, par);


        matrix = (cell_t *)malloc(cell_n * num_threads * sizeof(cell_t));
        cell = (cell_t *)malloc(cell_n * sizeof(cell_t));
        init_cell(cell, cell_n,num_threads, matrix);
        

        for (t = 0; t < time_steps; t++) {
            massCenter_each_cell(n_part, cell_n, par, cell, num_threads, matrix);
            gforce_each_part(n_part, grid_size, par, cell);
            newVelPos_each_part(n_part, grid_size, par);
            init_cell(cell, cell_n, num_threads, matrix);
        }
        printf("%.2f %.2f\n", par[0].x, par[0].y);
        total_center_of_mass(par, n_part);

        free(matrix);
        free(par);
        free(cell);
    } else { 
        printf("Wrong number of arguments!\n");
    }
    return 0;
}
