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

        //remove this from here
        par[i].c = (int)floor(par[i].x * ncside) + ((int)floor(par[i].y * ncside)) * ncside;
    }
}

void init_cell(cell_t *cell, long grid_size, int num_threads, cell_t *matrix) {
    int i,j;

    #pragma omp parallel for private(j)
        for (i = 0; i < grid_size * grid_size; i++) {
            cell[i].x = 0;
            cell[i].y = 0;
            cell[i].m = 0;
            //assert(i >= 0 && i < grid_size*grid_size);
            //printf("num_threads: %d\n",num_threads);
            for (j=0; j< num_threads; j++){
                //printf("%ld < %ld\n",i + grid_size * grid_size * j,grid_size * grid_size * num_threads);
                matrix[i + grid_size * grid_size * j].x = 0;
                matrix[i + grid_size * grid_size * j].y = 0;
                matrix[i + grid_size * grid_size * j].m = 0;
                //assert(i + grid_size * grid_size * j < grid_size * grid_size * num_threads);
                //printf("i: %d, grid: %ld, j: %d\n", i, grid_size * grid_size, j);
                //assert(i + grid_size * grid_size * j < grid_size * grid_size * num_threads);

            }
        }
}

/* determine the center of mass of each cell */
void massCenter_each_cell(int npar, int ncell, particle_t *par, cell_t *cell, int num_threads,cell_t *matrix) {
    int j, n, i, thread_num;


    #pragma omp parallel for private(n)
        for (i = 0; i < npar; i++) {
            n = par[i].c;
            thread_num = omp_get_thread_num();
            //printf("thread_num: %d\n",thread_num);
            //printf("thread_num: %i ; i: %i\n",thread_num,i);
            if(!matrix[n + ncell * ncell * thread_num].m){
                matrix[n + ncell * ncell * thread_num].m = par[i].m;
                matrix[n + ncell * ncell * thread_num].x = par[i].x*par[i].m;
                matrix[n + ncell * ncell * thread_num].y = par[i].y*par[i].m;
                //assert(n + ncell * ncell * thread_num < ncell * ncell * num_threads);
            }
            else{
                matrix[n + ncell * ncell * thread_num].m += par[i].m;
                matrix[n + ncell * ncell * thread_num].x += par[i].x*par[i].m;
                matrix[n + ncell * ncell * thread_num].y += par[i].y*par[i].m;
                //printf("n: %d\n",n);
                //assert(n + ncell * ncell * thread_num >= 0);

            }
        }
    
    
    for(j=0; j<num_threads;j++){
        #pragma omp parallel for
            for (n = 0; n < ncell*ncell; n++) {
                cell[n].m+=matrix[n + ncell * ncell * j].m;
                cell[n].x+=matrix[n + ncell * ncell * j].x;
                cell[n].y+=matrix[n + ncell * ncell * j].y;
                //assert(n + ncell * ncell * j < ncell * ncell * num_threads);
                // matrix[n + ncell * ncell * j].y=0;
                // matrix[n + ncell * ncell * j].m=0;
                // matrix[n + ncell * ncell * j].x=0;
                }
        }
}

/* compute the gravitational force applied to each particle */
void gforce_each_part(int npar, int ncell, particle_t *par, cell_t *cell) {
    double x, y, f, d, d2;
    int nn, c, nx, ny, vx, vy, i, n;

    #pragma omp parallel for private(nn, c, nx, ny, vx, vy, n, x, y, f, d, d2)
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

                //printf("cell[%i].y: %f\n",nn,cell[nn].y);

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

                d = pow(x, 2) + pow(y, 2);
                d2 = sqrt(d);
                if (d < EPSLON) {
                    f = 0;
                } else {
                    f = G * (par[i].m * cell[nn].m) / d;   

                    par[i].gforcex += x / (d2 / f);
                    par[i].gforcey += y / (d2 / f);
                }
            }
        }
}

/* calculate the new velocity and then the new position of each particle */
void newVelPos_each_part(int npar, int ncell, particle_t *par) {
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

            par[i].c = (int)floor(par[i].x * ncell) + ((int)floor(par[i].y * ncell)) * ncell;
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
    int num_threads = 4; //change this
    omp_set_num_threads(num_threads);
    //printf("num_threads: %d\n",omp_get_num_threads());
    
    /* input */
    if (argc == 5) {
        long rand_seed = strtol(argv[1], NULL, 10);
        long grid_size = strtol(argv[2], NULL, 10);
        long long n_part = strtol(argv[3], NULL, 10);
        int time_steps = strtol(argv[4], NULL, 10);

        printf("here\n");
        par = (particle_t *)malloc(n_part * sizeof(particle_t));
        printf("out\n"); 
        printf("00\n");
        init_particles(rand_seed, grid_size, n_part, par);
        printf("01\n"); 

        matrix = (cell_t *)malloc(grid_size * grid_size * num_threads * sizeof(cell_t));
        printf("02\n"); 
        cell = (cell_t *)malloc(grid_size * grid_size * sizeof(cell_t));
        printf("03\n"); 
        init_cell(cell, grid_size,num_threads, matrix);
        printf("04\n"); 
        
        for (t = 0; t < time_steps; t++) {
            printf("%d\n",t);
            fflush(stdout);
            massCenter_each_cell(n_part, grid_size, par, cell, num_threads, matrix);
            gforce_each_part(n_part, grid_size, par, cell);
            newVelPos_each_part(n_part, grid_size, par);
            init_cell(cell, grid_size, num_threads, matrix);
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
