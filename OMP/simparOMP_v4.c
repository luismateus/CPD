#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <time.h> //just to print the time spent in each function

#define RND0_1 ((double) random() / ((long long)1<<31))
#define G 6.67408e-11
#define EPSLON 0.00000025

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


void init_cell(long aux, double cellX[aux], double cellY[aux], double cellM[aux]) {
    int i;
    #pragma omp parallel for
    for (i = 0; i < aux; i++) {
        cellX[i] = 0;
        cellY[i] = 0;
        cellM[i] = 0;  
    }
}

/* determine the center of mass of each cell */
void massCenter_each_cell(long long npar, long ncell, particle_t *par, double cellX[ncell*ncell], double cellY[ncell*ncell], double cellM[ncell*ncell]) {

    //clock_t begin = clock();

    int n;
    long long i;
    double x,y,m,M;
    long aux = pow(ncell,2);
    
    #pragma omp parallel 
    {
        #pragma omp for private (n)
            for (i = 0; i < npar; i++) {
                n=par[i].c;
                if (!cellM[n]){ 
                    
                    x=par[i].x*par[i].m;
                    #pragma omp atomic
                        cellX[n] += x;
                    y=par[i].y*par[i].m;
                    #pragma omp atomic
                        cellY[n]+=y;
                    m=par[i].m;
                    #pragma omp atomic
                        cellM[n] += m;
                } else {
                    x=par[i].x*par[i].m;
                    #pragma omp atomic
                        cellX[n] += x;
                    y=par[i].y*par[i].m;
                    #pragma omp atomic
                        cellY[n] += y;
                    m=par[i].m;
                    #pragma omp atomic
                        cellM[n] += m; 
                }
            }   
    }
    
    

    // clock_t end = clock();
    // double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    //printf("massCenter_each_cell: %f seconds\n",time_spent);
}

/* compute the gravitational force applied to each particle */
void gforce_each_part(long long npar, long ncell, particle_t *par, double cellX[ncell*ncell], double cellY[ncell*ncell], double cellM[ncell*ncell]) {
    double x, y, f, d,d2;
    int nn, c, nx, ny, vx, vy, n;
    long long i;

    //clock_t begin = clock();
    #pragma omp parallel for private (nn,c,nx,ny,vx,vy,n,x,y,f,d,d2)
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
                x = cellX[nn] - par[i].x - 1;
            } else if (vx + nx > ncell) {
                x = cellX[nn] - par[i].x + 1;
            } else {
                x = cellX[nn] - par[i].x;
            }
            
            if (vy + ny < 0) {
                y = cellY[nn] - par[i].y - 1;
            } else if (vy + ny > ncell) {
                y = cellY[nn] - par[i].y + 1;
            } else {
                y = cellY[nn] - par[i].y;
            }

            d = pow(x, 2) + pow(y, 2);
            d2 = sqrt(d);
            if (d < EPSLON) {
                f = 0;
            } else {
                f = G * (par[i].m * cellM[nn]) / d;   

                par[i].gforcex += x / (d2 / f);
                par[i].gforcey += y / (d2 / f);
            }
        }
    }
    // clock_t end = clock();
    // double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    //printf("gforce_each_part: %f seconds\n",time_spent);
}

/* calculate the new velocity and then the new position of each particle */
void newVelPos_each_part(long long npar, long ncell, particle_t *par) {
    double ax, ay;
    long long i;

    //clock_t begin = clock();

    #pragma omp parallel for private(ax,ay)
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
    // clock_t end = clock();
    // double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    //printf("newVelPos_each_part: %f seconds\n",time_spent);
}

void total_center_of_mass(particle_t *par, long long npar) {
    double x = 0, y = 0, m = 0;
    long long i;
    #pragma omp parallel for reduction(+:x,y,m)
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
    /* input */

    int num_threads = 8; //try to change
    omp_set_num_threads(num_threads);
    if (argc == 5) {
        //printf("10\n");
        fflush(stdout);
        long rand_seed = strtol(argv[1], NULL, 10);
        long grid_size = strtol(argv[2], NULL, 10);
        long long n_part = strtol(argv[3], NULL, 10);
        int time_steps = strtol(argv[4], NULL, 10);

        particle_t *par = (particle_t *)malloc(n_part * sizeof(particle_t));
        if (par == NULL) {
            printf("malloc particle_t failed\n");
            exit(EXIT_FAILURE); 
        }

        init_particles(rand_seed, grid_size, n_part,(particle_t *) par);
        long aux = pow(grid_size,2);
        double cellX[aux], cellY[aux], cellM[aux];

        init_cell(aux, cellX, cellY, cellM);
        //printf("20\n");
        fflush(stdout);


        int t;
        for (t = 0; t < time_steps; t++) {
            //printf("%d\n",t);
            fflush(stdout);
            */
            massCenter_each_cell(n_part, grid_size,(particle_t *) par, cellX, cellY, cellM);
            gforce_each_part(n_part, grid_size,(particle_t *) par, cellX, cellY, cellM);
            newVelPos_each_part(n_part, grid_size,(particle_t *) par);
            init_cell(aux,cellX, cellY, cellM);
        }
        printf("%.2f %.2f\n", par[0].x, par[0].y);
        total_center_of_mass((particle_t *) par, n_part);
        free(par);

    } else { 
        printf("Wrong number of arguments!\n");
    }
    return 0;
}
