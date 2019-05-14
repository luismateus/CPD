#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <mpi.h>
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
    }
}

void init_cell(cell_t *cell, long cell_n, int num_threads, cell_t *matrix) {
    int i, j, matrix_pos;

    for (i = 0; i < cell_n; i++) {
        
        cell[i].x = 0;
        cell[i].y = 0;
        cell[i].m = 0;
        for (j=0; j< num_threads; j++){
            matrix_pos = i + cell_n * j;
            matrix[matrix_pos].x = 0;
            matrix[matrix_pos].y = 0;
            matrix[matrix_pos].m = 0;

        }
    }
}

/* determine the center of mass of each cell */
void massCenter_each_cell(int npar, int cell_n, particle_t *par, cell_t *cell, cell_t *matrix) {
    int n, i, rank, nprocs, matrix_pos;
     
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    //printf("nprocs: %d\n",nprocs);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    //dont think i need/makes cense to scatter
    //MPI_Scatter(matrix, sizeof(matrix)*cell_n*nprocs, MPI_BYTE, &aux ,sizeof(aux)*cell_n, MPI_BYTE, 0, MPI_COMM_WORLD);
    
    for (i = 0 + nprocs * rank; i < npar/nprocs * (rank + 1); i++) {
        printf("i: %d, max_proc: %d\n",i, npar/nprocs * (rank + 1));
        n = par[i].c;
        printf("n: %d\n",n); //problem with n, maybe because of initialization in different processors? idk
        matrix_pos = n + cell_n * rank;
        printf("rank: %d\n",rank);
        if(!matrix[matrix_pos].m){
            matrix[matrix_pos].m = par[i].m;
            matrix[matrix_pos].x = par[i].x*par[i].m;
            matrix[matrix_pos].y = par[i].y*par[i].m;
        }else{
            printf("matrix.m: %f\n",matrix[matrix_pos].m);
            matrix[matrix_pos].m += par[i].m;
            matrix[matrix_pos].x += par[i].x*par[i].m;
            matrix[matrix_pos].y += par[i].y*par[i].m;
        }
    
    }
    printf("HERE!\n");
    //MPI_Reduce(&matrix, &cell,sizeof(cell)*cell_n,MPI_BYTE,MPI_SUM,0,MPI_COMM_WORLD); //not working


    for (n = 0; n < cell_n; n++) {
        if(cell[n].m!=0){
            cell[n].x = cell[n].x/cell[n].m;
            cell[n].y = cell[n].y/cell[n].m;
        }
    }    
}

/* compute the gravitational force applied to each particle */
void gforce_each_part(int npar, int grid_size, particle_t *par, cell_t *cell) {
    double x, y, f, d, d2;
    int nn, c, nx, ny, vx, vy, i, n;

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

            //printf("cell[%i].y: %f\n",nn,cell[nn].y);

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

            d = pow(x, 2) + pow(y, 2);
            d2 = sqrt(d);
            if (d < EPSLON) {
                f = 0;
            } else {
                f = (G * cell[nn].m) / d;   

                par[i].gforcex += x / (d2 / f);
                par[i].gforcey += y / (d2 / f);
            }
        }
    }
}

/* calculate the new velocity and then the new position of each particle */
void newVelPos_each_part(int npar, int grid_size, particle_t *par) {
    int i;


    for(i = 0; i < npar; i++) {

        par[i].vx += par[i].gforcex;
        par[i].vy += par[i].gforcey;

        par[i].x += par[i].vx + par[i].gforcex / 2;
        par[i].x = fmod(par[i].x + 1 , 1);

        par[i].y += par[i].vy + par[i].gforcey / 2;
        par[i].y = fmod(par[i].y + 1 , 1);

        par[i].c = (int)floor(par[i].x * grid_size) + ((int)floor(par[i].y * grid_size)) * grid_size;
        //assert(par[i].c >= 0);
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
    cell_t *cell, *matrix, *aux;
    int t, nprocs, rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    //printf("nprocs: %d\n",nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* input */
    if (argc == 5) {
        long rand_seed = strtol(argv[1], NULL, 10);
        long grid_size = strtol(argv[2], NULL, 10);
        long long n_part = strtol(argv[3], NULL, 10);
        int time_steps = strtol(argv[4], NULL, 10);

        int cell_n = grid_size*grid_size;

        if(rank==0){
            par = (particle_t *)malloc(n_part * sizeof(particle_t));

            init_particles(rand_seed, grid_size, n_part, par);
            matrix = (cell_t *)malloc(cell_n * nprocs * sizeof(cell_t));
    
            cell = (cell_t *)malloc(cell_n * sizeof(cell_t));
            init_cell(cell, grid_size, nprocs, matrix);
        }
        
        MPI_Barrier(MPI_COMM_WORLD);

        for (t = 0; t < time_steps; t++) {
            massCenter_each_cell(n_part, cell_n, par, cell, matrix);
            

            gforce_each_part(n_part, grid_size, par, cell);
            newVelPos_each_part(n_part, grid_size, par);
            init_cell(cell, grid_size, nprocs, matrix);
        }
        printf("%.2f %.2f\n", par[0].x, par[0].y);

        total_center_of_mass(par, n_part);

        MPI_Finalize();
        free(matrix);
        free(par);
        free(cell);
       
    } else { 
        printf("Wrong number of arguments!\n");
    }
    return 0;
}