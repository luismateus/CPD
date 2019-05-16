#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <mpi.h>
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
   double gforcex;  
   double gforcey;
   int c; /* cell */
} particle_t;

typedef struct Cell_t {
   double x; /* x position */
   double y; /* y position */
   double m; /* mass */
} cell_t;

void mySum( cell_t *, cell_t *, int *, MPI_Datatype * );

MPI_Datatype MPI_particle_t;
MPI_Datatype MPI_cell_t;


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

void init_cell(cell_t *cell, long cell_n) { //, cell_t *matrix
    int i;

    #pragma omp parallel for
        for (i = 0; i < cell_n; i++) {
            cell[i].x = 0;
            cell[i].y = 0;
            cell[i].m = 0;
        }
}

/* determine the center of mass of each cell */
void massCenter_each_cell(int npar, int cell_n, particle_t *par, cell_t *cell, particle_t *par_aux, MPI_Datatype datatype,MPI_Op op,cell_t *cell_sum) {
    int n, i, matrix_pos, nprocs, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    for (i = 0; i < npar; i++) { // n_part%4!=0
        
        //printf("proc: %d, par_aux[i].c: %d\n",rank,par_aux[i].c);
        matrix_pos = par_aux[i].c;// + cell_n * rank;
        if(matrix_pos == -10){
            //printf("proc:%d; i: %d, par_aux[i].c: %d\n",rank,i,par_aux[i].c);
            continue;
        }
       
        if(!cell[matrix_pos].m){
            cell[matrix_pos].m = par_aux[i].m;
            cell[matrix_pos].x = par_aux[i].x*par_aux[i].m;
            cell[matrix_pos].y = par_aux[i].y*par_aux[i].m;
            //assert(par[i].m && par[i].x && par[i].y != 0);

        }else{
            //printf("matrix.m: %f\n",matrix[matrix_pos].m);
            cell[matrix_pos].m += par_aux[i].m;
            cell[matrix_pos].x += par_aux[i].x*par_aux[i].m;
            cell[matrix_pos].y += par_aux[i].y*par_aux[i].m;
            //assert(par[i].m && par[i].x && par[i].y != 0);
        }
    
    }
    MPI_Reduce(cell,cell_sum,cell_n,datatype,op,0,MPI_COMM_WORLD); 

    if (rank==0){
        #pragma omp parallel for
        for (n = 0; n < cell_n; n++) {
            if(cell_sum[n].m!=0){
                cell_sum[n].x = cell_sum[n].x/cell_sum[n].m;
                cell_sum[n].y = cell_sum[n].y/cell_sum[n].m;
            }
        }
        cell=cell_sum;
    }
    MPI_Bcast(cell, cell_n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

}

/* compute the gravitational force applied to each particle */
void gforce_each_part(int npar, int grid_size, particle_t *par, cell_t *cell) {
    double x, y, f, d, d2;
    int nn, c, nx, ny, vx, vy, i, n;
    //printf("2\n");

    #pragma omp parallel for private(x, y, f, d, d2,n, c, nx, ny, vx, vy)
    for (i = 0; i < npar; i++) {
        par[i].gforcex = 0;
        par[i].gforcey = 0;
        c = par[i].c;

        if(c == -10){
            continue;
        }
    

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
    //printf("3\n");
    #pragma omp parallel for
    for(i = 0; i < npar; i++) {

        if(par[i].c == -10){
            continue;
        }

        par[i].vx += par[i].gforcex;
        par[i].vy += par[i].gforcey;

        par[i].x += par[i].vx + par[i].gforcex / 2;
        par[i].x = fmod(par[i].x + 1 , 1);

        par[i].y += par[i].vy + par[i].gforcey / 2;
        par[i].y = fmod(par[i].y + 1 , 1);


        
        par[i].c = (int)floor(par[i].x * grid_size) + ((int)floor(par[i].y * grid_size)) * grid_size;
    }
}

void total_center_of_mass(particle_t *par, long npar,int rank) {
    double x = 0, y = 0, m = 0;
    int i;
    //printf("4\n");

    #pragma omp parallel for reduction(+:x,y,m)
    for (i = 0; i < npar; i++) {

        if(par[i].c == -10){
            continue;
        }

        m += par[i].m;
        x += par[i].x * par[i].m;
        y += par[i].y * par[i].m;
    }
    double x2,y2,m2;
    
    MPI_Reduce(&x,&x2,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&y,&y2,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&m,&m2,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    if (rank==0){
        x2 /= m2;
        y2 /= m2;
        printf("%.2f %.2f\n", par[0].x, par[0].y);
        printf("%.2f %.2f\n", x2, y2);
    }
}

void mySum(cell_t *invec, cell_t *inoutvec, int *len, MPI_Datatype *dtype)
{
    int i;
    #pragma omp parallel for
    for ( i=0; i<*len; i++ ) {
        inoutvec[i].m += invec[i].m;
        inoutvec[i].x += invec[i].x;
        inoutvec[i].y += invec[i].y;
    }
}

int main(int argc, char *argv[]) {
    particle_t *par, *par_aux;
    cell_t *cell;
    cell_t *cell_sum;
    int t, i, nprocs, rank;
    long malloc_count;


    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    //printf("nprocs: %d\n",nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Op CellSum;
    MPI_Op_create( (MPI_User_function *)mySum, 1, &CellSum );

    int items = 8;
    int blocklen[8]= {1,1,1,1,1,1,1,1};

    MPI_Datatype types[8] = {MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,\
    MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_INT};
    MPI_Datatype MPI_particle_t;
    MPI_Aint offsets[8] = {offsetof(particle_t,x), offsetof(particle_t,y), \
    offsetof(particle_t,vx), offsetof(particle_t,vy), offsetof(particle_t,m), \
    offsetof(particle_t,gforcex), offsetof(particle_t,gforcey), offsetof(particle_t,c)};
    MPI_Type_create_struct(items, blocklen, offsets, types, &MPI_particle_t);
	MPI_Type_commit(&MPI_particle_t);

    int items2 = 3;
    int blocklen2[3]= {1,1,1};
    MPI_Datatype types2[3] = {MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE};
    MPI_Datatype MPI_cell_t;
    MPI_Aint offsets2[3] = {offsetof(cell_t,x), offsetof(cell_t,y), offsetof(cell_t,m)};
    MPI_Type_create_struct(items2, blocklen2, offsets2, types2, &MPI_cell_t);
	MPI_Type_commit(&MPI_cell_t);



    /* input */
    if (argc == 5) {
        long rand_seed = strtol(argv[1], NULL, 10);
        long grid_size = strtol(argv[2], NULL, 10);
        long long n_part = strtol(argv[3], NULL, 10);
        int time_steps = strtol(argv[4], NULL, 10);

        int cell_n = grid_size*grid_size;
        
        

        if(n_part%nprocs != 0){
            malloc_count = 1 + ((n_part - 1) / nprocs);
        }else{
            malloc_count = n_part/nprocs;
        }

        MPI_Barrier(MPI_COMM_WORLD);

        if(rank==0){
            par = (particle_t *)malloc(n_part * sizeof(particle_t));

            par_aux = (particle_t *)malloc(malloc_count * sizeof(particle_t));

            init_particles(rand_seed, grid_size, n_part, par);
            //matrix = (cell_t *)malloc(cell_n * nprocs * sizeof(cell_t));
            cell_sum=(cell_t *)malloc(cell_n * sizeof(cell_t));
            cell = (cell_t *)malloc(cell_n * sizeof(cell_t));
            init_cell(cell_sum, cell_n);
            init_cell(cell, cell_n);
            //printf("calc: %lld / %d = %ld ; rest: %lld\n",n_part,nprocs,malloc_count,n_part%nprocs);


        }else{
            par_aux = (particle_t *)malloc(malloc_count * sizeof(particle_t));
            cell = (cell_t *)malloc(cell_n * sizeof(cell_t));
            cell_sum=cell;
            init_cell(cell, cell_n);
        }        
        
        MPI_Scatter(par,malloc_count, MPI_particle_t, par_aux , malloc_count, MPI_particle_t, 0, MPI_COMM_WORLD); //n_part%4!=0?
        
        for(i = 0; i< n_part%nprocs ; i++){
            //printf("proc: %d, par_aux[i].c: %d\n",rank,par_aux[i].c);
            if(rank == nprocs-1){
                par_aux[malloc_count-i-1].c = -10;
                //printf("proc: %d, i: %d,par_aux[i].c: %d\n",rank,i,par_aux[malloc_count-i-1].c);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);


        //printf("nprocs: %d\n",nprocs);
        for (t = 0; t < time_steps; t++) {
            massCenter_each_cell(malloc_count, cell_n, par, cell,par_aux,MPI_cell_t,CellSum,cell_sum);
            //MPI_Reduce(cell,cell,cell_n,MPI_cell_t,MPI_SUM,0,MPI_COMM_WORLD);
            gforce_each_part(malloc_count, grid_size, par_aux, cell);
            newVelPos_each_part(malloc_count, grid_size, par_aux);
            init_cell(cell, cell_n);
        }
        MPI_Barrier(MPI_COMM_WORLD);

        total_center_of_mass(par_aux, malloc_count, rank);
        
        if(rank == 0){
            free(par_aux);
            free(par);
            free(cell);
            free(cell_sum);

        }
        else
        {
            free(par_aux);
            free(cell);
        }
        
       
        MPI_Op_free( &CellSum );
        MPI_Finalize();
       
    } else { 
        printf("Wrong number of arguments!\n");
    }
    return 0;
}