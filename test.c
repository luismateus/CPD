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
   double gforcex;  
   double gforcey;
   int c; /* cell */
} particle_t;

typedef struct Cell_t {
   double x; /* x position */
   double y; /* y position */
   double m; /* mass */
} cell_t;

void mySum( void *, void *, int *, MPI_Datatype * );

MPI_Datatype MPI_particle_t;
MPI_Datatype MPI_cell_t;
//MPI_Datatype MPI_cell_list;
//MPI_Datatype MPI_particle_array;


void create_mpi_particle(){

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

    //printf("Created MPI_particle_t ...\n");
    // MPI_Type_contiguous(n_part, MPI_particle_t, &MPI_particle_array);
	// MPI_Type_commit(&MPI_particle_array);

}

void create_mpi_cell(double cell_n){

    int items = 3;
    int blocklen[3]= {1,1,1};

    MPI_Datatype types[3] = {MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE};
    
    MPI_Datatype MPI_cell_t;

    MPI_Aint offsets[3] = {offsetof(cell_t,x), offsetof(cell_t,y), offsetof(cell_t,m)};

    MPI_Type_create_struct(items, blocklen, offsets, types, &MPI_cell_t);
	MPI_Type_commit(&MPI_cell_t);

    // MPI_Type_contiguous(cell_n, MPI_cell_t, &MPI_cell_list);
	// MPI_Type_commit(&MPI_cell_list);


    //printf("Created MPI_cell_t ...\n");

}

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
        if(par[i].c < 0 || par[i].c >= 400){
            printf("init_particle n: %d\n", par[i].c);
        }
    }
}

void init_cell(double *cellM, double *cellY, double *cellX, long cell_n, int num_threads) { //, cell_t *matrix
    int i, j;//, matrix_pos;

    for (i = 0; i < cell_n; i++) {
        
        cellM[i] = 0;
        cellX[i] = 0;
        cellY[i] = 0;
        // for (j=0; j< num_threads; j++){
        //     matrix_pos = i + cell_n * j;
        //     matrix[matrix_pos].x = 0;
        //     matrix[matrix_pos].y = 0;
        //     matrix[matrix_pos].m = 0;

        // }
    }
}

/* determine the center of mass of each cell */
void massCenter_each_cell(int npar, int cell_n, particle_t *par, double *cellM, double *cellY, double *cellX, particle_t *par_aux, MPI_Datatype datatype,MPI_Op op) {
    int n, i, matrix_pos, nprocs, rank;
     

    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    
    //MPI_Barrier(MPI_COMM_WORLD);

    //MPI_Scatter(&matrix, cell_n, MPI_cell_t, &aux ,cell_n, MPI_BYTE, 0, MPI_COMM_WORLD); //n_part%4!=0?
    //MPI_cell_list


   // MPI_Bcast(&cell, cell_n, MPI_cell_list,0, MPI_COMM_WORLD); //without matrix?


    for (i = 0; i < npar/4; i++) { // n_part%4!=0
        //printf("i: %d, max_proc: %d\n",i, npar/nprocs * (rank + 1));
        n = par_aux[i].c;

        //printf("nprocs: %d\n",nprocs);
        //printf("rank: %d\n",rank);
        matrix_pos = n;// + cell_n * rank;
        
        if(!cellM[matrix_pos]){
            cellM[matrix_pos] = par_aux[i].m;
            cellX[matrix_pos] = par_aux[i].x*par_aux[i].m;
            cellY[matrix_pos] = par_aux[i].y*par_aux[i].m;
            //assert(par[i].m && par[i].x && par[i].y != 0);

        }else{
            //printf("matrix.m: %f\n",matrix[matrix_pos].m);
            cellM[matrix_pos] += par_aux[i].m;
            cellX[matrix_pos] += par_aux[i].x*par_aux[i].m;
            cellY[matrix_pos] += par_aux[i].y*par_aux[i].m;
            //assert(par[i].m && par[i].x && par[i].y != 0);
        }
    
    }
    printf("HERE!\n");

    /*
    if(rank==0){
        printf("123\n");
        MPI_Reduce(cell,cell,cell_n,datatype,op,0,MPI_COMM_WORLD); //not working
        printf("321\n");
    }
    */
    printf("lala - %f,%d\n",cellM[300],rank);
    printf("size %f, %d\n",(double) (sizeof(cellM) / sizeof(double)),cell_n);
    MPI_Reduce(&cellM,&cellM,cell_n,MPI_DOUBLE, MPI_SUM,0,MPI_COMM_WORLD); //not working
    //MPI_Reduce(cellY,cellY,cell_n,MPI_DOUBLE, MPI_SUM,0,MPI_COMM_WORLD); //not working
    //MPI_Reduce(cellX,cellX,cell_n,MPI_DOUBLE, MPI_SUM,0,MPI_COMM_WORLD); //not working
    printf("4321\n");
    MPI_Barrier(MPI_COMM_WORLD);
    //alterar
    if (rank==0){
        for (n = 0; n < cell_n; n++) {
            if(cellM[n]!=0){
                cellX[n] = cellX[n]/cellM[n];
                cellY[n] = cellY[n]/cellM[n];
            }
        }  
    }  
}

/* compute the gravitational force applied to each particle */
void gforce_each_part(int npar, int grid_size, particle_t *par, double *cellM, double *cellY, double *cellX) {
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
                x = cellM[nn] - par[i].x - 1;
            } else if (vx + nx > grid_size) {
                x = cellM[nn] - par[i].x + 1;
            } else {
                x = cellM[nn] - par[i].x;
            }
            
            if (vy + ny < 0) {
                y = cellY[nn] - par[i].y - 1;  
            } else if (vy + ny > grid_size) {
                y = cellY[nn] - par[i].y + 1;
            } else {
                y = cellY[nn] - par[i].y;
            }

            d = pow(x, 2) + pow(y, 2);
            d2 = sqrt(d);
            if (d < EPSLON) {
                f = 0;
            } else {
                f = (G * cellM[nn]) / d;   

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
        // if(par[i].c <= 0 || par[i].c >= 400){
        //     printf("par[i].c: %d\n",par[i].c);
        // }
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

void mySum(void *invec, void *inoutvec, int *len, MPI_Datatype *dtype)
{
    printf("123");
    cell_t *value = (cell_t*) invec;
    cell_t *res   = (cell_t*) inoutvec;
    int i;
    for ( i=0; i<*len; i++ ) {
        res[i].m += value[i].m;
        res[i].x += value[i].x;
        res[i].y += value[i].y;
    }
}

int main(int argc, char *argv[]) {
    particle_t *par, *par_aux;
    double cellM[400],cellY[400],cellX[400]; //*matrix, *aux;
    int t, nprocs, rank;


    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    //printf("nprocs: %d\n",nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Op CellSum;
    MPI_Op_create( (MPI_User_function *)mySum, 1, &CellSum );
    // MPI_Datatype MPI_particle_t;
    // MPI_Datatype MPI_cell_t;

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

    /* input */
    if (argc == 5) {
        long rand_seed = strtol(argv[1], NULL, 10);
        long grid_size = strtol(argv[2], NULL, 10);
        long long n_part = strtol(argv[3], NULL, 10);
        int time_steps = strtol(argv[4], NULL, 10);

        int cell_n = grid_size*grid_size;
        
        MPI_Barrier(MPI_COMM_WORLD);

        if(rank==0){
            par = (particle_t *)malloc(n_part * sizeof(particle_t));

            par_aux = (particle_t *)malloc(n_part/4 * sizeof(particle_t));

            init_particles(rand_seed, grid_size, n_part, par);
            //matrix = (cell_t *)malloc(cell_n * nprocs * sizeof(cell_t));
    
            //cell = (cell_t *)malloc(cell_n * sizeof(cell_t));
            init_cell(cellM, cellY, cellX, grid_size, nprocs);

        }else{
            par_aux = (particle_t *)malloc(n_part/4 * sizeof(particle_t));

            //cell = (cell_t *)malloc(cell_n * sizeof(cell_t));
            init_cell(cellM, cellY, cellX, grid_size, nprocs);


        }


        //create_mpi_particle(n_part);
        //create_mpi_cell(cell_n);

        MPI_Barrier(MPI_COMM_WORLD); 

        MPI_Scatter(par,n_part/4, MPI_particle_t, par_aux , n_part/4, MPI_particle_t, 0, MPI_COMM_WORLD); //n_part%4!=0?
        
        
        //printf("nprocs: %d\n",nprocs);
        for (t = 0; t < time_steps; t++) {
            printf("size %f, %i\n",(double) (sizeof(cellM) / sizeof(double)),t);
            //MPI_Reduce(cellM,cellM,cell_n,MPI_DOUBLE, MPI_SUM,0,MPI_COMM_WORLD); //not working
            massCenter_each_cell(n_part, cell_n, par, cellM, cellY, cellX,par_aux,MPI_cell_t,CellSum);
            //MPI_Reduce(cell,cell,cell_n,MPI_cell_t,MPI_SUM,0,MPI_COMM_WORLD);
            gforce_each_part(n_part/4, grid_size, par_aux, cellM, cellY, cellX);
            newVelPos_each_part(n_part/4, grid_size, par_aux);
            init_cell(cellM, cellY, cellX, grid_size, nprocs);
        }

        if(rank == 0){
            printf("%.2f %.2f\n", par[0].x, par[0].y);

            //fazer gather do par aqui!

            total_center_of_mass(par, n_part);
            free(par_aux);
            free(par);
            free(cellM);
            free(cellY);
            free(cellX);

        }
        else
        {
            free(par_aux);
            free(cellM);
            free(cellY);
            free(cellX);
        }
        
       
        MPI_Op_free( &CellSum );
        MPI_Finalize();
       
    } else { 
        printf("Wrong number of arguments!\n");
    }
    return 0;
}