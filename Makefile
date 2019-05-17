simpar: 
	mpicc -g -o simpar simpar-mpi.c -lm

clean:
	rm simpar