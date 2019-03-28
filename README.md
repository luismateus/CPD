# CPD


OpenMP command:
    compile: 
        "gcc -fopenmp simparOMP.C -lm -o simparOMP"
    run: 
        "./simparOMP 1 10 2000000 10"

OmpP command:
    compile:
        "kinst-ompp gcc -fopenmp simparOMP.C -o simparOMP -lm -lstdc++"
    run:
        "./simparOMP 1 10 2000000 10" and read outuput txt file 
    





