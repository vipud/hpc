# Hybrid MPI + OPENMP EP Benchmark

This code implements the random-number generator described in the
NAS Parallel Benchmark document RNR Technical Report RNR-94-007.
The code is "embarrassingly" parallel in that no communication is
required for the generation of the random numbers itself. 

The code is implemented using a combination of MPI and OPENMP code. 
We use a master thread to make MPI calls while the other threads carry the computations. 


# How to build

To build the benchmark simply type :

```
make CLASS="Class Name"
````

Class can be A, B, C, D, E, F, S or W without the quotes.


# How to run the benchmark

From the bin folder located at /hpc/nas/c/hybrid/mpi-openmp/bin
Run the executable by specifying the number of processes. For example:

```
mpiexec -n 12 ./ep.A.x
```

# Link to report in Google Document 
http://bit.ly/2qWIpxK

