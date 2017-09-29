# Makefile for MPI
ALL : hello

# Generates hello executable
hello: hello.c
	mpicc -o hello hello.c

# Clean exec files
clean:
	rm -rf hello *.o *.exe
