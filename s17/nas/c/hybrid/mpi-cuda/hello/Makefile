NVCC = /usr/local/cuda-7.5/bin/nvcc
CC = mpicc
LIBS = -L/usr/local/cuda-7.5/lib64 -lcuda -lcudart

all : hello.o

hello2.exec: hello2.cu
	$(CC) hello2.cu -o hello2.exec
run_hello2: hello2.exec
	mpiexec -np 2 hello2.exec
hello.o: hello.cu hello.c
	$(NVCC) --compiler-bindir mpicc hello.c hello.cu -o hello.exec
hello.exec: hello.o
	mpicc hello.o hello_cu.o $(LIBS) -o hello.exec
clean:
	rm -rf *~ *.o *.exec
