CC = mpicc -std=c99 
EXEC = mpiexec -n 4

all : hello.exec hello_add.exec

hello_add.exec : hello_add.c
	$(CC) hello_add.c -o hello_add.exec

hello.exec : hello.c
	$(CC) hello.c -o hello.exec

clean:
	rm -rf *~ *.exec
