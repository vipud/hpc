CC = mpicc 
EXEC = mpiexec

all: run

hello_add.exec : hello_add.c
	$(CC) hello_add.c -o hello_add.exec

hello.exec : hello.c
	$(CC) hello.c -o hello.exec

run: hello_add.exec
	$(EXEC) -n 1  hello_add.exec
clean:
	rm -rf *~ *.exec
