CC = mpicc
EXEC = mpiexec -n 4
CFLAGS = -std=c99

all : hello.exec

hello.exec : hello.c
	$(CC) $(CFLAGS) hello.c -o hello.exec
run_hello : hello.exec
	$(EXEC) hello.exec
clean:
	rm -rf *~ *.exec *.out *.prof
