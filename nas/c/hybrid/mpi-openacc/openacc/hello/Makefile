CC = mpicc 
EXEC = mpiexec
FLAGS =  -fast -acc -ta=nvidia:cc2,fastmath -Minfo=accel -Wl, -rpath,/opt/pgi/modulefiles/openmpi/1.10.2/2016 

all: hello_add.exec 

hello_add.exec : hello_add.c
	$(CC) $(FLAGS) hello_add.c -o hello_add.exec

hello.exec : hello.c
	$(CC) hello.c -o hello.exec

run: hello_add.exec
	$(EXEC) -n 1  hello_add.exec
clean:
	rm -rf *~ *.o *.exec
