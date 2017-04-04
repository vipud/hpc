CC = mpicc
FLAGS = -ta=tesla:cc60 

all : loop loop_acc

loop : loop.c
	$(CC) loop.c -o loop

loop_acc : loop_acc.c
	$(CC) $(FLAGS) loop_acc.c -o loop_acc

clean:
	rm -rf loop loop_acc *~ *.out
