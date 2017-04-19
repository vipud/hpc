CC = mpicc 
EXEC = srun -n 1 --gres=gpu:1
FLAGS =  -fast -ta=nvidia:cc60,fastmath -acc -Minfo=ccff  

all: hello_add.exec 

hello_seq.exec: hello_seq.c
	$(CC) $(FLAGS) hello_seq.c -o hello_seq.exec

hello_acc.exec: hello_acc.c
	$(CC) $(FLAGS) hello_acc.c -o hello_acc.exec

hello_add.exec : hello_add.c
	$(CC) $(FLAGS) hello_add.c -o hello_add.exec

hello.exec : hello.c
	$(CC) hello.c -o hello.exec

run: hello_add.exec hello_seq.exec
	$(EXEC) hello_add.exec 100000
	mpiexec -n 1 hello_seq.exec 100000
clean:
	rm -rf *~ *.o *.exec
