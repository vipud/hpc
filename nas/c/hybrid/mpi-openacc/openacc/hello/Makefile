CC = pgcc 
LIBS = 
CFLAGS = -Mmpi=mpich -Minfo=all -Mprof=ccff -ta=tesla:cc60

all : hello.exec hello_add.exec

hello_add.exec : hello_add.c
	$(CC) $(CFLAGS) hello_add.c -o hello_add.exec

hello.exec : hello.c
	$(CC) $(CFLAGS) hello.c -o hello.exec

clean:
	rm -rf *~ *.exec
