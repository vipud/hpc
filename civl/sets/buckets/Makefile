# GENERIC MAKEFILE
# CHANGE/DELETE AS NEEDED

CC = gcc
CIVL = civl

# compile - 'make linked_list'
# run - './linked_list'
# verify - 'make verify'
# clean extra files - 'make clean'

linked_list : linked_list.c
	$(CC) linked_list.c linked_list.h 
	./a.out
buckets: buckets.c
	$(CC) buckets.c linked_list.h
	./a.out

verify : linked_list.c
	$(CIVL) verify linked_list.c
clean:
	rm -rf *~ *.out *.trace *.log *.gch
	rm -rf example linked_list
