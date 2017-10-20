# GENERIC MAKEFILE
# CHANGE/DELETE AS NEEDED

CC = gcc
CIVL = civl

# replace 'example.c' with your file name!
# compile - 'make example'
# run - './example'
# verify - 'make verify'
# clean extra files - 'make clean'

example : example.c
	$(CC) example.c -o example

verify : example.c
	$(CIVL) verify example.c
clean:
	rm -rf *~ *.out *.trace *.log *.gch
