# GENERIC MAKEFILE
# CHANGE/DELETE AS NEEDED

CC = gcc
CIVL = civl

# replace 'example.c' with your file name!
# compile - 'make example'
# run - './example'
# verify - 'make verify'
# clean extra files - 'make clean'

example : BaseHash.c
	$(CC) BaseHash.c -o BaseHash

verify : BaseHash.c
	$(CIVL) verify BaseHash.c
clean:
	rm -rf *~ *.out *.trace *.log *.gch BaseHash
