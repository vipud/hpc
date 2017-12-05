# STRIPED HASH MAKEFILE

CC = gcc
CIVL = civl

# replace 'example.c' with your file name!
# compile - 'make example'
# run - './example'
# verify - 'make verify'
# clean extra files - 'make clean'

striped : striped_hash.cvl
	$(CC) striped_hash.cvl -o striped_hash

verify : striped_hash.cvl
	$(CIVL) verify striped_hash.cvl

clean:
	rm -rf *~ *.out *.trace *.log *.gch