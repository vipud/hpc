# GENERIC MAKEFILE
# CHANGE/DELETE AS NEEDED

CC = gcc
CIVL = civl

naive_2 : naive/naive_set.cvl driver/max_driver.cvl driver/perm_util.c \
          include/set.cvh
	$(CIVL) verify -inputVAL_B=2 -inputTHD_B=2 -inputSTP_B=2 \
                driver/max_driver.cvl naive/naive_set.cvl 

striped_2 : striped/striped_set.cvl driver/max_driver.cvl driver/perm_util.c \
          include/set.cvh
	$(CIVL) verify -inputVAL_B=2 -inputTHD_B=2 -inputSTP_B=2 \
                driver/max_driver.cvl striped/striped_set.cvl

coarse_2 : coarse/coarse_set.cvl driver/max_driver.cvl driver/perm_util.c \
          include/set.cvh
	$(CIVL) verify -inputVAL_B=2 -inputTHD_B=2 -inputSTP_B=2 \
                driver/max_driver.cvl coarse/coarse_set.cvl

cuckoo_2 : cuckoo/cuckoo_set.cvl driver/max_driver.cvl driver/perm_util.c \
          include/set.cvh
	$(CIVL) verify -inputVAL_B=2 -inputTHD_B=2 -inputSTP_B=2 \
                driver/max_driver.cvl cuckoo/cuckoo_set.cvl

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
