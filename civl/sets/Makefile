# make <algorithmname_type>
# For running test driver with your algorithm interfaced

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

cuckoo_2 : cuckoo/cuckoo_seq.cvl driver/max_driver.cvl driver/perm_util.c \
          include/set.cvh
	$(CIVL) verify -inputVAL_B=2 -inputTHD_B=2 -inputSTP_B=2 \
                driver/max_driver.cvl cuckoo/cuckoo_seq.cvl

striped_cuckoo_2 : cuckoo/striped_cuckoo.cvl driver/max_driver.cvl \
	          driver/perm_util.c include/set.cvh
		$(CIVL) verify -inputVAL_B=3 -inputTHD_B=2 -inputSTP_B=2 \
                driver/max_driver.cvl cuckoo/striped_cuckoo.cvl

coarse_list_2 : coarselist/coarse_list.cvl driver/list_max_driver.cvl \
		driver/perm_util.c include/set.cvh
		$(CIVL) verify -inputVAL_B=3 -inputTHD_B=2 -inputSTP_B=2 \
                driver/list_max_driver.cvl coarselist/coarse_list.cvl

naive_6v : naive/naive_set.cvl driver/max_driver.cvl driver/perm_util.c \
          include/set.cvh
	$(CIVL) verify -inputVAL_B=6 -inputTHD_B=2 -inputSTP_B=2 \
                driver/max_driver.cvl naive/naive_set.cvl 

striped_6v : striped/striped_set.cvl driver/max_driver.cvl driver/perm_util.c \
          include/set.cvh
	$(CIVL) verify -inputVAL_B=6 -inputTHD_B=2 -inputSTP_B=2 \
                driver/max_driver.cvl striped/striped_set.cvl

coarse_6v : coarse/coarse_set.cvl driver/max_driver.cvl driver/perm_util.c \
          include/set.cvh
	$(CIVL) verify -inputVAL_B=6 -inputTHD_B=2 -inputSTP_B=2 \
                driver/max_driver.cvl coarse/coarse_set.cvl

cuckoo_6v : cuckoo/cuckoo_seq.cvl driver/max_driver.cvl driver/perm_util.c \
          include/set.cvh
	$(CIVL) verify -inputVAL_B=6 -inputTHD_B=2 -inputSTP_B=2 \
                driver/max_driver.cvl cuckoo/cuckoo_seq.cvl

naive_3t : naive/naive_set.cvl driver/max_driver.cvl driver/perm_util.c \
          include/set.cvh
	$(CIVL) verify -inputVAL_B=2 -inputTHD_B=3 -inputSTP_B=2 \
                driver/max_driver.cvl naive/naive_set.cvl 

striped_3t : striped/striped_set.cvl driver/max_driver.cvl driver/perm_util.c \
          include/set.cvh
	$(CIVL) verify -inputVAL_B=2 -inputTHD_B=3 -inputSTP_B=2 \
                driver/max_driver.cvl striped/striped_set.cvl

coarse_3t : coarse/coarse_set.cvl driver/max_driver.cvl driver/perm_util.c \
          include/set.cvh
	$(CIVL) verify -inputVAL_B=2 -inputTHD_B=3 -inputSTP_B=2 \
                driver/max_driver.cvl coarse/coarse_set.cvl

cuckoo_3t : cuckoo/cuckoo_seq.cvl driver/max_driver.cvl driver/perm_util.c \
          include/set.cvh
	$(CIVL) verify -inputVAL_B=2 -inputTHD_B=3 -inputSTP_B=2 \
                driver/max_driver.cvl cuckoo/cuckoo_seq.cvl

naive_3s : naive/naive_set.cvl driver/max_driver.cvl driver/perm_util.c \
          include/set.cvh
	$(CIVL) verify -inputVAL_B=2 -inputTHD_B=2 -inputSTP_B=3 \
                driver/max_driver.cvl naive/naive_set.cvl 

striped_3s : striped/striped_set.cvl driver/max_driver.cvl driver/perm_util.c \
          include/set.cvh
	$(CIVL) verify -inputVAL_B=2 -inputTHD_B=2 -inputSTP_B=3 \
                driver/max_driver.cvl striped/striped_set.cvl

coarse_3s : coarse/coarse_set.cvl driver/max_driver.cvl driver/perm_util.c \
          include/set.cvh
	$(CIVL) verify -inputVAL_B=2 -inputTHD_B=2 -inputSTP_B=3 \
                driver/max_driver.cvl coarse/coarse_set.cvl

cuckoo_3s : cuckoo/cuckoo_seq.cvl driver/max_driver.cvl driver/perm_util.c \
          include/set.cvh
	$(CIVL) verify -inputVAL_B=2 -inputTHD_B=2 -inputSTP_B=3 \
                driver/max_driver.cvl cuckoo/cuckoo_seq.cvl


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
