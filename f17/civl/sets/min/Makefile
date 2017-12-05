# Settings for Makefile
# Path to your CIVL JAR file.
CIVL_BIN = /Users/wwh/Documents/workspace/CIVL/civl.jar
# Path to your set implementation
IMPL = min_set.cvl

# ==== ==== ==== ==== Change constants above ==== ==== ==== ====

# Command elements
JVM = java -jar 
VERF = verify
TRANS = -showTransitions

# Files required
HEADER = set.cvh
DRIVER = min_driver.cvl
UTILS = num_perms.c

verify: $(UTILS) $(HEADER) $(IMPL) $(DRIVER)
	$(JVM) $(CIVL_BIN) $(VERF) $(DRIVER)

show: $(UTILS) $(HEADER) $(IMPL) $(DRIVER)
	$(JVM) $(CIVL_BIN) $(VERF) $(TRANS) $(DRIVER)

clean:
	rm -rf *~ *.out *.trace *.log *.gch