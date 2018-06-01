import sys

dirname=sys.argv[1]

output = open("output.txt","w+")

print("bb_predict.dat")
output.write("bb_predict.dat")

filename="bb_predict.dat"
file1=open(filename,"r")
file2=open(dirname + "/" + filename, "r")

line1 = file1.readline()
line2 = file2.readline()

while line1 != '' and line2 != '':
	if line1 != line2:
		print(line1)
		output.write(line1)
		print(line2)
		output.write(line2)
		print("\n")
		output.write("\n")
	line1 = file1.readline()
	line2 = file2.readline()

print("Finished")
output.write("Finished")

print("proton_predict.dat")
output.write("proton_predict.dat")

filename="proton_predict.dat"
file1=open(filename,"r")
file2=open(dirname + "/" + filename, "r")

line1 = file1.readline()
line2 = file2.readline()

while line1 != '' and line2 != '':
	if line1 != line2:
		print(line1)
		output.write(line1)
		print(line2)
		output.write(line2)
		print("\n")
		output.write("\n")
	line1 = file1.readline()
	line2 = file2.readline()

print("Finished")
output.write("Finished")

print("cs_rmsd.dat")
output.write("cs_rmsd.dat")

filename="cs_rmsd.dat"
file1=open(filename,"r")
file2=open(dirname + "/" + filename, "r")

line1 = file1.readline()
line2 = file2.readline()

while line1 != '' and line2 != '':
	if line1 != line2:
		print(line1)
		output.write(line1)
		print(line2)
		output.write(line2)
		print("\n")
		output.write("\n")
	line1 = file1.readline()
	line2 = file2.readline()

print("Finished")
output.write("Finished")

print("bmrb_pre.dat")
output.write("bmrb_pre.dat")

filename="bmrb_pre.dat"
file1=open(filename,"r")
file2=open(dirname + "/" + filename, "r")

line1 = file1.readline()
line2 = file2.readline()

while line1 != '' and line2 != '':
	if line1 != line2:
		print(line1)
		output.write(line1)
		print(line2)
		output.write(line2)
		print("\n")
		output.write("\n")
	line1 = file1.readline()
	line2 = file2.readline()

print("Finished")
output.write("Finished")

