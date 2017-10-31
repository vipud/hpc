seq_bb_predict = open("seq_bb_predict.dat", "r")
seq_bmrb_pre = open("seq_bmrb_pre.dat", "r")
seq_cs_rmsd = open("seq_cs_rmsd.dat", "r")
seq_proton_predict = open("seq_proton_predict.dat", "r")

para_bb_predict = open("bb_predict.dat", "r")
para_bmrb_pre = open("bmrb_pre.dat", "r")
para_cs_rmsd = open("cs_rmsd.dat", "r")
para_proton_predict = open("proton_predict.dat", "r")

if seq_bb_predict.read() != para_bb_predict.read():
	print("bb_predict incorrect")
else:
	print("bb_predict correct")

if seq_bmrb_pre.read() != para_bmrb_pre.read():
	print("bmrb_pre incorrect")
else:
	print("bmrb_pre correct")

if seq_cs_rmsd.read() != para_cs_rmsd.read():
	print("cs_rmsd incorrect")
else:
	print("cs_rmsd correct")

if seq_proton_predict.read() != para_proton_predict.read():
	print("proton_predict incorrect")
else:
	print("proton_predict correct")


