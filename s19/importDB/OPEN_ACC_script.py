
from pymongo import MongoClient
import os
import csv
client = MongoClient()

db = client['hpc']

path = '/usa/dgbaum/vipgithub/hpc/s19/webscraper/data/openacc/raw-csvs-openacc1.2/output-current/'

collection = db.openAcc
for filename in os.listdir(path):
	with open(path + filename, 'r') as f:
		reader = csv.reader(f)
		for row in reader:
			collection.insert_one({
			'_id': row[9],
			'Hardware_Model': row[0],
			'SPECaccel_acc_base': row[1],
			'CPU Name': row[2],
			'CPU(s)_Enabled': row[3],
			'Accel Model Name': row[4],
			'Memory': row[5],
            'Operating System': row[6],
            'Other Software': row[7],
            'Compiler': row[8],
			})

