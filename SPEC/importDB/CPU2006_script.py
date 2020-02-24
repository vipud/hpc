from pymongo import MongoClient
import os
import csv
client = MongoClient()

db = client['hpc']

path = '/usa/dgbaum/vipgithub/hpc/s19/webscraper/data/cpu2006/output/'

collection = db.cpu2006
for filename in os.listdir(path):
	with open(path + filename, 'r') as f:
		reader = csv.reader(f)
		for row in reader:
			collection.insert_one({
			'_id': row[8],
			'Hardware_Model': row[0],
			'SPECint_base2006:': row[1],
			'CPU_Name': row[2],
			'CPU(s)_Enabled': row[3],
			'Memory': row[4],
			'Operating System': row[5],
			'Other Software': row[6],
            'Compiler': row[7],
			})
