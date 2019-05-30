
from pymongo import MongoClient
import os
import csv
client = MongoClient()

db = client['hpc']

path = '/usa/dgbaum/vipgithub/hpc/s19/webscraper/data/mpi2007/output/'

collection = db.mpi2007
for filename in os.listdir(path):
	with open(path + filename, 'r') as f:
		reader = csv.reader(f)
		for row in reader:
			collection.insert_one({
			'_id': row[9],
			'Hardware_Model': row[0],
			'SPECmpiL_base2007:': row[1],
			'Chips enabled': row[2],
			'Cores enabled': row[3],
			'Cores per chip': row[4],
			'Threads per core': row[5],
            'Memory': row[6],
            'Operating System': row[7],
            'Other Software': row[8],
			})
