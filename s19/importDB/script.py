import os
import csv
import dns

from pymongo import MongoClient
client = MongoClient('mongodb+srv://user1:1234@cluster0-qr76p.mongodb.net/test?retryWrites=true')
db = client['hpc']

#path = "/Users/max/Downloads/output/"
path = "/hpc/s19/webscraper/csvscraper/FileFolder/raw-csvs-openacc1.2/output"

collection = db.test

for filename in os.listdir(path):
    with open(path + filename, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            if collection.count({'CPU_name': row[1]}) == 0:
                collection.insert_one({'Score': [row[0]],
                                       'CPU_name': row[1],
                                       'CPU(s)_enabled': row[2],
                                       'ACCEL_Model_Name': row[3],
                                       'Memory': row[4],
                                       'OS': row[5]})
            else:
                collection.update_one({'CPU_name': row[1]}, {'$push': {'Score': row[0]}})

for c in collection.find():
    print(c)
