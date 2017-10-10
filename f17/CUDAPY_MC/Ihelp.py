import csv

with open('input.csv', mode='r') as infile:
    reader = csv.reader(infile)
    input_dict = {rows[0]:rows[1] for rows in reader}

SP = input_dict['StockPrice']
K = input_dict['StrikePrice']
Beta = input_dict['Beta']
IR = input_dict['InterestRate']
Maturity = float(input_dict['Maturity'])
Number_of_Paths = input_dict['Paths']
Number_of_Steps = input_dict['Steps']
