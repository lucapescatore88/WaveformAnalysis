import csv
import os


data_dir = os.path.dirname("/home/analysis/Analysis_waveforms/H2015/H2015_sequence/55.0V/")

conv_dir = os.path.dirname("/home/analysis/Analysis_waveforms/H2015/H2015_sequence/55.0V/Conv/")
with open(data_dir+"/C1H00000.csv",'rb') as csvfile:
	
	spamreader = csv.reader(csvfile)
	for i in range(0,24):
		next(spamreader, None)
	
	with open(conv_dir+"/C1H00000.csv","wb") as f:
		for row in spamreader:
			writer = csv.writer(f)
			writer.writerows(range(25,2002))	
		



	 
