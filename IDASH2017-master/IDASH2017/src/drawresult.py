# coding=utf8
# 2019-05-14 01:17 a.m. GMT +08：00
'''
+++++++++++++++++++++++++++++++++++++++++++++++++++++
+       NumPy version 1.10.2    Python 2.7.11       +
+++++++++++++++++++++++++++++++++++++++++++++++++++++
'''
import os
import math
import time
import random

from copy import deepcopy
from math import log, exp, pow, sqrt

import matplotlib
#matplotlib.use('pdf')
import matplotlib.pyplot as plt
import numpy as np

import csv



pathNesterov                 = '../data/testPlainResultNesterov_'
pathNesterovWithG            = '../data/testPlainResultNesterovWithG_'
pathBonteSFH                 = '../data/testPlainResultBonteSFH_'
pathBonteSFHwithLearningRate = '../data/testPlainResultBonteSFHwithLearningRate_'

whichtype                        = 'MLE.csv' # 'AUC.csv', 'MLE.csv', 'TIME.csv'

#     Step 1. Extract data from a csv file
with open(pathNesterov + whichtype,'r') as csvfile:
	reader = csv.reader(csvfile)
	#reader.next() # leave behind the first row
	data = []
	for row in reader:
		# reader.next() return a string
		row = [float(x) for x in row[1:]]
		data.append(row)
csvfile.close()
if whichtype == 'TIME.csv' :
	for rowidx in range(len(data)):
		pretime = 0.0
		for idx in range(len(data[0])):
			data[rowidx][idx] = pretime + data[rowidx][idx]
			pretime = data[rowidx][idx]
if whichtype == 'AUC.csv' :
	for rowidx in range(len(data)):
		data[rowidx][0] = data[rowidx][1]

Result_Nesterov = [0.0 for i in data[0] ]
for rowidx in range(len(data)):
	for idx in range(len(data[0])):
		Result_Nesterov[idx] += data[rowidx][idx]
for idx in range(len(data[0])):
	Result_Nesterov[idx] = Result_Nesterov[idx]/len(data)



#     Step 1. Extract data from a csv file
with open(pathNesterovWithG + whichtype,'r') as csvfile:
	reader = csv.reader(csvfile)
	#reader.next() # leave behind the first row
	data = []
	for row in reader:
		# reader.next() return a string
		row = [float(x) for x in row[1:]]
		data.append(row)
csvfile.close()
if whichtype == 'TIME.csv' :
	for rowidx in range(len(data)):
		pretime = 0.0
		for idx in range(len(data[0])):
			data[rowidx][idx] = pretime + data[rowidx][idx]
			pretime = data[rowidx][idx]
if whichtype == 'AUC.csv' :
	for rowidx in range(len(data)):
		data[rowidx][0] = data[rowidx][1]

Result_NesterovWithG = [0.0 for i in data[0] ]
for rowidx in range(len(data)):
	for idx in range(len(data[0])):
		Result_NesterovWithG[idx] += data[rowidx][idx]
for idx in range(len(data[0])):
	Result_NesterovWithG[idx] = Result_NesterovWithG[idx]/len(data)



#     Step 1. Extract data from a csv file
with open(pathBonteSFH + whichtype,'r') as csvfile:
	reader = csv.reader(csvfile)
	#reader.next() # leave behind the first row
	data = []
	for row in reader:
		# reader.next() return a string
		row = [float(x) for x in row[1:]]
		data.append(row)
csvfile.close()
if whichtype == 'TIME.csv' :
	for rowidx in range(len(data)):
		pretime = 0.0
		for idx in range(len(data[0])):
			data[rowidx][idx] = pretime + data[rowidx][idx]
			pretime = data[rowidx][idx]
if whichtype == 'AUC.csv' :
	for rowidx in range(len(data)):
		data[rowidx][0] = data[rowidx][1]

Result_BonteSFH = [0.0 for i in data[0] ]
for rowidx in range(len(data)):
	for idx in range(len(data[0])):
		Result_BonteSFH[idx] += data[rowidx][idx]
for idx in range(len(data[0])):
	Result_BonteSFH[idx] = Result_BonteSFH[idx]/len(data)



#     Step 1. Extract data from a csv file
with open(pathBonteSFHwithLearningRate + whichtype,'r') as csvfile:
	reader = csv.reader(csvfile)
	#reader.next() # leave behind the first row
	data = []
	for row in reader:
		# reader.next() return a string
		row = [float(x) for x in row[1:]]
		data.append(row)
csvfile.close()
if whichtype == 'TIME.csv' :
	for rowidx in range(len(data)):
		pretime = 0.0
		for idx in range(len(data[0])):
			data[rowidx][idx] = pretime + data[rowidx][idx]
			pretime = data[rowidx][idx]
if whichtype == 'AUC.csv' :
	for rowidx in range(len(data)):
		data[rowidx][0] = data[rowidx][1]

Result_BonteSFHwithRate = [0.0 for i in data[0] ]
for rowidx in range(len(data)):
	for idx in range(len(data[0])):
		Result_BonteSFHwithRate[idx] += data[rowidx][idx]
for idx in range(len(data[0])):
	Result_BonteSFHwithRate[idx] = Result_BonteSFHwithRate[idx]/len(data)


#label = ["Bonte's SFH",'Nesterov', 'Method1 (SFH with rate)', 'Method4 (Nesterov +G by .25XTX)']
label = ['Nesterov','Nesterov +G by .25XTX',"Bonte's SFH","Bonte's SFH with Rate"]
plt.plot(range(len(Result_Nesterov)),        Result_Nesterov,          's-k')
plt.plot(range(len(Result_NesterovWithG)),   Result_NesterovWithG,     'v-b')
plt.plot(range(len(Result_BonteSFH)),        Result_BonteSFH,          'p-g')
plt.plot(range(len(Result_BonteSFHwithRate)),Result_BonteSFHwithRate,  'o-c')
#plt.axis("equal")
plt.title(whichtype)
plt.legend(label, loc = 0, ncol = 1)  
plt.show()
#plt.savefig("TIME4.png")
#plt.close()





#     Step 1. Extract data from a csv file
with open(pathNesterovWithG + whichtype,'r') as csvfile:
	reader = csv.reader(csvfile)
	#reader.next() # leave behind the first row
	data = []
	for row in reader:
		# reader.next() return a string
		row = [float(x) for x in row[1:]]
		data.append(row)
csvfile.close()
