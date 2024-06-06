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



N = 104*1579
a = 1014

x0 = 1./N
listx =  [x0]
for idx in range(12):
	xk = listx[-1]
	xk1 = 2*xk - a*xk*xk
	listx.append(xk1)
	print xk

print listx
print 'a = ',a
print '1/a = ',1./a
plt.plot(range(len(listx)),        listx,          's-k')
#plt.axis("equal") 
plt.show()
#plt.savefig("TIME4.png")
#plt.close()

print log(0.1)/log(0.99)