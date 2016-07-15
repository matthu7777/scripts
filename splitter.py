#!/usr/bin/env python

import mgutils as mg
import numpy as np
import matplotlib.pyplot as plt
from sys import argv
from os import remove
from os.path import isfile



fname = argv[1]
datname = argv[2]
tempnames = ["a.mod", "b.mod", "c.mod"]

names = ["t1", "temp_disc", "temp_spot"]
reps = [15, 1, 1]





for i in range(len(names)):
	labels = np.loadtxt(fname, usecols=(0,), dtype=str)
	col1 = np.genfromtxt(fname, usecols=(1,), dtype=str, skip_footer=34)
	vals = np.loadtxt(fname, usecols=(2,), dtype=str)
	col3 = np.genfromtxt(fname, usecols=(3,), dtype=str, skip_footer=34)
	col4 = np.genfromtxt(fname, usecols=(4,), dtype=str, skip_footer=34)
	col5 = np.genfromtxt(fname, usecols=(5,), dtype=str, skip_footer=34)
	col6 = np.genfromtxt(fname, usecols=(6,), dtype=str, skip_footer=34)


	
	tempname = tempnames[i]
	
	
	for j in range(len(names)):
		if j != i:
			vals[labels==names[j]] = reps[j]
			print labels[labels==names[j]], reps[j]
	
	if isfile(tempname):
		remove(tempname)
	
	temp = open(tempname, 'a')
	
	for a,b,c,d,e,f,g in zip(labels[:len(col1)], col1, vals[:len(col1)], col3, col4, col5, col6):
		temp.write(a+" "+b+" "+c+" "+d+" "+e+" "+f+" "+g+"\n")
	
	for a,b,c in zip(labels[len(col1):], col1[:34], vals[len(col1):]):
		temp.write(a+" "+b+" "+c+"\n")
	
	print "seperated", names[i]
	temp.close()
	
	





