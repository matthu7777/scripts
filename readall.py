#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from trm import dnl, subs, sla
from scipy.optimize import curve_fit
from sys import argv, exit
import stars, runs
from subprocess import call
from os import listdir
import mgutils as mg

par = argv[1]

num = "1"

fnames = []


skip = ["a.mod","b.mod","c.mod","d.mod"]

for f in listdir("."):
	if f[-4:] == ".mod" and f not in skip and "gaia" not in f and "_%s"%(num) in f:
		fnames += [f]
		
# find params that varied
for line in open(fnames[0][:-4]+".log"):
	if line.startswith('##'):
		oi = line.split()[1:]

params = np.zeros((90,len(fnames)))

# extract best models and read 
for i in range(len(fnames)):
	fname = fnames[i]
	
	labels = np.loadtxt(fname, usecols=(0,), dtype=str)
	vals = np.loadtxt(fname, usecols=(2,), dtype=str)
	
	labels = labels[(labels!="limb1")&(labels!="limb2")]
	vals = vals[vals!="Poly"].astype(float)
	
	params[:,i] = vals



		

colours = np.zeros(len(fnames), dtype=str)

for i in range(len(fnames)):
	f = fnames[i]
	if "-r-" in f:
		colours[i] = "r"
	if "-g-" in f:
		colours[i] = "g"
	if "-b-" in f:
		colours[i] = "b"

eclnum = np.zeros(len(fnames))
for i in range(len(fnames)):
	eclnum[i] = runs.eclnums[fnames[i][:-6]]


for c in ["r", "g", "b"]:
	y = params[(labels==par)][0][colours==c]
	x = eclnum[colours==c]
	
	plt.plot(x, y, (c+"."))

plt.vlines([3.5,7.5,12.5,15.5,21.5],params[(labels==par)][0].min(),params[(labels==par)][0].max())

mg.formatGraph(xlabel="Eclipse Number",ylabel=par)
mg.saveAsTemp()


plt.show()










