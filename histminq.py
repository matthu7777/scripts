#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from trm import dnl, subs, sla, roche
import stars, runs
import mgutils as mg

""" Script to plot histogram of min q for roche-filling stars. Reads output of format put out by readqi.py
"""

fname = "qi.dat"

names = np.loadtxt(fname, usecols=(0,), dtype=str)
data = np.loadtxt(fname, usecols=(1,2,3,4,5,6,7,8,))

q = data[:,0]
qerr = data[:,1]
i = data[:,2]
ierr = data[:,3]
phi = data[:,4]
phierr = data[:,5]
minq = data[:,6]
minqerr = data[:,7]

nz = (q != 0)
red = np.zeros(len(names[nz]), dtype=bool)
gre = np.zeros(len(names[nz]), dtype=bool)
blu = np.zeros(len(names[nz]), dtype=bool)
for name in names[nz]:
	if "-r" in name:
		red[names[nz]==name] = True
	if "-g" in name:
		gre[names[nz]==name] = True
	if "-b" in name:
		blu[names[nz]==name] = True

## Format graph
mg.formatGraph("Occurrence", "Minimum mass ratio $q$")


## Plot histogram
plt.hist(minq[nz], bins=20)


print names
print minq
print names[nz][minq[nz] < 0.016]
print names[minq < 0.016]

plt.show()
