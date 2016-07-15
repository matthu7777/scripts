#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from trm import dnl, subs, sla
from scipy.optimize import curve_fit
from sys import argv, exit
import stars, runs
import mgutils as mg

""" Plots a histogram of results from a .log file. Usage: hist.py logfile.log (param) (cburn). Burns points where chisq < cburn
"""


if len(argv) > 1:
	fname = argv[1]
else:
	fname = "r4.log"

if len(argv) > 2:
	param = argv[2]
else:
	param = None

if len(argv) > 3:
	cburn = float(argv[3])
else:
	cburn = None


for line in open(fname):
	if line.startswith('##'):
		labels = np.array(line.split()[1:])

data = np.loadtxt(fname)

chisq = data[:,labels=='chisq'].T[0]
pprob = data[:,labels=='pprob'].T[0]
lp = data[:,labels=='lp'].T[0]
if cburn is not None:
	ok = chisq < cburn
else:
	ok = np.ones(len(chisq),dtype=bool)
	
print "Using", len(ok[ok==True]), "out of", len(ok)

for i in range(len(labels)-2):
	if param==None or param==labels[i]:
		oi = data[:,i]
		print len(oi)
		oi = oi[ok]
		plt.figure(i)
		n, bins, patches = plt.hist(oi, bins=20)
		plt.suptitle(labels[i]+"    ("+fname+")")
		
		print labels[i], "Mean", oi.mean(), "Std", oi.std()
		
		try:
			mids = (bins[1:] + bins[:-1]) / 2
			m0 = oi.mean()
			popt, pcov = curve_fit(mg.gauss, mids, n, p0=[0,m0,1])
			
			a, m, s = popt
			
			plt.plot(mids, mg.gauss(mids, *popt), 'r')
			
			print labels[i], "Gaussian fit", m, "Sigma", s
		except RuntimeError:
			print labels[i], "did not converge"


mg.saveAsTemp()		
plt.show()
	
