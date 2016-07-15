#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from trm import dnl, subs, sla, roche
from scipy.optimize import curve_fit
from sys import argv, exit
import stars, runs
from subprocess import call
from os import listdir
import mgutils as mg

""" Script to compare the results of multiple mcmc runs on different data sets. 
"""

files = "outfiles4.lis"

plotHists = False	#Plot histograms of individual mcmcs?
testData = False		#Test normality of data distribution?

if len(argv) > 1:
	par = argv[1]
else:
	par = "wdwarf"

def readMcmcLog(fname, par):
	""" Read data points from mcmc .log file and return those matching the requested parameter
	"""
	# find params that varied
	for line in open(fname):
		if line.startswith('##'):
			labels = line.split()[1:]
			labels = np.array(labels, dtype=str)
			break
	#read data and select the desired parameter
	data = np.loadtxt(fname)
	oi = data[:, (labels==par)]
	if len(oi) == 0:
		print "!!ERROR finding that parameter in log file!!"
	return oi

def mcmcHist(data, plotHists=True, label=''):
	"""	Returns the mean and std of a dataset, as well correlation with a Gaussian distribution. Can optionally plot a histogram
	"""
	## Get histogram of results
	n, bins, patches = plt.hist(data, bins=20)
	
	## Fit with Gaussian
	try:
		mids = (bins[1:] + bins[:-1]) / 2
		a0 = data.max()
		m0 = data.mean()
		s0 = data.std()
		popt, pcov = curve_fit(mg.gauss, mids, n, p0=[a0,m0,s0])
		
		a, m, s = popt
		
		rms = np.sqrt( ((n - mg.gauss(mids, *popt))**2).mean() )
		
		print label, "RMS from Gaussian distribution", rms
		
		plt.plot(mids, mg.gauss(mids, *popt), 'r')
	except RuntimeError:
		print label, "did not converge"
	
	if plotHists is not True:
		plt.clf()
	else:
		plt.show()
	return data.mean(), data.std(), rms
	
	

################# Read in data and find means
names = []
ms, ss, ecls = np.zeros(75), np.zeros(75), np.zeros(75)
r,g,b = np.zeros(75, dtype=bool), np.zeros(75, dtype=bool), np.zeros(75, dtype=bool)
with open(files, 'r') as lis:
	for j, fname in enumerate(lis):
		fname = fname[:-1]		#trim the newline character
		name = fname.split('_')[0]
		data = readMcmcLog(fname, par)
		
		## If desired, report goodness of fit of data to a Gaussian
		if testData is True:
			m, s, rms = mcmcHist(data, plotHists=plotHists,label=name)
		else:
			m, s = data.mean(), data.std()
		eclnum = runs.eclnums[name]
		
		ecls[j] = eclnum
		ms[j] = m
		ss[j] = s
		
		if "r" in name:
			r[j] = True
		if "g" in name:
			g[j] = True
		if "b" in name:
			b[j] = True

############# Analysis

avg = np.average( ms[g], weights=((1/ss[g])**2) )
print "Average", avg
chisq = mg.chisq(avg, ms[g], ss[g])
print "Chis sqared", chisq

iteravg, pcov, mask = mg.iterativeFit(mg.flat, np.arange(len(ms[g])), ms[g], ss[g], 10, 3, p0=[avg])
print "Iterative average", iteravg[0], "masked", len(mask[mask==False])


## Individual nights? Jan / May?





################ Plotting


## Make array of x axis tick labels
labels = []
for run in runs.runoblist:		#for nights
	label = run.name.title()
	labels += [label]
#for run in runs.runoblist:			#for individual eclipses
	#for j in range(1,run.numec+1):
		#label = run.name.title()+" eclipse %d"%(j)
		#labels += [label]

## Format graph
mg.formatGraph(ylabel="Flux (mJy)")
#plt.xticks(np.arange(1,25), labels, rotation=17)
plt.xticks([2,5.5,10,14,18.5,23.5],labels, rotation=17)
plt.gca().set_xlim(0,26)
plt.gca().set_ylim(0.04, 0.18)

## Plot data
plt.errorbar(ecls[r], ms[r], ss[r], fmt='r.')
plt.errorbar(ecls[g], ms[g], ss[g], fmt='g.')
plt.errorbar(ecls[b], ms[b], ss[b], fmt='b.')

## Lines to separate days
plt.vlines([3.5,7.5,12.5,15.5,21.5], 0.08, 0.16, linestyle='dashed')

plt.show()






