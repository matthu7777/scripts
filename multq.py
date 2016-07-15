#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from sys import argv, exit
from subprocess import call
from os import listdir
from scipy.optimize import curve_fit
import runs
import mgutils as mg

""" Script to compare the results of multiple mcmc runs on different data sets, tailored to eclipse depths of gaia14aae but can plot other parameters too
"""

files = "outfiles4.lis"
datname = "meanq.dat"

refreshData = False		#Re-import data from log files? (else reads means from a text file)
testData = False		#Test normality of data distribution?
plotHists = False	#Plot histograms of individual mcmcs?
saveData = True

if len(argv) > 1:
	par = argv[1]
else:
	par = "q"

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
	

def chisq(model, data, err):
	return np.sum( (data-model)**2 / err**2 )



################# Read in data and find means
names = []
ms, ss, ecls = np.zeros(75), np.zeros(75), np.zeros(75)
i,r,g,b = np.zeros(75, dtype=bool), np.zeros(75, dtype=bool), np.zeros(75, dtype=bool), np.zeros(75, dtype=bool)
with open(files, 'r') as lis:
	for j, fname in enumerate(lis):
		fname = fname[:-1]		#trim the newline character
		name = fname.split('_')[0]
		eclnum = runs.eclnums[name]
		
		if refreshData:
			## read in data and find means
			data = readMcmcLog(fname, par)
			m, s = data.mean(), data.std()
			
			print m 
			
			ecls[j] = eclnum
			ms[j] = m
			ss[j] = s

			## If desired, report goodness of fit of data to a Gaussian
			if testData is True:
				m, s, rms = mcmcHist(data, plotHists=plotHists,label=name)
			else:
				m, s = data.mean(), data.std()
		
		
		
		if "r" in name and eclnum <= 3:
			i[j] = True
		if "r" in name and eclnum >= 4:
			r[j] = True
		if "g" in name:
			g[j] = True
		if "b" in name:
			b[j] = True
if not refreshData:
	## load means from file
	data = np.loadtxt(datname)
	ecls = data[:,0]
	ms = data[:,1]
	ss = data[:,2]
elif saveData:
	## Save means to dat file
	tosave = np.column_stack((ecls, ms, ss, i.astype(int), r.astype(int), g.astype(int), b.astype(int)))
	np.savetxt(datname, tosave)
	print "Saved to", datname



############# Analysis
colours = ['m','r','g','b']
for j, c in enumerate([i,r,g,b]):
	col = colours[j]

	avg = np.average( ms[c], weights=((1/ss[c])**2) )
	av,averr = mg.weightedAv(ms[c],ss[c])
	avgchisq = chisq(avg, ms[c], ss[c])
	print "\nAverage %s"%(col), avg, "Error", averr, "Chi squared", avgchisq
	#plt.hlines(avg,0,26,linestyle='dashed', color=col)
	
	
	iteravg, pcov, mask = mg.iterativeFit(mg.flat, np.arange(len(ms[c])), ms[c], ss[c], 10, 3, p0=[avg])
	print "Iterative average", iteravg[0], "masked", len(mask[mask==False])

## Find average for each eclipse
avs, averrs = np.zeros(len(ecls[g])), np.zeros(len(ecls[g]))
for j, ecl in enumerate(ecls[g]):
	vals = ms[ecls==ecl]
	errs = ss[ecls==ecl]
	av, averr = mg.weightedAv(vals, errs)
	avs[j] = av
	averrs[j] = averr



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
fig, (ax1,ax2,ax3,ax4) = mg.formatSubplots((4,1),ylabel=par, sharex=True, sharey=True, figsize=(8,11))
ymin, ymax = 0.01, 0.07
for ax in [ax1,ax2,ax3,ax4]:
	plt.sca(ax)
	plt.xticks([2,5.5,10,14,18.5,23.5],labels, rotation=17)
	#plt.xticks(np.arange(1,25), labels, rotation=17)
	plt.gca().set_xlim(0,26)
	plt.gca().set_ylim(ymin, ymax)
	
	## Lines to separate days
	plt.vlines([3.5,7.5,12.5,15.5,21.5], -ymin, ymax*3, linestyle='solid')

## Plot data
ax1.errorbar(ecls[i], ms[i], ss[i], fmt='m.')
ax1.errorbar(ecls[r], ms[r], ss[r], fmt='r.')
ax2.errorbar(ecls[g], ms[g], ss[g], fmt='g.')
ax3.errorbar(ecls[b], ms[b], ss[b], fmt='b.')
ax4.errorbar(ecls[g], avs, averrs, fmt='k.')

## plot average
ax4.hlines(np.average(avs, weights=1/averrs**2), 0, 26, linestyle='dashed')


mg.saveAsTemp()

plt.show()






