#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from sys import argv
import getopt
from subprocess import call
from os import listdir
from scipy.optimize import curve_fit
import runs
import mgutils as mg

help = """ Script to compare the results of multiple mcmc runs on different data sets, tailored to eclipse depths of gaia14aae but can plot other parameters too

    If refreshData, recalculates mean parameter values and stds from log files; otherwise, reads from text file "datname".
    If testData, calculates goodness of fit for each mcmc result against a normal distribution.
    If plotHists, plots histograms one-by-one so you can check normality by eye.
    If saveData, saves data to text file "datname".
    If doAnalysis, performs some analysis specified lower down in this script -- eg averages. 

    Arguements:
    Parameter of interest should be specified as argv[1], or will default to wdwarf depth.
    Text file to save to/read from should be named in "datname"
    Set of log files should must be listed in a separate file, which is named in "files" variable
    Chisq above which to ignore should be specifed as "cburn"
    Boolean flags should be specified either as arguments or below.
"""
usage = "multmcmcresults.py -l <listfile> -t|d <textfile> -c <cburn> --refreshdata= --testdata= --plothists= --savedata= --doanalysis="


files = "outfiles4.lis"
datname = "meandepths.dat"

refreshData = True		#Re-import data from log files? (else reads means from a text file)
testData = True		#Test normality of data distribution?
plotHists = True	#Plot histograms of individual mcmcs?
saveData = True
doAnalysis = False
cburn = None


##  Handle arguments

def evalBool(string):
	if string == "True":
		return True
	elif string == "False":
		return False
	else:
		print "Input not recognised as True or False, defaulting to False"
		return False

if len(argv) > 1:
	if argv[1] == "-h":
		print help
		print usage
	par = argv[1]
else:
	par = "wdwarf"

if len(argv) > 2:
	try:
		opts, args = getopt.getopt(argv[2:],"hl:t:d:c:",["list=","text=","refreshdata=","testdata=","plothists=","savedata=","doanalysis=","cburn="])
	except getopt.GetoptError:
		print usage
		exit(2)

	for opt, arg in opts:
		if opt == '-h':
			print help
			print usage
			exit()
		elif opt in ['-l','--list']:
			files = arg
		elif opt in ['-t','--textfile','-d']:
			datname = arg
		elif opt.lower() in ['--refreshdata']:
			refreshData = evalBool(arg)
		elif opt.lower() in ['--testdata']:
			testData = evalBool(arg)
		elif opt.lower() in ['--plothists']:
			plotHists = evalBool(arg)
		elif opt.lower() in ['--savedata']:
			saveData = evalBool(arg)
		elif opt.lower() in ['--doanalysis']:
			doAnalysis = evalBool(arg)
		elif opt.lower() in ['-c','--cburn']:
			cburn = float(arg)


if plotHists:
	testData = True		# Otherwise it won't be possible.
	refreshData = True
if testData:
	refreshData = True



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
	print "Read %d points from %s"%(len(oi),fname)
	return oi

def mcmcHist(data, plotHists=True, label='', chisq=None, cburn=None):
	"""	Returns the mean and std of a dataset, as well correlation with a Gaussian distribution. Can optionally plot a histogram
	"""
	## Allow for some high chisq points to be ignored
	if cburn is not None and chisq is not None:
		ok = chisq < cburn
		if not ok.any():
			print "No points above cburn"
		else:
			print "%s points above cburn"%(np.sum(ok))
	else:
		ok = np.ones(len(data),dtype=bool)

	## Get histogram of results
	if plotHists:
		# This one plots histogram and returns values
		n, bins, patches = plt.hist(data[ok], bins=20)
	else:
		# This one just returns values, without plotting
		n, bins = np.histogram(data[ok], bins=20)
	
	## Fit with Gaussian
	try:
		mids = (bins[1:] + bins[:-1]) / 2
		a0 = data[ok].max()
		m0 = data[ok].mean()
		s0 = data[ok].std()
		popt, pcov = curve_fit(mg.gauss, mids, n, p0=[a0,m0,s0])
		
		a, m, s = popt
		
		rms = np.sqrt( ((n - mg.gauss(mids, *popt))**2).mean() )
		
		print label, "RMS from Gaussian distribution", rms
		
		if plotHists:
			# Don't bother plotting gaussian, saves a bit of time
			plt.plot(mids, mg.gauss(mids, *popt), 'r')
	except RuntimeError:
		print label, "did not converge"
	
	if plotHists:
		if label:
			plt.suptitle(label)
		mg.formatGraph()
		plt.show()
	return data[ok].mean(), data[ok].std(), rms
	

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
		if name in runs.eclnums:
			eclnum = runs.eclnums[name]		# if we're looking at gaia14aae single eclipses, it's useful to know which one is which
		else:
			eclnum = j 					# otherwise just number them in order
		
		if refreshData:
			## read in data and find means
			print "Reading data from", fname
			data = readMcmcLog(fname, par)
			c2 = readMcmcLog(fname, "chisq")

			## If desired, report goodness of fit of data to a Gaussian
			if testData is True:
				m, s, rms = mcmcHist(data, plotHists=plotHists, label=name, chisq=c2, cburn=cburn)
			else:
				if cburn:
					ok = c2 < cburn
				else:
					ok = np.ones(len(data),dtype=None)
				m, s = data[ok].mean(), data[ok].std()
			
			print m 
			
			ecls[j] = eclnum
			ms[j] = m
			ss[j] = s
		
		
		
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
	print "Loading data from", datname
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
if doAnalysis:
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



		if c is not i:
			## Individual nights? Jan / May?
			jan = (ecls <= 15)
			late = (ecls > 15)
			avgjan = np.average( ms[(c)&(jan)], weights=((1/ss[(c)&(jan)])**2))
			chisqjan = chisq(avgjan, ms[(c)&(jan)], ss[(c)&(jan)])
			avglate = np.average( ms[(c)&(late)], weights=((1/ss[(c)&(late)])**2))
			chisqlate = chisq(avglate, ms[(c)&(late)], ss[(c)&(late)])

			print "January average %s"%(col), avgjan, "Chisquared", chisqjan
			print "Late average %s"%(col), avglate, "Chisquared", chisqlate
			print "Total split chisquared", chisqjan+chisqlate, "\n"


			iterjan, pcov, maskj = mg.iterativeFit(mg.flat, np.arange(len(ms[(c)&(jan)])), ms[(c)&(jan)], ss[(c)&(jan)], 10, 3, p0=[avg])
			iterlate, pcov, maskl = mg.iterativeFit(mg.flat, np.arange(len(ms[(c)&(late)])), ms[(c)&(late)], ss[(c)&(late)], 10, 3, p0=[avg])
			print "Iterative average jan", iterjan[0], "masked", len(maskj[maskj==False])
			print "Iterative average late", iterlate[0], "masked", len(maskl[maskl==False])

			#plt.hlines(avgjan, 0, 15, linestyle='dashed', color=col)
			#plt.hlines(avglate, 16, 25, linestyle='dashed', color=col)







################ Plotting


mg.formatGraph(ylabel=par)

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
#plt.xticks(np.arange(1,25), labels, rotation=17)
#ymin, ymax = 0.04, 0.18
if par == "wdwarf":
	ymin, ymax = 0.01, 0.07
else:
	ymin, ymax  = mg.adjLimits(ms)
plt.xticks([2,5.5,10,14,18.5,23.5],labels, rotation=17)
plt.gca().set_xlim(0,26)
plt.gca().set_ylim(ymin, ymax)

## Plot data
if i.any() or r.any() or g.any() or b.any():
	# If colour masks have been defined, plot data according to colour
	plt.errorbar(ecls[i], ms[i], ss[i], fmt='m.')
	plt.errorbar(ecls[r], ms[r], ss[r], fmt='r.')
	plt.errorbar(ecls[g], ms[g], ss[g], fmt='g.')
	plt.errorbar(ecls[b], ms[b], ss[b], fmt='b.')
else:
	# Otherwise plot it all the same
	plt.errorbar(ecls, ms, ss, fmt='k.')

## Lines to separate days
plt.vlines([3.5,7.5,12.5,15.5,21.5], ymin-np.abs(ymin)*3, np.abs(ymax)*3, linestyle='solid')

mg.saveAsTemp()

plt.show()






