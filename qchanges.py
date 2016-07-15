#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from sys import argv, exit
import shlex, subprocess, re
from os import listdir
import runs
import mgutils as mg


""" Plots the starting q value of levmarq/mcmc runs versus the final value, as well as chi squared for final fits
"""

qlis = "q.lis"		# list of original model files
if len(argv) > 1:
	datfile = argv[1]		# lightcurve data
else:
	datfile = "may-r-2.dat"
mcmcDataFile = "q.dat"			# text file of mcmc results
levmarq = True
mcmc = True

mcmcDict = {"qmin":0,"q002":1,"q003":2,"q004":3,"q005":4,"q006":5,"q007":6,"q008":7,"q009":8,"q010":9}

def readMod(fname):
	""" Read in an lroche model file (.mod). Returns columns, as well as the extra stuff at the bottom.
	"""
	with open(fname) as lis:
		numlines = sum(1 for line in lis)
	labels, vals, a, b, c, d = np.genfromtxt(fname, usecols=(0,2,3,4,5,6),dtype=str,unpack=True,skip_footer=34)
	footer = np.loadtxt(fname, skiprows=(numlines-34),dtype=str)
	vals = vals.astype(float)

	return labels, vals, a, b, c, d, footer

def writeMod(fname,labels, vals, a, b, c, d, footer):
	""" Write a set of columns to a model file. Could take direct output from readMod if so desired.
	"""
	eqs = np.array(['=']*len(labels))
	top = np.column_stack((labels, eqs, vals, a, b, c, d))

	np.savetxt(fname, top, fmt=['%s','%s','%s','%s','%s','%s','%s'])
	with open(fname,'a') as f:
		for line in footer:
			line = line[0]+' '+line[1]+' '+line[2]+'\n'
			f.write(line)
	return True

def findChisq(modfile, datfile):
	""" Read the chisq of an lroche model to a data file, by calling lroche 
	"""
	## Run lroche to get output using subprocess module
	cmd = "/storage/astro1/phsaap/software/bin/lcurve/lroche %s %s device=null scale=yes nfile=0 \\\\"%(modfile, datfile)
	process = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE)		# runs lroche
	output, err = process.communicate()		# output is the stdout of lroche
	exit_code = process.wait()		# exit code, in case needed

	## find chisq in output with regex
	reg = re.compile(r"Weighted chi\*\*2\s\=\s.+")
	line = reg.search(output)
	if line is not None:
		chisq = float(line.group(0).split()[3][:-1])	#convert to string, split line, take number, trim comma, convert to float
	else:
		chisq = -1
	return chisq

################   Read q values 

qorigs, qlevs, qmcmcs, emcmcs, chisqlevs, chisqmcmcs = np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([])
fnames = []
with open(qlis) as lis:
	for fname in lis:
		fname = fname[:-1]	#trim newline
		fnames += [fname]

mcmcData = np.loadtxt(mcmcDataFile)

for fname in fnames:

	## Read starting q
	print "Reading", fname
	labels, vals, a, b, c, d, footer = readMod(fname)
	qorig = vals[labels=='q']

	## Read finishing q
	if levmarq:
		lname = fname[:-4] + "_lev.mod"
		try:
			labels, vals, a, b, c, d, footer = readMod(lname)
			qlev = vals[labels=='q']
			print "Read q"
			chisqlev = findChisq(lname,datfile)
		except (IOError, OSError):
			print "Unable to find file", lname
			qlev = -1
			chisqlev = -1

	if mcmc:
		mnum = mcmcDict[fname[:-4]]
		qmcmc, emcmc = mcmcData[mnum,1:3]

	## Pass out of loop
	qorigs = np.append(qorigs,qorig)
	if levmarq:
		qlevs = np.append(qlevs,qlev)
		chisqlevs = np.append(chisqlevs,chisqlev)
	if mcmc:
		qmcmcs = np.append(qmcmcs,qmcmc)
		emcmcs = np.append(emcmcs,emcmc)

#####################  Plot graphs


## Format
lowlim, highlim = 0.01, 0.11
fig, axs = mg.formatSubplots((2,1))
axs[0].set_xlabel("Starting value q")
axs[0].set_ylabel("Final value q")
axs[0].set_xlim(lowlim, highlim)
axs[0].set_ylim(lowlim, highlim)
axs[1].set_xlabel("Final value q")
axs[1].set_ylabel(r"$\chi ^2$")
axs[1].set_xlim(lowlim, highlim)
axs[1].set_ylim(mg.adjLimits(chisqlevs[chisqlevs!=-1]))

## Plot
if levmarq:
	axs[0].plot(qorigs, qlevs, 'r.')
	axs[1].plot(qlevs, chisqlevs, 'r.')
if mcmc:
	axs[0].errorbar(qorigs, qmcmcs, emcmcs, fmt='b.')

## Diagonal line y=x to guide eye
axs[0].plot([0,1],[0,1],'k--')

mg.saveAsTemp()
#plt.show()



