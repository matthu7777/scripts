#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from sys import argv, exit
from subprocess import call
from os import listdir
import runs
from trm import roche

if len(argv) > 1:
	origname=argv[1]
else:
	origname = "orig.mod"

""" Creates a set of lroche model files (.mod) with a range of mass ratios, and corresponding inclinations.
"""


def readMod(fname):
	""" Read in an lroche model file (.mod). Returns columns, as well as the extra stuff at the bottom.
	"""
	labels, vals, a, b, c, d = np.genfromtxt(fname, usecols=(0,2,3,4,5,6),dtype=str,unpack=True,skip_footer=34)
	footer = np.loadtxt(fname, skiprows=58,dtype=str)
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



## Read in original model file
labels, valsorig, a, b, c, d, footer = readMod(origname)
qorig = valsorig[labels=='q']
iorig = valsorig[labels=='iangle']
phi = roche.findphi(qorig, iorig)
minq = roche.findq(90,phi)

## Loop over q values wanted
qlist = np.append([minq],np.linspace(0.02,0.1,9))
for q in qlist:
	## find i and write model
	iangle = roche.findi(q,phi)

	vals = np.copy(valsorig)
	vals[labels=='q'] = q
	vals[labels=='iangle'] = iangle

	if q == minq:
		outname = "qmin.mod"
	else:
		outname = "q%03d.mod"%(q*100)

	writeMod(outname, labels, vals, a, b, c, d, footer)

	print "Written model with q", q, "to", outname

