#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from trm import dnl, subs, sla, roche
import stars, runs
import mgutils as mg
from os import listdir
import pickle

""" Script for calculating q and i from bright spot contact points
"""

avphi = 0.0373289888593
avphierr = 0.000198442858877

ilowerlim = 80

class Eclipse(object):
	""" A class for handling times of an eclipse.
	"""
	def __init__(self, name):
		self.name = name
		
def iestream(q, iangle, step=1./1000, n=200):
	""" Returns ingress, egress, x and y positions for stream. A wrapper around trm.roche.stream
	"""
	xstream, ystream = roche.stream(q, step, n)
	ings, egs = toingeg(xstream,ystream)
	
	return ings, egs, xstream, ystream
	
def toingeg(x,y):
	""" Converts x,y coords to ingress, egress phases. Wrapper around trm.roche.ineg that does the iteration for you
	"""
	ings,egs = [],[]

	for j in xrange(len(x)):
		ing, eg = roche.ineg(q, iangle, x[j], y[j])
		ings += [ing]
		egs += [eg]
		
	return ings, egs
#############  Read eclipses in

eclipses, eings, eegs, eingerrs, eegerrs = [], [], [], [], []
for fname in listdir("bseclipses"):
	with open("bseclipses/"+fname, 'rb') as f:
		eclipse = pickle.load(f)
		eclipses += [eclipse]

er,eg,eb,efr,efg,efb,cautions = np.zeros(len(eclipses),dtype=bool),np.zeros(len(eclipses),dtype=bool),np.zeros(len(eclipses),dtype=bool),np.zeros(len(eclipses),dtype=bool),np.zeros(len(eclipses),dtype=bool),np.zeros(len(eclipses),dtype=bool),np.zeros(len(eclipses),dtype=bool)
for j in range(len(eclipses)):
		eclipse = eclipses[j]
		eings += [(eclipse.p2 + eclipse.p1) / 2]
		eegs += [(eclipse.p4 + eclipse.p3) / 2]
		eingerrs += [(eclipse.p2 - eclipse.p1) / 2]
		eegerrs += [(eclipse.p4 - eclipse.p3) / 2]
		
		if "-r" in eclipse.name and "f" not in eclipse.name:
			er[j] = True
		if "-g" in eclipse.name and "f" not in eclipse.name:
			eg[j] = True
		if "-b" in eclipse.name and "f" not in eclipse.name:
			eb[j] = True
		if "r" in eclipse.name and "f" in eclipse.name:
			efr[j] = True
		if "r" not in eclipse.name and "b" not in eclipse.name and "f" in eclipse.name:		#why did i use this naming convention?
			efg[j] = True
		if "b" in eclipse.name and "f" in eclipse.name:
			efb[j] = True
		
		if hasattr(eclipse, 'caution') and eclipse.caution:
			cautions[j] = True

eings = np.array(eings)
eegs = np.array(eegs)
eingerrs = np.array(eingerrs)
eegerrs = np.array(eegerrs)

################## Find q and i 
qs, iangles, rbss, ok = [], [], [], []
for j in xrange(len(eings)):
	try:
		q, iangle, rbs = roche.qirbs(avphi,eings[j], eegs[j])
		qs += [q]
		iangles += [iangle]
		rbss += [rbs]
		ok += [True]

	except roche.RocheError:
		print "Skipping", eclipses[j].name
		qs += [0]
		iangles += [0]
		rbss += [0]
		ok += [False]
		continue
	
rbss = np.array(rbss)
qs = np.array(qs)
iangles = np.array(iangles)

########## Create model q/i relationship based on mean phi
modelq, modelplus, modelminus = [], [], []
modeli = np.linspace(ilowerlim, 90, 1000)
for j in modeli:
	modelq += [roche.findq(j, avphi)]
	modelplus += [roche.findq(j, avphi+avphierr)]
	modelminus += [roche.findq(j, avphi-avphierr)]



## Create and format figures
fig, (ax1, ax2, ax3) = mg.formatSubplots((3,1),sharex=True,sharey=True, xlabel="Mass ratio $q$", ylabel="Inclination $i$, deg", figsize=(8,11))
for ax in (ax1,ax2,ax3):
	ax.invert_yaxis()
	ax.set_ylim(90,80)


	
## Plot model q/i relationship
for ax in (ax1,ax2,ax3):
	ax.plot(modelq, modeli, 'k-')
	ax.plot(modelplus, modeli, 'k--')
	ax.plot(modelminus, modeli, 'k--')


## Plot measured phases and inclinations
ax3.plot(qs[eb],iangles[eb],'b.')
ax2.plot(qs[eg],iangles[eg],'g.')
ax1.plot(qs[er],iangles[er],'r.')

mg.saveAsTemp()

plt.show()











