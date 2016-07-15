#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from trm import dnl, subs, sla, roche
import stars, runs
import mgutils as mg
from os import listdir
import pickle

""" Script for plotting things in ingress/egress space
"""

avphi = 0.0373289888593

ianglemax=90
qmin = roche.findq(ianglemax, avphi)

q = 0.024
iangle = roche.findi(q, avphi)
q2 = 0.03
iangle2 = roche.findi(q2, avphi)

plotFoldedEclipses = True
plotUnfoldedEclipses = False

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
	
	## Play around with wrapping
	ings[ings>0.5] = ings[ings>0.5] - 1
	egs[egs>1] = egs[egs>1] -1
	
	
	
	return ings, egs, xstream, ystream
	
def toingeg(x,y):
	""" Converts x,y coords to ingress, egress phases. Wrapper around trm.roche.ineg that does the iteration for you
	"""
	ings,egs = [],[]

	for j in xrange(len(x)):
		ing, eg = roche.ineg(q, iangle, x[j], y[j])
		ings += [ing]
		egs += [eg]
		
	return np.array(ings), np.array(egs)


def circle(q):
	""" Returns coordinates of a circle at the circularisation radius for a given q
	"""
	rcirc = roche.rcirc(q)
	angle = np.linspace(0,2*np.pi,200)
	
	x = rcirc * np.sin(angle)
	y = rcirc * np.cos(angle)
	
	return x, y

############# Maths and models	

## Find coordinates for streams in the two cases we're interested in

ings, egs, xstream, ystream = iestream(q, iangle, 1./1000, 1000)
ings2, egs2, xstream2, ystream2 = iestream(q2, iangle2, 1./1000, 1000)
ingsl, egsl, xstreaml, ystreaml = iestream(qmin, ianglemax, 1./1000, 1000)


## Find Roche lobe of primary
xrl1, yrl1 = roche.lobe1(q)
xrl2, yrl2 = roche.lobe1(q2)
xrll, yrll = roche.lobe1(qmin)

## Find circularisation radius of disc
xc1, yc1 = circle(q)
xc2, yc2 = circle(q2)
xcl, ycl = circle(qmin)
inc1, egc1 = toingeg(xc1[(xc1>0) & (yc1>0.2)], yc1[(xc1>0) & (yc1>0.2)])
inc2, egc2 = toingeg(xc2[(xc2>0) & (yc2>0.2)], yc2[(xc2>0) & (yc2>0.2)])
incl, egcl = toingeg(xcl[(xcl>0) & (ycl>0.2)], ycl[(xcl>0) & (ycl>0.2)])


############## Measured data 

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



################ Plotting

## Make graphs pretty
plt.style.use('seaborn-colorblind')
mg.formatGraph()
ax1 = plt.gca()
ax1.invert_xaxis()
ax1.set_xlabel("Ingress phase")
ax1.set_ylabel("Egress phase")
ax1.set_xlim(0.05,0)
ax1.set_ylim(0.02,0.12)

## Plot phase space graph
ax1.plot(ings, egs, 'k--')	# plot streams
ax1.plot(ings2, egs2, 'k--')
ax1.plot(ingsl, egsl, 'k-')

#ax1.plot(inc1, egc1, 'g--')	# circ radii
#ax1.plot(inc2, egc2, 'b--')
#ax1.plot(incl, egcl, 'k--')

## Plot eclipses
if plotUnfoldedEclipses:
	ax1.errorbar(eings[er], eegs[er], eegerrs[er], eingerrs[er], 'r.')	#unfolded
	ax1.errorbar(eings[eg], eegs[eg], eegerrs[eg], eingerrs[eg], 'g.')
	ax1.errorbar(eings[eb], eegs[eb], eegerrs[eb], eingerrs[eb], 'b.')
if plotFoldedEclipses:
	ax1.errorbar(eings[efr], eegs[efr], eegerrs[efr], eingerrs[efr], 'r.')	#folded
	ax1.errorbar(eings[efg], eegs[efg], eegerrs[efg], eingerrs[efg], 'g.')
	ax1.errorbar(eings[efb], eegs[efb], eegerrs[efb], eingerrs[efb], 'b.')



## Plot real space graph
plt.figure(2)
plt.plot(xstream, ystream, 'g-')
plt.plot(xstream2, ystream2, 'b-')
plt.plot(xstreaml, ystreaml, 'k-')
plt.plot(xrl1, yrl1, 'g--')
plt.plot(xrl2, yrl2, 'b--')
plt.plot(xrll, yrll, 'k--')
plt.plot(xc1, yc1, 'g--')
plt.plot(xc2, yc2, 'b--')
plt.plot(xcl, ycl, 'k--')


plt.figure(1)
mg.saveAsTemp()

plt.show()
