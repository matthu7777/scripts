#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from sys import argv, exit
from subprocess import call
from os import listdir
import runs
import mgutils as mg



""" Plots colour-colour plot of white dwarf. Reads from a dat file created by multmcmcresults.py when the flag in that program refreshData=True
"""

datname = "meandepths.dat"


def weightedAv(x, err):
	av = np.average(x, weights=1/err**2)
	averr = np.sqrt(1 / (np.sum(1/err**2)))		# Equation from Tom H's book p50
	return av, averr
	
def magsFromMJy(mj, mjerr):
	""" Converts milliJanskys to AB magnitudes
	"""
	mags = 23.9 - 2.5*np.log10(1000*mj)
	magerr = np.abs(23.9 - 2.5*np.log10(1000*(mj+mjerr)) - mags)
	return mags, magerr

def plotColour(a,aerr,b,berr,c,cerr,fmt='',ax=None):
	""" Function to plot points on colour-colour plot, saving on typing. Plots c-b against b-a
	"""
	if not ax:
		ax = plt.gca()
	if aerr is not None and berr is not None and cerr is not None:
		return ax.errorbar(b-a, c-b, np.sqrt(cerr**2 + berr**2), np.sqrt(berr**2 + aerr**2), fmt=fmt)
	else:
		return ax.plot(b-a, c-b, fmt)
	
def closest(x,y,xs,ys):
	distx = xs - x
	disty = ys - y
	distSq = distx**2 + disty**2
	return ( distSq==distSq.min() )


############# Load measured data

data = np.loadtxt(datname)
ecls = data[:,0]
ms = data[:,1]
ss = data[:,2]

## Masks
ok = np.ones(len(ecls),dtype=bool)
ok[ecls==10] = False
i = data[:,3].astype(bool)
r = data[:,4].astype(bool)
g = data[:,5].astype(bool)
b = data[:,6].astype(bool)
jan = ecls <= 15
late = ecls > 15
i1 = (i)&(ok)
r1 = (r)&(jan)&(ok)
r2 = (r)&(late)&(ok)
g1 = (g)&(jan)&(ok)
g2 = (g)&(late)&(ok)
b1 = (b)&(jan)&(ok)
b2 = (b)&(late)&(ok)






## Find averages
avi, avierr = magsFromMJy(*weightedAv(ms[i1], ss[i1]))
avr1, avr1err = magsFromMJy(*weightedAv(ms[r1], ss[r1]))
avr2, avr2err = magsFromMJy(*weightedAv(ms[r2], ss[r2]))
avg1, avg1err = magsFromMJy(*weightedAv(ms[g1], ss[g1]))
avg2, avg2err = magsFromMJy(*weightedAv(ms[g2], ss[g2]))
avb1, avb1err = magsFromMJy(*weightedAv(ms[b1], ss[b1]))
avb2, avb2err = magsFromMJy(*weightedAv(ms[b2], ss[b2]))
avr, avrerr = magsFromMJy(*weightedAv(ms[r], ss[r]))	#currently 'ok' mask not included
avg, avgerr = magsFromMJy(*weightedAv(ms[g], ss[g]))
avb, avberr = magsFromMJy(*weightedAv(ms[b], ss[b]))

avi -= 0.28 	## Roughly compensate for calibration problem -- TEMP


############# Read in Bergeron models
teff, logg, mass, mbol, uberg, gberg, rberg, iberg, age = np.loadtxt("bergeron_models.dat", usecols=(0,1,2,3,13,14,15,16,27), unpack=True, skiprows=2)




############# Find closest model
cl = closest(avg2-avr2, avb2-avg2, gberg-rberg, uberg-gberg)

print "Closest point has:"
print "Temperature", teff[cl]
print "Log g", logg[cl]
print "Mass", mass[cl]
print "MBol", mbol[cl]
print "Age", age[cl]

#print 

print "\nExpected u/g,g/r flux ratios", mg.mJyFromMags(uberg[cl],0)[0] / mg.mJyFromMags(gberg[cl],0)[0], mg.mJyFromMags(gberg[cl],0)[0] / mg.mJyFromMags(rberg[cl],0)[0]
print "Measured ratios", avb2/avg2, avg2/avr2


############# Plot
## Format
fig, (ax1, ax2) = mg.formatSubplots((2,1), figsize=(8,10))
ax1.set_xlabel("g'-r'")
ax1.set_ylabel("u'-g'")
ax2.set_xlabel("g'-i'")
ax2.set_ylabel("u'-g'")
ax1.invert_yaxis()
ax2.invert_yaxis()
ax1.set_ylim(0.2,-0.2)
ax2.set_ylim(0.2,-0.2)
ax1.set_xlim(-0.4,-0.2)
ax2.set_xlim(-0.6,-0.4)


## Data
# Individual eclipses
plotColour(*magsFromMJy(ms[r2],ss[r2])+magsFromMJy(ms[g2],ss[g2])+magsFromMJy(ms[b2],ss[b2]), ax=ax1, fmt='m,')

# combined
plotColour(avr1,avr1err, avg1,avg1err, avb1,avb1err, ax=ax1, fmt='b,')
plotColour(avr2,avr2err, avg2,avg2err, avb2,avb2err, ax=ax1, fmt='g,')
#plotColour(avr,avrerr, avg2,avg2err, avb,avberr, ax=ax1, fmt='r,')### temp
plotColour(avr,avrerr, avg,avgerr, avb,avberr, ax=ax1, fmt='r,')### temp?
plotColour(avi,avierr, avg1,avg1err, avb1,avb1err, ax=ax2)




## Model
plotColour(rberg,None, gberg,None, uberg,None, 'k.', ax=ax1)
plotColour(rberg[cl],None, gberg[cl],None, uberg[cl],None, 'r.', ax=ax1)
plotColour(iberg,None, gberg,None, uberg,None, 'k.', ax=ax2)
plotColour(iberg[cl],None, gberg[cl],None, uberg[cl],None, 'r.', ax=ax2)

## Plot favourite model (probably temp)
plt.figure("Fave")
plt.plot([3543,4770,6231,7625],[uberg[cl],gberg[cl],rberg[cl],iberg[cl]])
plt.gca().invert_yaxis()


mg.saveAsTemp(fig=fig)
plt.show()




