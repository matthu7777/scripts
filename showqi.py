#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from trm import dnl, subs, sla, roche
from sys import argv
import stars, runs
import mgutils as mg

""" Script to plot q/i relations for roche-filling stars. Reads output of format put out by readqi.py
"""
if len(argv) > 1:
	fname = argv[1]
else:
	fname = "qi.dat"

ilowerlim = 80

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



## Find weighted average min q and do some statistics

avmq, avmqerr = mg.weightedAv(minq[nz],minqerr[nz])

print "\nStraight mean min q is", minq[nz].mean()
print "Weighted mean min q is", avmq, "+/-", avmqerr
print "Standard deviation of all min q is", minq[nz].std()

avphi, avphierr = mg.weightedAv(phi[nz], phierr[nz])
print "\nWeighted average phase width is", avphi, "+/-", avphierr
print "Min q from this is", roche.findq(90, avphi), "+/-", np.fabs(roche.findq(90, avphi+avphierr) - avphi)


## Create model q/i relationship based on mean phi
modelq, modelplus, modelminus = [], [], []
modeli = np.linspace(ilowerlim, 90, 1000)
for j in modeli:
	modelq += [roche.findq(j, avphi)]
	modelplus += [roche.findq(j, avphi+avphierr)]
	modelminus += [roche.findq(j, avphi-avphierr)]


## Find statistics of points from model
s = []
chisq, chisqr, chisqg, chisqb = 0, 0, 0, 0
for j in xrange(len(q[nz])):
	s += [(i[nz][j] - roche.findi(q[nz][j], avphi)) ** 2]		#find difference squared for each point
	chisq += ((i[nz][j] - roche.findi(q[nz][j], avphi)) / ierr[nz][j])**2 + ((q[nz][j] - roche.findq(i[nz][j], avphi)) / qerr[nz][j])**2
	if "-r" in names[j]:
		chisqr += ((i[nz][j] - roche.findi(q[nz][j], avphi)) / ierr[nz][j])**2 + ((q[nz][j] - roche.findq(i[nz][j], avphi)) / qerr[nz][j])**2
	if "-g" in names[j]:
		chisqg += ((i[nz][j] - roche.findi(q[nz][j], avphi)) / ierr[nz][j])**2 + ((q[nz][j] - roche.findq(i[nz][j], avphi)) / qerr[nz][j])**2
	if "-b" in names[j]:
		chisqb += ((i[nz][j] - roche.findi(q[nz][j], avphi)) / ierr[nz][j])**2 + ((q[nz][j] - roche.findq(i[nz][j], avphi)) / qerr[nz][j])**2
	
rms = np.sqrt(np.array([s]).mean())						# take mean and square root
print "\nRMS is", rms
print "Chisq is", chisq
print "Reduced chisq is", chisq / len(q[nz])
print "\nChisq of red points is", chisqr, "reduced is", chisqr / len(q[nz][red])
print "Chisq of green points is", chisqg, "reduced is", chisqg / len(q[nz][gre])
print "Chisq of blue points is", chisqb, "reduced is", chisqb / len(q[nz][blu])



## Create and format figures
fig, (ax1, ax2, ax3) = mg.formatSubplots((3,1),sharex=True,sharey=True, xlabel="Mass ratio $q$", ylabel="Inclination $i$, deg", figsize=(8,11))
for ax in (ax1,ax2,ax3):
	ax.invert_yaxis()


	
## Plot model q/i relationship
for ax in (ax1,ax2,ax3):
	ax.plot(modelq, modeli, 'k-')
	ax.plot(modelplus, modeli, 'k--')
	ax.plot(modelminus, modeli, 'k--')


## Plot measured phases and inclinations
ax3.errorbar(q[nz][blu],i[nz][blu],ierr[nz][blu],xerr=qerr[nz][blu],fmt='b.')
ax2.errorbar(q[nz][gre],i[nz][gre],ierr[nz][gre],xerr=qerr[nz][gre],fmt='g.')
ax1.errorbar(q[nz][red],i[nz][red],ierr[nz][red],xerr=qerr[nz][red],fmt='r.')

mg.saveAsTemp()


print "\n"
plt.show()



