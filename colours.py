#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs
from trm import dnl, subs, sla
from scipy.optimize import curve_fit
from sys import argv, exit
from subprocess import call
from os import listdir
import stars, runs
import mgutils as mg

""" Reads from a csv containing eclipse depths and plots them as well as colour-colour plot"""

# read data
nights = np.loadtxt("ecldepths.csv", delimiter=',', usecols=(0,), dtype=str)
eclnums = np.loadtxt("ecldepths.csv", delimiter=',', usecols=(1,), dtype=str)
data = np.loadtxt("ecldepths.csv", delimiter=',', usecols=(2,3,4,5,6,7))


# set plot colours
cols = {"jan14":'r.',"jan15":'y.',"jan16":'g.',"jan17":'c.',"may":'b.',"june":'m.'}


# plot settings
plt.figure(1,figsize=(8,11), dpi=160)		# A4 ratio
grid = gs.GridSpec(3,1,wspace=0,hspace=0.15)
subgrid = gs.GridSpecFromSubplotSpec(2,1,subplot_spec=grid[1:3,0],wspace=0,hspace=0)
ax1 = plt.subplot(grid[0,0])
ax2 = plt.subplot(subgrid[0,0])
ax3 = plt.subplot(subgrid[1,0])
# font
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=16)
# axis labels
ax1.set_xticklabels([])
ax1.set_xlabel("Eclipse number")
ax1.set_ylabel("WD flux (mJy)")
ax1.set_xlim(0,26)
ax2.set_ylabel("u' - g' Colour")
ax2.set_xticklabels([])
ax2.set_xlim(-0.39,-0.20)
ax2.set_ylim(0.16,-0.08)
ax3.set_ylabel("u' - g' Colour")
ax3.set_xlabel("g' - r' Colour")
ax3.set_xlim(-0.39,-0.20)
ax3.set_ylim(0.16,-0.08)



r = data[:,0]
rerr = data[:,1]
g = data[:,2]
gerr = data[:,3]
b = data[:,4]
berr = data[:,5]

temp = np.core.defchararray.add(nights, "-")
labels = np.core.defchararray.add(temp, eclnums)

# plot mjys
for i in range(len(r)):
		
	ax1.errorbar(runs.eclnums[labels[i]],r[i],rerr[i],fmt='r.')
	ax1.errorbar(runs.eclnums[labels[i]],g[i],gerr[i],fmt='g.')
	ax1.errorbar(runs.eclnums[labels[i]],b[i],berr[i],fmt='b.')


ax1.vlines([3.5,7.5,12.5,15.5,21.5],0.08,0.16, linestyles='dashed')	#separate nights by lines


# plot colours
for i in range(len(r)):
	
	rm, rmerr = mg.magsFromMJy(r[i], rerr[i])
	gm, gmerr = mg.magsFromMJy(g[i], gerr[i])
	bm, bmerr = mg.magsFromMJy(b[i], berr[i])
	
	# one figure with error bars
	ax2.errorbar(gm-rm, bm-gm, np.sqrt(gmerr**2 + bmerr**2), np.sqrt(rmerr**2 + gmerr**2), fmt=cols[nights[i]])
	
	# one figure without error bars, annotated
	ax3.plot(gm-rm, bm-gm, cols[nights[i]])
	ax3.annotate(labels[i],(gm-rm, bm-gm), color=cols[nights[i]][0], size=8)
	
#ax2.invert_yaxis()
#ax3.invert_yaxis()

plt.savefig("/storage/astro2/phulbz/gaia14aae/figures/colours.pdf",dpi='figure',bbox_inches='tight')
plt.savefig("/storage/astro2/phulbz/gaia14aae/figures/colours.png",dpi='figure',bbox_inches='tight')

plt.show()

