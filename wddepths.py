#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from trm import dnl, subs, sla
from scipy.optimize import curve_fit
from sys import argv, exit
from subprocess import call
from os import listdir
import stars, runs
import mgutils as mg

""" Reads white dwarf depths straight from mod files -- probably outdated, use colours.py 
"""



skip = []#["may-g-6", "may-b-5", "may-r-5", "june-g-2", "may-b-6", "may-g-4", "june-g-4", "may-g-3", "may-r-1", "may-b-1", "may-g-2", "jan14-b-3", "may-r-2", "may-g-1", "may-b-2", "may-r-3", "may-b-3"]


if len(argv) > 1 and argv[1] == "update":
	depths = []
	fnames = []
	for f in sorted(listdir(".")):
		if f[-4:] == ".log" and f[:-6] not in skip:
			logname = f
			modname = f[:-4]+".mod"
			initname = "gaia1.mod"
			datname = f[:-6]+".dat"
			method = "p"
			
			call(["extbmod.py",initname,logname, modname, method])
			
			call(["splitter.csh", modname, datname])
			
			a = np.loadtxt("a.out", usecols=(2,))
			
			depth = a.max()
			
			depths = np.append(depths, [depth])
			fnames = np.append(fnames, [f])
	
			print depths
			
			if "-r" in f:
				fmt = "r."
			if "-g" in f:
				fmt = "g."
			if "-b" in f:
				fmt = "b."
			
			num = runs.eclnums[f[:-6]]
			plt.figure(1)
			plt.plot(num, depth, fmt)
	
	np.savetxt( "wddepths.txt", np.vstack((fnames, depths)).T, fmt="%s" )
else:
	fnames = np.loadtxt("wddepths.txt", usecols=(0,), dtype=str)
	depths = np.loadtxt("wddepths.txt", usecols=(1,))
	for i in range(len(fnames)):
		f = fnames[i]
		
		if f[:-6] not in skip:
			if "-r" in f:
				fmt = "r."
			if "-g" in f:
				fmt = "g."
			if "-b" in f:
				fmt = "b."
			
			num = runs.eclnums[f[:-6]]
			plt.figure(1)
			plt.plot(num, depths[i], fmt)

plt.vlines([3.5,7.5,12.5,15.5,21.5],0.08,0.16)

mags, merr = mg.magsFromMJy(depths, np.zeros(len(depths)))

nums = []
colours = []
for f in fnames:
		nums += [runs.eclnums[f[:-6]]]
		if "-r" in f:
			colours += ["r"]
		if "-g" in f:
			colours += ["g"]
		if "-b" in f:
			colours += ["b"]

nums = np.array(nums)
colours = np.array(colours)
enums = range(1,26)
b = np.array([])
r = np.array([])
g = np.array([])
for n in enums:
	if not mags[(nums==n) & (colours=="r")]:
		r = np.append(r, [0])
	else:
		r = np.append(r, mags[(nums==n) & (colours=="r")])
	if not mags[(nums==n) & (colours=="g")]:
		g = np.append(g, [0])
	else:
		g = np.append(g, mags[(nums==n) & (colours=="g")])
	if not mags[(nums==n) & (colours=="b")]:
		b = np.append(b, [0])
	else:
		b = np.append(b, mags[(nums==n) & (colours=="b")])

np.savetxt( "wdmags.txt", np.vstack((enums, r, g, b)).T, fmt="%s" )

plt.figure(2)
plt.plot(b-g, g-r, 'k.')
plt.gca().invert_yaxis()
for n in enums:
	plt.annotate(n, ((b-g)[n-1], (g-r)[n-1]))

plt.show()
		
		
		
