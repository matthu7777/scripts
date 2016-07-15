#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs
from trm import dnl, subs, sla
from scipy.optimize import curve_fit
from sys import argv, exit
import stars, runs
import mgutils as mg


grid = gs.GridSpec(25,1,wspace=0,hspace=0)

colours = ['r','g','b']

address = "/storage/astro2/phulbz/gaia14aae/combined\ fitting/"

for run in runs.runoblist:
	
	rx, rxerr, ry, ryerr, gx, gxerr, gy, gyerr, bx, bxerr, by, byerr = mg.ulg2mjysWithComp(address+run.name+".log", stars.gaia, stars.comp2, run, comp=2)
	
	for num in range(1,run.numec+1):
			name = run.name + "-" + colour + "-" + str(num)
			
			#data = np.loadtxt(name+".dat")
			
			bins=np.linspace(-0.5,0.5,100+1)	#doesn't matter
			btdb = mg.correctTimes(rx, mg.coords, mg.obsname)
			phase, phaseCont, pherr, mids, miderrs, means, errs, counts = mg.phaseFold(btdb, rxerr, ry, mg.t0, mg.t0err, mg.period, mg.perr, bins, yerr)
			
			phmin = np.around((phaseCont+0.5).min(),0) - 0.5
			phmax = np.around((phaseCont+0.5).max(),0) - 0.5
			
			
			
			
			
			
			
			
			
			i = runs.eclnums[name] - 1
			ax = plt.subplot(grid[i,0])
			mg.lightcurve(data[:,0], data[:,2], data[:,3], fmt=colour+'.', sub=ax)
			
			



	
	





plt.show()






