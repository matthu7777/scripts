#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs
import mgutils as mg
import stars, runs

""" Separates runs into different eclipses, and can also plot stacked figures of different eclipses (for residual stacked figure, use subtracteclipses.py in ../residuals folder
"""


runlist = [runs.jan14, runs.jan15, runs.jan16, runs.jan17, runs.may, runs.june]
	
colours = ['r', 'g', 'b']

plot = False
plotall = False


if plot:
	if not plotall:			#plots on multiple figures that can be printed on different pages
		for i in range(9):
			plt.figure(i,figsize=(8,11), dpi=160)
		grid = gs.GridSpec(3,1,wspace=0,hspace=0)
	else:
		plt.figure(1,figsize=(8,80), dpi=160)
		grid = gs.GridSpec(25,1,wspace=0,hspace=0)
	

for run in runlist:
	loadname = run.locn+".log"
	
	runname = run.name		#format will be runname-colour-number.dat
	
	r, g, b = mg.readulog(loadname)		#This just to create object - data not used
	r.x.data, r.x.errors, r.y.data, r.y.errors, g.x.data, g.x.errors, g.y.data, g.y.errors, b.x.data, b.x.errors, b.y.data, b.y.errors \
									= mg.ulg2mjysWithComp(loadname, stars.gaia, stars.comp3, runs.may, oi=1, comp=3)
	
	for i in range(len([r,g,b])):
		l = [r,g,b][i]
			
		btdb = mg.correctTimes(l.x.data, mg.coords, mg.obsname)
	
		bins = np.linspace(-0.5,0.5,1000+1)	#doesn't matter
		
		phase, phaseCont, pherr, mids, miderrs, means, errs, counts \
							= mg.phaseFold(btdb, l.x.errors, l.y.data, mg.t0, mg.t0err, mg.period, mg.perr, bins, l.y.errors)	# the only thing we need from here is phaseCont and its error
		
		phmin = np.around((phaseCont+0.5).min(),0) - 0.5
		phmax = np.around((phaseCont+0.5).max(),0) - 0.5
		
		print phaseCont.min(), phmin
		print phaseCont.max(), phmax
		
		for centre in range(int(phmin+0.5), int(phmax+0.5)):
			oi = (phaseCont > (centre - 0.5)) & (phaseCont < (centre + 0.5))
			
			j = centre - phmin + 1
			
			phase = phaseCont[oi] - centre
			
			weights = np.ones(len(btdb[oi]))			# some error came up on weights line, dunno what it means 
			weights[(phase < -0.0235)] = 0.1
			weights[(phase > 0.025)] = 0.1
			weights[(phase < -0.09)] = 0.001
			weights[(phase > 0.11)] = 0.001
			#if runname == "may":
				#weights[(phase > -0.025)&(phase < -0.015)] = 100
				#weights[(phase > 0.015)&(phase < 0.025)] = 100
				#weights = weights / 100
			
			savename = runname + "-" + colours[i] + "-%d.dat"%(j)
			
			mg.saveForLcurve(btdb[oi], l.x.errors[oi], l.y.data[oi], l.y.errors[oi], savename, weights)
			
			print "Saved to", savename
			
			if plot:
				k = runs.eclnums[runname+"-"+colours[i]+"-"+"%d"%(j)] - 1
				if not plotall:
					plt.figure(k//3)			#ax = plt.subplot(grid[k,0])
					ax = plt.subplot(grid[k%3,0])
					mg.lightcurve(phase, l.y.data[oi], l.y.errors[oi], sub=ax, fignum=k//3, fmt=colours[i]+'.', phase=True, mj=True, xleft=-0.12, xright=0.2, ylow=0, yhigh=0.22)
					plt.figtext(0.15,0.11+(2-k%3.)/3.75,run.date+" eclipse %d"%(j))
				else:
					ax = plt.subplot(grid[k,0])
					#plt.vlines([0,0.011],-1,1)
					mg.lightcurve(phase, l.y.data[oi], l.y.errors[oi], sub=ax, fmt=colours[i]+'.', phase=True, mj=True, xleft=-0.12, xright=0.2, ylow=0, yhigh=0.22)

if plot:
	if plotall is not True:				
		for i in range(9):
			plt.figure(i)
			plt.savefig("stacked"+str(i)+".pdf",dpi='figure',bbox_inches='tight')
	else:
		plt.savefig("stacked.pdf",dpi='figure',bbox_inches='tight')

plt.show()
		
