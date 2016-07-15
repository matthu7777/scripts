#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from trm import dnl, subs, sla, roche
from scipy.optimize import curve_fit
from sys import argv, exit
import stars, runs
from subprocess import call
from os import listdir
import mgutils as mg

""" My own simplified version of Tom Marsh's visualise program, plotting using pyplot and Tom's roche library. Always shows top-down.
    Takes a .mod file as its argument. Optionally takes a second argument, phase, for purposes of drawing shadow.
"""


fname = argv[1]

if len(argv) > 2:
	phase = float(argv[2])	# For shadow
else:
	phase = None 		# For shadow

fillShapes = True

def readModFile(fname):
	""" Read all parameters from a parameter file used for lroche.
	"""
	labels = np.loadtxt(fname, usecols=(0,), dtype=str)
	vals = np.loadtxt(fname, usecols=(2,), dtype=str)
	
	# weird parameters that cause numpy errors -- don't need them anyway
	labels = labels[(labels!="limb1")&(labels!="limb2")]
	vals = vals[vals!="Poly"].astype(float)
		
	return labels, vals

def spotProfile(x,blen,bexn,bexm):
	""" Returns the bright spot profile as a function of distance along spot (x).
		x should be in units of angular separation with x=0 at the peak of the profile (x=0 at radius_spot)
	"""
	xmax = (bexn / bexm) ** (1 / bexm) * blen	# the distance from flux=0 to flux=max
	profile = (((x+xmax)/blen)**bexn)*np.exp(-((x+xmax)/blen)**bexm)
	return profile / profile.max()

labels, vals = readModFile(fname)

q = vals[labels=='q'][0]
iangle = vals[labels=='iangle'][0]
brad = vals[labels=='radius_spot'][0]
blen = vals[labels=='length_spot'][0]
bang = vals[labels=='angle_spot'][0]
byaw = vals[labels=='yaw_spot'][0]
bexn = vals[labels=='expon_spot'][0]
bexm = vals[labels=='epow_spot'][0]
wdrad = vals[labels=='r1'][0]
drad = vals[labels=='rdisc2'][0]

print "q", q
print "iangle", iangle
print "spot peak rad", brad
print "spot length", blen
print "spot angle", bang
print "spot yaw", byaw
print "disc rad", drad
print "wd rad", wdrad

################## System geometry

## Roche lobes
xl1, yl1 = roche.lobe1(q)
xl2, yl2 = roche.lobe2(q)

## Stream
xs, ys = roche.stream(q, 1./200)

## Disc radius
xd, yd = mg.circle(drad,n=1000)

## Circularisation radius
xc, yc = mg.circle(roche.rcirc(q))

## White dwarf radius
xw, yw = mg.circle(wdrad)

## Bright spot central position
xb, yb, vxb, vyb = roche.bspot(q, brad)

## Shadow
if phase is not None:
	xsh, ysh, sh = roche.shadow(q, iangle, phase, 10000)
	xsh[np.arange(len(xsh)/2)] = xsh[np.arange(len(xsh)/2)][::-1]
	xsh[np.arange(len(xsh)/2,len(xsh))] = xsh[np.arange(len(xsh)/2,len(xsh))][::-1]
	ysh[np.arange(len(ysh)/2)] = ysh[np.arange(len(ysh)/2)][::-1]
	ysh[np.arange(len(ysh)/2,len(ysh))] = ysh[np.arange(len(ysh)/2,len(ysh))][::-1]
	shmask = xsh**2 + ysh**2 < drad**2  	# only plot points within the disc

##  Spot profile
nmax = (bexn / bexm) ** (1 / bexm)	# the number of length scales from flux=0 to flux=max
l = np.linspace(-nmax,20,500) * blen	# offset by lmax so that l=0 corresponds to flux=max
profile = spotProfile(l,blen,bexn,bexm)
halfpoints = l[profile>profile.max()/2][[0,-1]] / blen

## Bright spot ends
theta = np.abs(np.arctan(yb/xb)) + (bang * np.pi / 180 + np.pi / 2) # Azimuthal angle of spot (defined from "East") + angle of spot from radial direction
lx = blen * np.cos(theta)
ly = blen * np.sin(theta)
xb1, yb1 = xb + lx * halfpoints[0], yb + ly  * halfpoints[0]
xb2, yb2 = xb + lx * halfpoints[1], yb + ly * halfpoints[1]

## Beaming direction
scale1, scale2 = 0.1, 0.3		# arbitrary numbers to control start and end points of arrow
xa1 = xb + scale1 * np.sin(theta + byaw*np.pi/180)
ya1 = yb - scale1 * np.cos(theta + byaw*np.pi/180)
xa2 = xb + scale2 * np.sin(theta + byaw*np.pi/180)
ya2 = yb - scale2 * np.cos(theta + byaw*np.pi/180)

##################  Plotting

## Format
plt.figure(1, figsize=(8,6))
mg.formatGraph()
plt.suptitle(mg.escapeForLatex(fname[:-4]))
## Make it square
xhigh = xl2.max() + 0.2
xlow = xl1.min() - 0.2
ylim = (xhigh - xlow) / 2 * 6/8
plt.gca().set_xlim(xlow, xhigh)
plt.gca().set_ylim(-ylim, ylim)

## Filled regions	-- do these first so that the shading is in the background
if fillShapes:
	plt.fill(xd, yd, 'k', alpha=0.5)	# Disc

	## Fill shadow -- This is hacky. Fills in shadow out to a region slightly wider than the disc radius, then fills in  
	## with white the region between the disc and the roche lobe to cover up the overhang.
	## May need to adjust the "fac" parameter if problems occur.
	if phase is not None:
		fac = 1.05		# how far outside the disc radius to paint the shadow.
		plt.fill(xsh[xsh**2 + ysh**2 < (drad*fac)**2],ysh[xsh**2 + ysh**2 < (drad*fac)**2],'b', alpha=0.5)	#shadow
		plt.fill(np.append(xl1,xd),np.append(yl1,yd), 'w', edgecolor='none')		#cover-up

	plt.fill(xw, yw, 'w')		# White dwarf
	plt.fill(xl2, yl2, 'c', alpha=0.5)		# Secondary


## Plot lines and points
if not fillShapes:	# do shadow first, so it doesn't cover other points.
	plt.plot(xsh[shmask], ysh[shmask], 'b-')	# Outline of shadow, if shadow not filled in
plt.plot(xl1,yl1,'k')	# Roche lobes
plt.plot(xl2,yl2,'k')
plt.plot(xs, ys, 'k--')		# Stream
plt.plot(xd, yd, 'k-')		# Disc
plt.plot(xc, yc, 'k--')		# Rcirc
plt.plot(xw, yw, 'k-')		# White dwarf
plt.plot(xb, yb, 'rD')		# Bright spot
plt.plot([xb1,xb2],[yb1,yb2],'r-')	# Bright spot end
plt.arrow(xa1,ya1,xa2-xa1,ya2-ya1, color='r')	# Bright spot beaming direction





## Plot spot profile (figure 2)
plt.figure(2)
ax2 = plt.gca()
ax2.plot(l, profile, 'k')						# profile
ax2.plot(0, profile.max()/2, 'rD')					# centre for comparison with fig 1
ax2.plot(halfpoints*blen, [profile.max()/2]*2, 'r-')	# end points
ax2.set_xlim(l.min(),l.max())

## Show

mg.saveAsTemp(fignum=1)
plt.show()

