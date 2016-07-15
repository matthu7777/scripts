#!/usr/bin/env python

import mgutils as mg
import numpy as np
import matplotlib.pyplot as plt
from sys import argv
from os import remove
from os.path import isfile


modname = argv[1]
fname = argv[2]

t, terr, y, yerr = mg.readdat(fname)


readLcurve = mg.readLcurve
	
x, a = readLcurve("a.out")
x, b = readLcurve("b.out")
x, c = readLcurve("c.out")
x, d = readLcurve("d.out")

toff = np.around(x.min(),1)
time = x - toff			#tidying up time display

###################   Plot ###########

mg.formatGraph()
plt.suptitle(modname[:-6].encode('string-escape'))

mg.lightcurve(t, y, yerr, btdb=True, mj=True, fmt='k.')
mg.lightcurve(t, y-d, yerr, btdb=True, mj=True, fmt='k.')
mg.lightcurve(x, a, btdb=True, mj=True, fmt='b')
mg.lightcurve(x, b, btdb=True, mj=True, fmt='g')
mg.lightcurve(x, c, btdb=True, mj=True, fmt='c')
mg.lightcurve(x, d, btdb=True, mj=True, fmt='r', ylow=-0.025, yhigh=0.225, xleft=np.median(time)-0.02, xright=np.median(time)+0.02)

print "Saving to", modname[:-4]+".pdf"

mg.saveAsTemp()
plt.savefig(modname[:-4]+".pdf")

plt.show()
