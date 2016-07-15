#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from trm import dnl, subs, sla
from scipy.optimize import curve_fit
from sys import argv, exit
import stars, runs


sets = [runs.jan14, runs.jan15, runs.jan16, runs.jan17, runs.may, runs.june]
colours = ["r", "g", "b"]
number = "4"
modflags = np.ones(75, dtype=str)
#modflags[45] = '2'
#modflags[46] = '2'#'t'
#modflags[47] = '2'
#modflags[49] = '2'
#modflags[50] = 't'
#modflags[57] = '2'
#modflags[58] = 't'
#modflags[59] = 't'
#modflags[60] = 't'
#modflags[61] = '2'
#modflags[62] = 't'

i = 0
for set in sets:
	for colour in colours:
		for num in np.arange(set.numec)+1:
			sname = set.name + "-" + colour + "-" + str(num) + "_" + number + ".pbs"
			
			data = set.name + "-" + colour + "-" + str(num) + ".dat"
			
			
			
			
			model = set.name + "-" + colour + "-" + str(num) + "_" + modflags[i]+ ".mod"   #set+colour+".mod"
			i += 1
			prior = "gaia1.py"
			chain = set.name + "-" + colour + "-" + str(num) + "_" + number + ".log"
			append = False
			
			file = open(sname, 'w')
			file.write("#!/bin/bash\n#PBS -l nodes=1:ppn=4,pvmem=1000mb,walltime=50:00:00\n#PBS -V\n\ncd $PBS_O_WORKDIR\n\n./wdmcmc.py data=%s model=%s init=gaia.ini prior=%s chain=%s append=%s lmarq=gaia.lmq \\\\"%(data, model, prior, chain, append))
			file.close()
	
	
file = open("cow.bsh", 'w')
file.write("#!/bin/bash \n\n")
for set in sets:
	for colour in colours:
		for num in np.arange(set.numec)+1:
			sname = set.name + "-" + colour + "-" + str(num) + "_" + number + ".pbs"
			file.write("qsub %s \n"%(sname))
file.close()	
	
	
