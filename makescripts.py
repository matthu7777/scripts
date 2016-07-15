#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from trm import dnl, subs, sla
from scipy.optimize import curve_fit
from sys import argv, exit
import stars, runs


sets = [runs.jan14, runs.jan15, runs.jan16, runs.jan17, runs.may, runs.june]
colours = ["r", "g", "b"]
number = "8"

for set in sets:
	for colour in colours:
		for num in np.arange(set.numec)+1:
			sname = set.name + "-" + colour + "-" + str(num) + "_" + number + ".pbs"
			
			data = set.name + "-" + colour + "-" + str(num) + ".dat"
			
			
			
			
			model = "gaia2.mod"   #set+colour+".mod"
			prior = "gaia.py"
			chain = set.name + "-" + colour + "-" + str(num) + "_" + number + ".log"
			append = False
			
			file = open(sname, 'w')
			file.write("#!/bin/bash\n#PBS -l nodes=1:ppn=4,pvmem=1000mb,walltime=30:00:00\n#PBS -V\n\ncd $PBS_O_WORKDIR\n\nlcmcmc.py data=%s model=%s init=gaia.ini prior=%s chain=%s append=%s lmarq=gaia.lmq \\\\"%(data, model, prior, chain, append))
			file.close()
	
	
	
	
	
