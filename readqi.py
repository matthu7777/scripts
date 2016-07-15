#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from trm import dnl, subs, sla, roche
import stars, runs
import mgutils as mg
import re
from sys import argv

""" Script to find q and i from levmarq or mcmc outputs. Determines which from names of files. Saves to a text file.
    If file names are *.pbs.o*, uses levmarq, if *.log* uses mcmc, else doesn't know what to do.
	If lev, greps through levmarq output to find something of form q = ??? and equivalent. 
"""

#### File list -- list of files to read q and i from
if len(argv) > 1:
	files = argv[1]
else:
	files = "outfiles.lis"			# files containing levmarq output. Should include output of form q = ???

####  Check file names to find out whether levmarq or mcmc results
with open(files, 'r') as lis:
	line = lis.readline()
	if ".log" in line:
		mcmc = True
	elif ".pbs.o" in line:
		mcmc = False
	else:
		print "File types not recognised -- should be either .log for mcmc or .pbs.o* (cow output) for levmarq"
		exit()

##### Output file -- text file to save output to
if len(argv) > 2:
	outfile = argv[2]
elif mcmc:
	outfile = "qimcmc.dat"
else:
	outfile = "qilev.dat"
outfile = "qi.dat"		

qreg = re.compile("q\s\=\s.+")
ireg = re.compile("iangle\s\=\s.+")

nlist = []
numrows = 75 #82
qlist, qerrlist, ilist, ierrlist, philist, phierrlist,  minqlist, minqerrlist = np.zeros(numrows), np.zeros(numrows), np.zeros(numrows), np.zeros(numrows), np.zeros(numrows), np.zeros(numrows), np.zeros(numrows), np.zeros(numrows)

#def functionErrors(funct, a, aerr, b=None


j=0
with open(files, 'r') as lis:
	for fname in lis:
		fname = fname[:-1]		#trim the newline character
		## Read levmarq results
		if not mcmc:			#good old double negative. 
			print "Reading levmarq results"
			with open(fname, 'r') as f:
				text = f.read()
				
				r1 = qreg.search(text)
				r2 = ireg.search(text)
				
				if r1 is not None and r2 is not None:
					q = float(r1.group(0).split()[2])
					qerr = float(r1.group(0).split()[4])
					
					i = float(r2.group(0).split()[2])
					ierr = float(r2.group(0).split()[4])
				
					phi = roche.findphi(q,i)
					phierr = np.sqrt( (roche.findphi(q+qerr,i) - phi)**2 + (roche.findphi(q, i+ierr) - phi)**2 )
					
					minq = roche.findq(90, phi)
					minqerr = np.fabs(roche.findq(90, phi+phierr) - minq)
					
				else:
					q, qerr, i, ierr, phi, phierr, minq, minqerr = 0,0,0,0,0,0,0,0

		## Read mcmc results
		else:	# if mcmc
			print "Reading mcmc results"
			dataq = mg.readMcmcLog(fname, "q")
			q, qerr = dataq.mean(), dataq.std()

			datai = mg.readMcmcLog(fname, "iangle")
			i, ierr = datai.mean(), datai.std()

			dataphi, dataminq = np.array([]), np.array([])
			for tempq, tempi in zip(dataq, datai):

				tempphi = roche.findphi(tempq,tempi)
				tempminq = roche.findq(90,tempphi)

				dataphi = np.append(dataphi,tempphi)
				dataminq = np.append(dataminq,tempminq)

			phi, phierr = dataphi.mean(), dataphi.std()
			minq, minqerr = dataminq.mean(), dataminq.std()


		name = fname.split('_')[0]
		
		nlist += [name]
		qlist[j] = q
		qerrlist[j] = qerr
		ilist[j] = i
		ierrlist[j] = ierr
		philist[j] = phi
		phierrlist[j] = phierr
		minqlist[j] = minq
		minqerrlist[j] = minqerr
		
		j += 1
print nlist
print qlist
print ilist
print philist
print minqlist


out = np.column_stack((nlist, qlist, qerrlist, ilist, ierrlist, philist, phierrlist, minqlist, minqerrlist))
np.savetxt(outfile, out, fmt='%s')






