#!/usr/bin/env python

from subprocess import call
from os import listdir


for f in listdir("."):
	if f[-4:] == ".log":
		logname = f
		modname = f[:-4]+".mod"
		initname = "gaia1.mod"
		datname = f[:-6]+".dat"
		method = "p"
		
		call(["extbmod.py",initname,logname, modname, method])
		
		call(["splitter.csh", modname, datname])
		
		
