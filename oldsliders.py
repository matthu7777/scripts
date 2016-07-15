#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
from os.path import isfile
import pickle

fname = "may-r-2.dat"
fwd = "a.out"

data = np.loadtxt(fname)
time = data[:,0]
flux = data[:,2]
err = data[:,3]

wdorig = np.loadtxt(fwd)[:,2]


period = 0.034519576212191858
t0 = 57166.116142215891

offset = 0.01
dlscale = 1

scale = 1



class Eclipse(object):
	def __init__(self, name):
		self.name = name

name = "may-r-1"
if isfile("eclipses/%s.pk1"%(name)):
	with open("eclipses/%s.pk1"%(name), 'rb') as f:
		eclipse = pickle.load(f)
	p1 = eclipse.p1
	p2 = eclipse.p2
	p3 = eclipse.p3
	p4 = eclipse.p4
	print "Loaded previous contact points"

else:
	p1 = 0
	p2 = 0
	p3 = 0
	p4 = 0

xmin = -0.1
xmax = 0.2
ymin = -0.04
ymax = 0.1


#time = data[:,1]
#counts = data[:,14]/data[:,28]
#errors = np.sqrt((data[:,15]/data[:,14])**2 + (data[:,29]/data[:,28])**2) * counts

phase = ((time-t0)/period) % 1
phase[phase>0.5] -= 1

wd = wdorig*scale

f = flux - wd


def differentiate(f, time):
	#dt = np.gradient(time)
	#df = np.gradient(f, dt)
	c = np.zeros(len(f))
	for i in xrange(1):
		c[i+1] += 1
		c[-i] -= 1
	#c = [-1,1]
	df = np.convolve(f, c)
	dfr = np.convolve(df, [1]*10)
	return df, dfr
df, dfr = differentiate(f, time)








fig, ax = plt.subplots()
plt.subplots_adjust(left=0.1, bottom=0.25, right=0.75)
plt.axis([xmin, xmax, ymin, ymax])

oi = (phase > xmin) & (phase < xmax)

l, = plt.plot(phase[oi], f[oi], 'b.')
dl, = plt.plot(phase[oi], df[oi]*dlscale - offset, 'g.')
dlr, = plt.plot(phase[oi], dfr[oi]*dlscale - offset, 'g-')

lp1, = plt.plot([p1,p1,], [ymin, ymax], 'k-')
lp2, = plt.plot([p2,p2,], [ymin, ymax], 'k-')
lp3, = plt.plot([p3,p3,], [ymin, ymax], 'k-')
lp4, = plt.plot([p4,p4,], [ymin, ymax], 'k-')

axwd = plt.axes([0.1,0.18,0.65,0.03])
axp1 = plt.axes([0.1, 0.13, 0.65, 0.03])
axp2 = plt.axes([0.1, 0.1, 0.65, 0.03])
axp3 = plt.axes([0.1, 0.07, 0.65, 0.03])
axp4 = plt.axes([0.1, 0.04, 0.65, 0.03])
axdl = plt.axes([0.8, 0.3, 0.1, 0.03])
axsc = plt.axes([0.8, 0.5, 0.1, 0.03])
axsp = plt.axes([0.8, 0.4, 0.1, 0.03])

swd = Slider(axwd, 'WD', 0, 2*np.max(wdorig), valinit=np.max(wdorig))
sp1 = Slider(axp1, 'P1', xmin, xmax, valinit=p1)
sp2 = Slider(axp2, 'P2', xmin, xmax, valinit=p2)
sp3 = Slider(axp3, 'P3', xmin, xmax, valinit=p3)
sp4 = Slider(axp4, 'P4', xmin, xmax, valinit=p4)

twd = fig.text(0.8,0.8,str(np.max(wdorig)))
tp1 = fig.text(0.8,0.75,str(p1))
tp2 = fig.text(0.8,0.7,str(p2))
tp3 = fig.text(0.8,0.65,str(np.max(p3)))
tp4 = fig.text(0.8,0.6,str(np.max(p4)))

def updatewd(val):
	valwd = swd.val
	scale = valwd / np.max(wdorig)
	wd = wdorig*scale
	global f
	f = flux-wd
	l.set_ydata(f[oi])
	#dl.set_ydata()
	twd.set_text(str(valwd))
	fig.canvas.draw_idle()

def updatep1(val):
	p1 = sp1.val
	lp1.set_xdata([p1,p1])
	tp1.set_text(str(p1))
	fig.canvas.draw_idle()

def updatep2(val):
	p2 = sp2.val
	lp2.set_xdata([p2,p2])
	tp2.set_text(str(p2))

def updatep3(val):
	p3 = sp3.val
	lp3.set_xdata([p3,p3])
	tp3.set_text(str(p3))

def updatep4(val):
	p4 = sp4.val
	lp4.set_xdata([p4,p4])
	tp4.set_text(str(p4))

swd.on_changed(updatewd)
sp1.on_changed(updatep1)
sp2.on_changed(updatep2)
sp3.on_changed(updatep3)
sp4.on_changed(updatep4)




button = Button(axdl, 'Differentiate', hovercolor='0.975')


def refreshdl(event):
    df, dfr = differentiate(f, time)
    dl.set_ydata(df[oi]*dlscale - offset)
    dlr.set_ydata(dfr[oi]*dlscale - offset)
    fig.canvas.draw_idle()
button.on_clicked(refreshdl)

bsavec = Button(axsc, 'Save contact points', hovercolor='0.975')
def savec(event):
	eclipse = Eclipse(name)
	eclipse.p1 = sp1.val
	eclipse.p2 = sp2.val
	eclipse.p3 = sp3.val
	eclipse.p4 = sp4.val
	with open("bseclipses/%s.pk1"%(name), 'wb') as f:
		pickle.dump(eclipse,f,pickle.HIGHEST_PROTOCOL)
bsavec.on_clicked(savec)

bsavep = Button(axsp, 'Save picture', hovercolor='0.975')
def savep(event):
	plt.savefig("figures/%s.pdf"%(name))
	plt.savefig("figures/%s.png"%(name))
bsavep.on_clicked(savep)

plt.show()
