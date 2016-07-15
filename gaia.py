
def prior(p):
	lp = 0.
	lp += ((p['length_spot']-0.037)/0.003)**2
	lp += ((p['yaw_spot']-45)/20)**2
	lp += ((p['angle_spot']-15)/20)**2
	if p['rdisc2'] > 0.7:
		lp += ((p['rdisc2']-0.7)/0.01)**2
	
	return -lp/2
