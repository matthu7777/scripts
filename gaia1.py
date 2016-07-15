
def prior(p):
	lp = 0.
	lp += ((p['length_spot']-0.038)/0.002)**2
	lp += ((p['rdisc2']-0.54)/0.005)**2
	lp += ((p['radius_spot']-0.410)/0.004)**2
	lp += ((p['angle_spot']-15)/2)**2
	lp += ((p['yaw_spot']-45)/2)**2
	
	return -lp/2
