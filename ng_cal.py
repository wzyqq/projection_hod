import numpy as np
subsample = 100
observeall=[]
for i in range(subsample):	
	observe_length  = "/home/wzy/Desktop/cfht_w1/jack_no_stellar/whole/whole_sample_position_length/whole_sample_position_length_%d"%(i)
	observelength = np.loadtxt(observe_length)
#	observe = np.array([len(observelength)])
#	for k in range(len(observe)):
	observeall.append(observelength)
observe = np.array(observeall)
ng = observe/(80284132361.29527*16.4992*87.0/100.0/41252.96125/100)
ng_average = np.mean(ng)
ng_std = np.std(ng)
#np.savetxt("ng_average",ng_average,fmt='%14.8lf')
#np.savetxt("ng_std",ng_std,fmt='%14.8lf')
