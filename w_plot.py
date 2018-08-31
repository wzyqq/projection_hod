import numpy as np
import matplotlib.pyplot as plt

parallel_max = 5
perpendicular_max = 28

x1 = np.logspace(-1,1,21)
x2 = np.logspace(-1,1.7,28)

y1 = np.loadtxt('chain1_5964').reshape([3,21])
y2 = np.loadtxt('w_mean')

covariance_matrix  = np.loadtxt('covariance_matrix')
covariance_old     = covariance_matrix.reshape((parallel_max*perpendicular_max,parallel_max*perpendicular_max))
covariance 	       = covariance_old
bar = np.zeros(parallel_max*perpendicular_max)
for i in range(parallel_max*perpendicular_max):
	bar[i] = np.sqrt(covariance[i][i])
plt.errorbar(x2, y2[0*perpendicular_max:1*perpendicular_max], yerr=bar[0*perpendicular_max:1*perpendicular_max], capsize=2, marker='s', markersize=1, color='r', label='whole100')
plt.errorbar(x2, y2[1*perpendicular_max:2*perpendicular_max], yerr=bar[1*perpendicular_max:2*perpendicular_max], capsize=2, marker='s', markersize=1, color='g', label='whole200')
plt.errorbar(x2, y2[2*perpendicular_max:3*perpendicular_max], yerr=bar[2*perpendicular_max:3*perpendicular_max], capsize=2, marker='s', markersize=1, color='b', label='whole300')
#plt.errorbar(x2, y2[3*perpendicular_max:4*perpendicular_max], yerr=bar[3*perpendicular_max:4*perpendicular_max], capsize=2, marker='s', markersize=1, color='m', label='whole400')
#plt.errorbar(x2, y2[4*perpendicular_max:5*perpendicular_max], yerr=bar[4*perpendicular_max:5*perpendicular_max], capsize=2, marker='s', markersize=1, color='c', label='whole500')
plt.plot(x1,y1[0],'r',marker='o',label='bestfit-whole100')
plt.plot(x1,y1[1],'g',marker='o',label='bestfit-whole200')
plt.plot(x1,y1[2],'b',marker='o',label='bestfit-whole300')
plt.xlabel("$r_{p}$(Mpc/h)",fontsize=30)
plt.ylabel('$w_{p}$',fontsize=30)
plt.title('projected correlation functions',fontsize=30)
plt.xscale('log')
plt.yscale('log',nonposy='clip')
plt.legend()
plt.show()