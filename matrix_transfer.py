import numpy as np
from numpy.linalg import inv
points = 21
parallel_max = 3
perpendicular_max = 28
covariance = np.zeros((points*parallel_max,points*parallel_max))
covariance_inverse = np.zeros((points*parallel_max,points*parallel_max))
covariance_old = np.loadtxt('covariance_matrix')
for i in range(points*parallel_max):
	for j in range(points*parallel_max):
		covariance[i][j] = covariance_old[int(i/points)*perpendicular_max+0+i%points,int(j/points)*perpendicular_max+0+j%points]
covariance_inverse = inv(covariance)
np.savetxt('inv_matrix21',covariance_inverse,fmt='%.8f')
