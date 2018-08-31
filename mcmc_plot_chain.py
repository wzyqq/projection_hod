import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(1,10000,10000)

data = np.loadtxt('foutput')

data = data.reshape(8,30,10000)

for i in range(30):
	plt.subplot(421)
	plt.plot(x,data[0,i])
	plt.title('$logM_{min}$')	
	plt.subplot(422)
	plt.plot(x,data[1,i])
	plt.title('$logM_{1}$')
	plt.subplot(423)
	plt.plot(x,data[2,i])
	plt.title('$\sigma$')
	plt.subplot(424)
	plt.title('$n_{g}$')
	plt.plot(x,data[3,i])
	plt.subplot(425)
	plt.plot(x,data[4,i])
	plt.title('$M_{h}$')
	plt.subplot(426)
	plt.plot(x,data[5,i])
	plt.title('$f_{sat}$')
	plt.subplot(427)
	plt.plot(x,data[6,i])
	plt.title('$b_{eff}$')	
	plt.subplot(428)
	plt.plot(x,data[7,i])
	plt.ylim(0,250)
	plt.title('${\chi}^{2}$')			
	plt.show()