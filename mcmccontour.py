import numpy as np
import matplotlib.pyplot as plt
import corner

ndim = 7
data = np.loadtxt('foutput')
data = data.reshape(8,30,10000)
draw = data[:,:,4000:]
use = list(range(30))
use.remove(6)
#sample = (draw[0:ndim,[0,1,2,3,4,5,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29],:].reshape([ndim,-1])).T
sample = (draw[0:ndim,use,:].reshape([ndim,-1])).T

#value1 = np.array([50,1.0])
value2 = np.mean(sample, axis=0)

figure = corner.corner(sample,bins=[101,101,101,101,101,101,101],levels=(1-np.exp(-0.5),1-np.exp(-2),1-np.exp(-4)),
labels=["$logM_{min}$", "$logM_{1}$", "$\sigma$","$n_{g}$","$M_{h}$","$f_{sat}$","$b_{eff}$"],plot_contours=True,quantiles=[0.0013, 0.0228, 0.1587, 0.5, 0.8413, 0.9772, 0.997],
show_titles=True,title_kwargs={"fontsize": 15},label_kwargs={"fontsize": 20},title_fmt=".2f")
axes = np.array(figure.axes).reshape((ndim, ndim))
# Loop over the diagonal
for i in range(ndim):
    ax = axes[i, i]
#    ax.axvline(value1[i], color="g")
    ax.axvline(value2[i], color="r")

# Loop over the histograms
for yi in range(ndim):
    for xi in range(yi):
        ax = axes[yi, xi]
#        ax.axvline(value1[xi], color="g")
        ax.axvline(value2[xi], color="r")
#        ax.axhline(value1[yi], color="g")
        ax.axhline(value2[yi], color="r")
#        ax.plot(value1[xi], value1[yi], "sg")
        ax.plot(value2[xi], value2[yi], "sr")

#sigam_best = np.percentile(sample[:,0],50)
#sigam_less = np.percentile(sample[:,0],16)
#sigma_more = np.percentile(sample[:,0],84)
#bias_best = np.percentile(sample[:,1],50)
#bias_less = np.percentile(sample[:,1],16)
#bias_more = np.percentile(sample[:,1],84)
#print(sigam_best,sigam_less,sigma_more,bias_best,bias_less,bias_more)
plt.show()