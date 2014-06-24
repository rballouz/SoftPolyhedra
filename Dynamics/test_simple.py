import matplotlib.pyplot as plt
from matplotlib.figure import *
import numpy as np
import pylab
data=np.genfromtxt("test.dat")
for i in np.arange(len(data[:,0])):
  plt.xlim(-20,20)
  plt.ylim(-20,20)
  plt.scatter(data[i,1],data[i,3],s=20,color='k',marker='o')
  plt.scatter(data[i,2],data[i,4],s=20,color='r',marker='o')


  pylab.savefig("test"+str(i).zfill(3)+".png")
  plt.clf()



