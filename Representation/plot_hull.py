import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


points=np.genfromtxt("hull_3d.in")
hull=np.genfromtxt("hull_edges.txt")


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.scatter(points[:,0],points[:,1],points[:,2],color='b',marker='o')
for i in np.arange(hull.shape[0]):
    ax.plot(hull[i,(0,3)],hull[i,(1,4)],hull[i,(2,5)],color='r')
plt.show()


