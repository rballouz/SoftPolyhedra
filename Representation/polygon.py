import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


#Random Points from -20 to 20 with integer Coordinates
def square(x):
    return x*x
def mag(x,y,z):
    return np.sqrt(square(x)+square(y)+square(z))
N=50
lat=(np.random.rand(N)*180)*np.pi/180
lon=(np.random.rand(N)*360)*np.pi/180
Rad=(np.random.rand(N)*5)+95
x=Rad*np.sin(lat)*np.cos(lon)
y=Rad*np.sin(lat)*np.sin(lon)
z=Rad*np.cos(lat)

for i in np.arange(N):
    print x[i].astype(int),y[i].astype(int),z[i].astype(int)#,mag(x[i],y[i],z[i])

