#Picking Random Points on a Sphere
#Method from NRiC, Section 21.5.1
#Originally by Marsaglia (1972),

import numpy as np
import sys

def square(x):
    return x*x
def mag(x,y,z):
    return np.sqrt(square(x)+square(y)+square(z))

i=0


if (len(sys.argv) != 2):
    print "Usage: python random_sphere.py number_vertices"
    exit()

N=int(sys.argv[1])
Rad=1


while i < N:
    u0=(np.random.rand()*2)-1
    u1=(np.random.rand()*2)-1
    if mag(u0,u1,0)<=1:
        x=Rad*2*u0*np.sqrt(1-square(u0)-square(u1))
        y=Rad*2*u1*np.sqrt(1-square(u0)-square(u1))
        z=Rad*(1-2*(square(u0)+square(u1)))
        print x,y,z
        i=i+1
