import numpy as np
from numpy import *
import sys

def square(x):
    return x*x
def mag(x,y,z):
    return np.sqrt(square(x)+square(y)+square(z))

#Pre-amble text
print "#include \"colors.inc\""
print "#include \"shapes.inc\""
print "#include \"polyhedron.inc\""
print "#include \"tagsam.inc\""

print "camera {"
print "\t location <0, 20, 10>"
print "\t look_at <0, 0, 10>"
print "}"
print "light_source { <40, 40, 20> color White shadowless}"
print "object {tagsam"
print "\t translate <0,0,5>"
print "\t }"

#List of polyhedron types
list=("polyfaces_006","polyfaces_007","polyfaces_008","polyfaces_009","polyfaces_010","polyfaces_011","polyfaces_012")

if (len(sys.argv) != 2):
    print "Usage: python poly_generate.py number_polys"
    exit()


N=int(sys.argv[1])
R=20
i=0

#Generate N polyhedron within bounds of sphere of Radius R
"""
while i < N:
    x=np.random.rand()*R-(R/2)
    y=np.random.rand()*R-(R/2)
    z=np.random.rand()*R-(R/2)
    index=int(np.random.rand()*(len(list)-1))

    if mag(x,y,z) <= 10:
        rx=np.random.rand()*360
        ry=np.random.rand()*360
        rz=np.random.rand()*360

        print "object {"+list[index]
        print "\t rotate <"+str(rx)+","+str(ry)+","+str(rz)+">"
        print "\t translate <"+str(x)+","+str(y)+","+str(z)+">"
        print "\t }"
        print " "

    i=i+1
"""
"""
#Generate shapes in a line
x=2*np.arange(len(list))-5
for i in np.arange(len(list)):
    
    print "object {"+list[i]
    print "\t translate x*"+str(x[i])
    print "\t }"
    print " "
"""
#Generate N polyhedron within bounds of sphere of Radius R

while i < N:
    x=np.random.rand()*3*R-(1.5*R)
    y=np.random.rand()*2*R-(R)
    z=np.random.rand()*0.5*R-0.25*R
    index=int(np.random.rand()*(len(list)-1))

    rx=np.random.rand()*360
    ry=np.random.rand()*360
    rz=np.random.rand()*360
    s=np.random.rand()*0.4+0.1
    if (z<4.2):
        print "object {"+list[index]
        print "\t rotate <"+str(rx)+","+str(ry)+","+str(rz)+">"
        print "\t scale "+str(s)
        print "\t translate <"+str(x)+","+str(y)+","+str(z)+">"

        print "\t }"
        print " "

    i=i+1







