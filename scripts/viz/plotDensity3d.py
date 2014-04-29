# -*- coding: utf-8 -*-
"""
Spyder Editor

This temporary script file is located here:
/home/henrik/.spyder2/.temp.py
"""

import numpy as np
from mayavi import mlab


# Read dimension (shape) of density array
f = open('density.dat','r')
xstart = float(f.readline())
ystart = float(f.readline())
zstart = float(f.readline())
xend = float(f.readline())
yend = float(f.readline())
zend = float(f.readline())
dx = float(f.readline())
dy = float(f.readline())
dz = float(f.readline())
nx = float(f.readline())
ny = float(f.readline())
nz = float(f.readline())

# Read array from file and reshape
density= np.loadtxt('density.dat', skiprows=12)
density = np.reshape(density, (nx,ny,nz))

# meshgrid
#x = np.linspace(xstart, xend, nx)
#y = np.linspace(ystart, yend, ny)
#z = np.linspace(zstart, zend, nz)
#X, Y, Z = np.meshgrid(x,y,z)

x = np.linspace(xstart, xend, nx)
y = np.linspace(ystart, yend, ny)
z = np.linspace(zstart, zend, nz)
X, Y, Z = np.mgrid[xstart :xend:  complex(nx), ystart: yend:complex(ny) ,zstart:zend:complex(nz)]

mlab.contour3d(X,Y,Z,log10(density))
mlab.axes()
mlab.show()

#fig, ax = plt.subplots()
#p = ax.contour(Y, X, density, 300, cmap = cm.RdBu, vmin = density.min(),
#               vmax = density.max())
#
#density = sum(density, axis=1)
#
#fig2, ax2 = plt.subplots()
#ax2.plot(x, density)