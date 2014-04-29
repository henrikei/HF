# -*- coding: utf-8 -*-
"""
Spyder Editor

This temporary script file is located here:
/home/henrik/.spyder2/.temp.py
"""

import numpy as np
import matplotlib.pyplot as plt


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

# Project density onto xy-plane
density = sum(density, axis=2)

# meshgrid
x = np.linspace(xstart, xend, nx)
y = np.linspace(ystart, yend, ny)
X, Y = np.meshgrid(x,y)

fig, ax = plt.subplots()
p = ax.contour(Y, X, density, 300, cmap = cm.RdBu, vmin = density.min(),
               vmax = density.max())

density = sum(density, axis=1)

fig2, ax2 = plt.subplots()
ax2.plot(x, density)