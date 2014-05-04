import numpy as np
from mayavi import mlab

A = np.loadtxt('min_config.dat')
    
atoms_x = A[:,0]
atoms_y = A[:,1]
atoms_z = A[:,2]

a1 = A[1,:]
a2 = A[2,:]
a3 = A[3,:]

f = 180.0/pi;
angle1 = np.arccos(np.dot(a1,a2)/(norm(a1)*norm(a2)))*f
angle2 = np.arccos(np.dot(a2,a3)/(norm(a2)*norm(a3)))*f
angle3 = np.arccos(np.dot(a3,a1)/(norm(a3)*norm(a1)))*f

angle1 = np.arccos(np.dot(a1,a2))


N = mlab.points3d(atoms_x[0:1], atoms_y[0:1], atoms_z[0:1],
                  scale_factor=0.4,
                  resolution=20,
                  color=(0, 0, 0),
                  scale_mode='none')

H = mlab.points3d(atoms_x[1:2], atoms_y[1:2], atoms_z[1:2],
                   scale_factor=0.4,
                   resolution=20,
                   color=(1, 0, 0),
                   scale_mode='none')

H = mlab.points3d(atoms_x[2:3], atoms_y[2:3], atoms_z[2:3],
                   scale_factor=0.4,
                   resolution=20,
                   color=(1, 0, 0),
                   scale_mode='none')
                   
H = mlab.points3d(atoms_x[3:4], atoms_y[3:4], atoms_z[3:4],
                   scale_factor=0.4,
                   resolution=20,
                   color=(1, 0, 0),
                   scale_mode='none')
bonds_x = zeros(6)                   
bonds_y = zeros(6)
bonds_z = zeros(6)

for i in range(6):
    if np.mod(i,2) == 0:
        bonds_x[i] = atoms_x[0]
        bonds_y[i] = atoms_y[0]
        bonds_z[i] = atoms_z[0]
    else:
        bonds_x[i] = atoms_x[1 + i/2]
        bonds_y[i] = atoms_y[1 + i/2]
        bonds_z[i] = atoms_z[1 + i/2]
        
                   
mlab.plot3d(bonds_x, bonds_y, bonds_z, [1, 2, 1, 2, 1, 2],
            tube_radius=0.1, colormap='Reds')
            
mlab.show()