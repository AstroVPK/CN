import math
import cmath
import numpy as np
import scipy.sparse as sparse
import scipy.sparse.linalg as linalg
import os
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import pdb

import UserFuncs

#plt.ion()
pwd = os.environ['PWD']
cmapChoice = 'plasma'

x0 = 0.0
y0 = 0.0
z0 = 0.0
x1 = 1.0
y1 = 1.0
delta_l = 1.0e-2                                 # block_size (m)
delta_t = 1.0e-3                                 # timestep (s)
X = int((x1 - x0)/delta_l)  # Number of blocks in X
Y = int((y1 - y0)/delta_l)  # Number of blocks in Y
Z = 1                       # Number of blocks in Z
z1 = z0 + delta_l
T = 500  # Number of timesteps

nblocks = X*Y*Z

def offset(i, j, k, X, Y, Z):
	return i + j*X + k*X*Y

def ijk(offset, X, Y, Z):
	k = offset/(X*Y)
	offset -= k*(X*Y)
	j = offset/X
	offset -= j*X
	i = offset
	return i,j,k

def xyz(i, j, k, x0, y0, z0, delta_l):
	x = x0 + i*delta_l
	y = y0 + j*delta_l
	z = z0 + k*delta_l
	return x, y, z

def initAll(lhsM, rhsM, grid, X, Y, Z):
	for k in xrange(Z):
		for j in xrange(Y):
			for i in xrange(X):

				x, y, z = xyz(i, j, k, x0, y0, z0, delta_l)
				grid[offset(i, j, k, X, Y, Z)] = UserFuncs.T(x, y, z)
				R = UserFuncs.r(x, y, z)
				loc_Block = offset(i, j, k, X, Y, Z)
				lhsM[loc_Block, loc_Block] = 6*R + 1
				rhsM[loc_Block, loc_Block] = 1.0 - 6*R

				loc_Nbr = offset(i + 1, j, k, X, Y, Z)
				if 0 <= loc_Nbr < X*Y*Z:
					x, y, z = xyz(i + 1, j, k, x0, y0, z0, delta_l)
					R = UserFuncs.r(x, y, z)
					lhsM[loc_Block, loc_Nbr] = -1.0*R
					rhsM[loc_Block, loc_Nbr] = R
				loc_Nbr = offset(i - 1, j, k, X, Y, Z)
				if 0 <= loc_Nbr < X*Y*Z:
					x, y, z = xyz(i - 1, j, k, x0, y0, z0, delta_l)
					R = UserFuncs.r(x, y, z)
					lhsM[loc_Block, loc_Nbr] = -1.0*R
					rhsM[loc_Block, loc_Nbr] = R
				loc_Nbr = offset(i, j + 1, k, X, Y, Z)
				if 0 <= loc_Nbr < X*Y*Z:
					x, y, z = xyz(i, j + 1, k, x0, y0, z0, delta_l)
					R = UserFuncs.r(x, y, z)
					lhsM[loc_Block, loc_Nbr] = -1.0*R
					rhsM[loc_Block, loc_Nbr] = R
				loc_Nbr = offset(i , j - 1, k, X, Y, Z)
				if 0 <= loc_Nbr < X*Y*Z:
					x, y, z = xyz(i, j - 1, k, x0, y0, z0, delta_l)
					R = UserFuncs.r(x, y, z)
					lhsM[loc_Block, loc_Nbr] = -1.0*R
					rhsM[loc_Block, loc_Nbr] = R
				loc_Nbr = offset(i, j, k + 1, X, Y, Z)
				if 0 <= loc_Nbr < X*Y*Z:
					x, y, z = xyz(i, j, k + 1, x0, y0, z0, delta_l)
					R = UserFuncs.r(x, y, z)
					lhsM[loc_Block, loc_Nbr] = -1.0*R
					rhsM[loc_Block, loc_Nbr] = R
				loc_Nbr = offset(i, j, k - 1, X, Y, Z)
				if 0 <= loc_Nbr < X*Y*Z:
					x, y, z = xyz(i, j, k - 1, x0, y0, z0, delta_l)
					R = UserFuncs.r(x, y, z)
					lhsM[loc_Block, loc_Nbr] = -1.0*R
					rhsM[loc_Block, loc_Nbr] = R

def readGrid3D(grid, grid3D, X, Y, Z):
	vmin = np.min(grid)
	vmax = np.max(grid)
	for offset in xrange(grid.shape[0]):
		i, j, k = ijk(offset, X, Y, Z)
		grid3D[i,j,k] = grid[offset]
	return vmin, vmax

LHSMatrix = sparse.csr_matrix((nblocks, nblocks))
RHSMatrix = sparse.csr_matrix((nblocks, nblocks))
grid = np.zeros(X*Y*Z)
grid3D = np.zeros((X, Y, Z))
x = np.linspace(x0, x1 + delta_l, X + 1)
y = np.linspace(y0, y1 + delta_l, Y + 1)
z = np.linspace(z0, z1 + delta_l, Z + 1)
initAll(LHSMatrix, RHSMatrix, grid, X, Y, Z)
vmin, vmax = readGrid3D(grid, grid3D, X, Y, Z)
fig2 = plt.figure()
ims = []
ims.append((plt.pcolor(x, y, grid3D[:,:,0], vmin=vmin, vmax=vmax, cmap=cmapChoice),))

for t in xrange(T):
	RHSMatrix.multiply(grid)
	grid = linalg.spsolve(LHSMatrix, grid)
	readGrid3D(grid, grid3D, X, Y, Z)
	ims.append((plt.pcolor(x, y, grid3D[:,:,0], vmin=vmin, vmax=vmax, cmap=cmapChoice),))

im_ani = animation.ArtistAnimation(fig2, ims, interval=50, repeat_delay=3000, blit=True)
# To save this second animation with some metadata, use the following command:
im_ani.save('Si_Si3N4_Plane.mp4', metadata={'artist':'VPK'})
plt.show()


