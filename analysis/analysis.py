#############################################################################
#
# Simple analysis file for xraypol
#
#############################################################################

# import packages #
import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#############################################################################

# read data #
data = h5py.File("data.hdf5",'r+')    

# assign #
theta = data['theta0'][:]
stokes = data['stokes'][:].T

# assign #
isc_total = stokes[0,:]
qsc_total = stokes[1,:]
usc_total = stokes[2,:]

# assign #
r_grid = data['r-dir'][:].T
th_grid = data['th-dir'][:].T
phi_grid = data['phi-dir'][:].T

# assign #
nebar = data['electron'][:].T
vol = data['volume'][:].T
gammam1 = data['gammam1'][:].T

# assign #
angular = data['angular'][:].T
ang_isc_min = angular[0,0,:,:,:]
ang_qsc_min = angular[0,1,:,:,:]
ang_usc_min = angular[0,2,:,:,:]
ang_isc_max = angular[1,0,:,:,:]
ang_qsc_max = angular[1,1,:,:,:]
ang_usc_max = angular[1,2,:,:,:]

# shape #
nx = nebar.shape[0]
ny = nebar.shape[1]
nz = nebar.shape[2]

#############################################################################

# plot #
plt.clf()
plt.plot(theta, qsc_total/isc_total)
plt.grid()
plt.show()

#############################################################################

#mesh grid#
X_half = r_grid[:,:,0] * np.sin(th_grid[:,:,0])
Z_half = r_grid[:,:,0] * np.cos(th_grid[:,:,0])

# plot #
plt.clf()
z = nebar*vol
zmin = z.min()
zmax = z.max()
plt.contourf(X_half, Z_half, z[:,:,0], np.linspace(zmin, zmax, 100), cmap='plasma')
plt.contourf(-X_half, Z_half, z[:,:,nz//2], np.linspace(zmin, zmax, 100), cmap='plasma')
plt.colorbar()
plt.show()

# plot #
plt.clf()
z = ang_isc_min*nebar*vol
zmin = z.min()
zmax = z.max()
plt.contourf(X_half, Z_half, z[:,:,0], np.linspace(zmin, zmax, 100), cmap='magma')
plt.contourf(-X_half, Z_half, z[:,:,nz//2], np.linspace(zmin, zmax, 100), cmap='magma')
plt.colorbar()
plt.show()

# plot #
plt.clf()
z = ang_qsc_min*nebar*vol
zmin = z.min()
zmax = z.max()
plt.contourf(X_half, Z_half, z[:,:,0], np.linspace(zmin, zmax, 100), cmap='viridis')
plt.contourf(-X_half, Z_half, z[:,:,nz//2], np.linspace(zmin, zmax, 100), cmap='viridis')
plt.colorbar()
plt.show()

# plot #
plt.clf()
z = ang_usc_min*nebar*vol
zmin = z.min()
zmax = z.max()
plt.contourf(X_half, Z_half, z[:,:,0], np.linspace(zmin, zmax, 100), cmap='inferno')
plt.contourf(-X_half, Z_half, z[:,:,nz//2], np.linspace(zmin, zmax, 100), cmap='inferno')
plt.colorbar()
plt.show()

#############################################################################