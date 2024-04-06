#############################################################################
#
# Simple analysis file for xraypol
#
#############################################################################

# import packages #
import sys
import h5py
import math
import pyharm
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.offsetbox as offsetbox
from matplotlib.colors import TwoSlopeNorm

#############################################################################
# plotting function #
def plot(x, y, z, phi_in, var, cmap1, cmap2):
  plt.clf()
  zmin = np.nanmin(z)
  zmax = np.nanmax(z)
  if(zmin*zmax < 0):
    norm = TwoSlopeNorm(vmin=zmin, vcenter=0, vmax=zmax)
    tick1 = np.linspace(zmin, 0, 5, endpoint=False)
    tick2 = np.linspace(0, zmax, 5, endpoint=True)
    tick = np.concatenate((tick1, tick2), axis=0)
    cmap = cmap1
  else:
    norm = mpl.colors.Normalize(vmin=zmin, vmax=zmax)
    tick = np.linspace(zmin, zmax, 10,  endpoint=True)
    cmap = cmap2
  fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(5*2, 4*2), sharex=True, \
                          sharey=True,gridspec_kw={'hspace': 0, 'wspace': 0})
  ((ax1, ax2), (ax3, ax4)) = axs

  ax1.contourf(x, y, z[:,:,0], 100, norm=norm, cmap=cmap)
  textstr = '\n'.join((r'$\phi$ = %.2f' % (phi_in[0,0,0]*180/math.pi),))
  props = dict(color = 'black', size = 15) 
  ob = offsetbox.AnchoredText(textstr, loc='lower right', prop=props)
  ax1.add_artist(ob) 
  ax1.contourf(-x, y, z[:,:,0+z.shape[2]//2], 100, norm=norm, cmap=cmap)
  textstr = '\n'.join((r'$\phi$ = %.2f' % (phi_in[0,0,0+z.shape[2]//2]*180/math.pi),))
  props = dict(color = 'black', size = 15) 
  ob = offsetbox.AnchoredText(textstr, loc='lower left', prop=props)
  ax1.add_artist(ob) 

  ax2.contourf(x, y, z[:,:,z.shape[2]//8], 100, norm=norm, cmap=cmap)
  textstr = '\n'.join((r'$\phi$ = %.2f' % (phi_in[0,0,z.shape[2]//8]*180/math.pi),))
  props = dict(color = 'black', size = 15) 
  ob = offsetbox.AnchoredText(textstr, loc='lower right', prop=props)
  ax2.add_artist(ob) 
  ax2.contourf(-x, y, z[:,:,z.shape[2]//8+z.shape[2]//2], 100, norm=norm, cmap=cmap)
  textstr = '\n'.join((r'$\phi$ = %.2f' % (phi_in[0,0,z.shape[2]//8+z.shape[2]//2]*180/math.pi),))
  props = dict(color = 'black', size = 15) 
  ob = offsetbox.AnchoredText(textstr, loc='lower left', prop=props)
  ax2.add_artist(ob) 

  ax3.contourf(x, y, z[:,:,z.shape[2]//4], 100, norm=norm, cmap=cmap)
  textstr = '\n'.join((r'$\phi$ = %.2f' % (phi_in[0,0,z.shape[2]//4]*180/math.pi),))
  props = dict(color = 'black', size = 15) 
  ob = offsetbox.AnchoredText(textstr, loc='lower right', prop=props)
  ax3.add_artist(ob) 
  ax3.contourf(-x, y, z[:,:,z.shape[2]//4+z.shape[2]//2], 100, norm=norm, cmap=cmap)
  textstr = '\n'.join((r'$\phi$ = %.2f' % (phi_in[0,0,z.shape[2]//4+z.shape[2]//2]*180/math.pi),))
  props = dict(color = 'black', size = 15) 
  ob = offsetbox.AnchoredText(textstr, loc='lower left', prop=props)
  ax3.add_artist(ob) 

  ax4.contourf(x, y, z[:,:,3*z.shape[2]//8], 100, norm=norm, cmap=cmap)
  textstr = '\n'.join((r'$\phi$ = %.2f' % (phi_in[0,0,3*z.shape[2]//8]*180/math.pi),))
  props = dict(color = 'black', size = 15) 
  ob = offsetbox.AnchoredText(textstr, loc='lower right', prop=props)
  ax4.add_artist(ob) 
  ax4.contourf(-x, y, z[:,:,3*z.shape[2]//8+z.shape[2]//2], 100, norm=norm, cmap=cmap)
  textstr = '\n'.join((r'$\phi$ = %.2f' % (phi_in[0,0,3*z.shape[2]//8+z.shape[2]//2]*180/math.pi),))
  props = dict(color = 'black', size = 15) 
  ob = offsetbox.AnchoredText(textstr, loc='lower left', prop=props)
  ax4.add_artist(ob) 

  fig.text(0.01, 0.5, r'$y$', va='center', rotation='vertical')
  fig.text(0.45, 0.02, r'$x$', va='center', rotation='horizontal')
  fig.text(0.45, 0.96, var, va='center', rotation='horizontal')
  cbar_ax = fig.add_axes([0.85, 0.15, 0.03, 0.7])
  cbar = mpl.colorbar.ColorbarBase(cbar_ax, cmap=cmap, norm=norm, ticks=tick)
  plt.subplots_adjust(
    top=0.94,
    bottom=0.1,
    left=0.11,
    right=0.84,
    hspace=0.2,
    wspace=0.2)
  plt.show()
  plt.close()

#############################################################################

# read data #
data = h5py.File(sys.argv[1],'r+')

# assign #
theta = data['theta0'][:]
stokes = data['stokes'][:].T

# assign #
isc_total = stokes[0,:]
qsc_total = stokes[1,:]
usc_total = stokes[2,:]

#############################################################################

# plot #
plt.clf()
fig, axs = plt.subplots(nrows=3,ncols=1,sharex=True,gridspec_kw={'hspace': 0, 'wspace': 0})
(ax1, ax2, ax3) = axs
ax1.plot(theta*180/math.pi, qsc_total/isc_total, label=r'$Q/I$')
ax1.plot(theta*180/math.pi,np.sqrt(usc_total**2 + qsc_total**2)/isc_total, label=r'$\sqrt{U^{2} + Q^{2}}/I$')
ax1.legend()
ax1.grid()
ax2.hlines(y = 90, xmin = theta.min()*180/math.pi, xmax = theta.max()*180/math.pi, color='black', linestyles='-.')
ax2.plot(theta*180/math.pi, 0.5*(math.pi - np.arctan(usc_total/qsc_total))*180/math.pi, label=r'EVPA')
ax2.set_ylim(0, 180)
ax2.legend()
ax2.grid()
ax3.hlines(y = 2, xmin = theta.min()*180/math.pi, xmax = theta.max()*180/math.pi, color='red', linestyles='-.', label='threshold')
ax3.plot(theta*180/math.pi, np.log10(isc_total))
ax3.set_ylim(0, 2)
ax3.legend()
ax3.grid()
fig.text(0.01, 0.5, r'$y$', va='center', rotation='vertical')
fig.text(0.45, 0.05, r'$\theta$', va='center', rotation='horizontal')
plt.show()
plt.close()

1/0

# load map #
ang_fac_isc = angfile['ang_fac_isc'][()].T
ang_fac_qsc = angfile['ang_fac_qsc'][()].T
ang_fac_usc = angfile['ang_fac_usc'][()].T
nx = ang_fac_isc.shape[0]
n3 = ang_fac_isc.shape[2]

#########################################################################
# open HARM snapshot #
input = 'torus.h5'
dump = pyharm.load_dump(input)

# load harm #
dx1 = dump['dx1']
dx2 = dump['dx2']
dx3 = dump['dx3']
r = dump['r'][0:nx, :]  
th = dump['th'][0:nx, :]  
phi = dump['phi'][0:nx, :, :]  
rho = dump['rho'][0:nx, :, :]  
gdet = dump['gdet'][0:nx, :,]  
sigma = dump['sigma'][0:nx, :, :]  
vol = gdet*dx1*dx2*dx3
r = np.repeat(r, n3, axis=2)
th = np.repeat(th, n3, axis=2)
vol = np.repeat(vol, n3, axis=2)

#mesh grid#
X_half = r[:,:,0] * np.sin(th[:,:,0])
Z_half = r[:,:,0] * np.cos(th[:,:,0])

# plot #
z = ang_fac_isc
plot(X_half, Z_half, z, phi, '$_{e}$', 'RdGy', 'magma')
z = ang_fac_qsc
plot(X_half, Z_half, z, phi, '$n_{e}$', 'RdGy', 'magma')
z = ang_fac_usc
plot(X_half, Z_half, z, phi, '$n_{e}$', 'RdGy', 'magma')
1/0

#############################################################################

# assign #
r_grid = data['r-dir'][:].T
th_grid = data['th-dir'][:].T
phi_grid = data['phi-dir'][:].T

# assign #
nebar = data['electron'][:].T
vol = data['volume'][:].T
gam_fac = data['gam_fac'][:].T

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

#mesh grid#
X_half = r_grid[:,:,0] * np.sin(th_grid[:,:,0])
Z_half = r_grid[:,:,0] * np.cos(th_grid[:,:,0])

#############################################################################

# plot #
z = np.log10(nebar+np.max(nebar)/10**10) 
plot(X_half, Z_half, z, phi_grid, '$n_{e}$', 'RdGy', 'plasma')

# plot #
z = nebar*vol
plot(X_half, Z_half, z, phi_grid, '$n_{e}\Delta V$', 'RdGy', 'plasma')

# plot #
z = ang_isc_min*nebar*vol
plot(X_half, Z_half, z, phi_grid, '$f(\\theta, \phi)n_{e}\Delta V, PD = PD_{min}$', 'bwr', 'magma')

# plot #
z = ang_qsc_min*nebar*vol
plot(X_half, Z_half, z, phi_grid, '$g(\\theta, \phi)n_{e}\Delta V, PD = PD_{min}$', 'coolwarm', 'viridis')


# plot #
z = ang_usc_min*nebar*vol
plot(X_half, Z_half, z, phi_grid, '$h(\\theta, \phi)n_{e}\Delta V, PD = PD_{min}$', 'seismic', 'inferno')

# plot #
z = ang_isc_max*nebar*vol
plot(X_half, Z_half, z, phi_grid, '$f(\\theta, \phi)n_{e}\Delta V, PD = PD_{max}$', 'bwr', 'magma')

# plot #
z = ang_qsc_max*nebar*vol
plot(X_half, Z_half, z, phi_grid, '$g(\\theta, \phi)n_{e}\Delta V, PD = PD_{max}$', 'coolwarm', 'viridis')


# plot #
z = ang_usc_max*nebar*vol
plot(X_half, Z_half, z, phi_grid, '$h(\\theta, \phi)n_{e}\Delta V, PD = PD_{max}$', 'seismic', 'inferno')


#############################################################################