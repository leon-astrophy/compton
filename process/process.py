########################################################
#
# Pre-process rad GRMHD snapshots to fed into xraypol
#
########################################################

#import#
import sys
import h5py
import math
import pyharm
import numpy as np
import matplotlib.pyplot as plt

########################################################

###############################################
# convert 4-vector from FMKS to KS coordinate #
def ummks2uks(x1_in, x2_in, x3_in, u0_in, u1_in, u2_in, u3_in):
          
  # find derivative
  n_fac = 0.5*math.pi*((poly_alpha+1)*poly_xt**poly_alpha)/((poly_alpha+1)*poly_xt**poly_alpha+1)
  dthgdx2 = math.pi*(1 + (1 - hslope)*np.cos(2*math.pi*x2_in))
  dthjdx2 = 2*n_fac*(1 + ((2*x2_in - 1)/(poly_xt))**poly_alpha)
  drdx1 = np.exp(x1_in)
  dthdx1 = -(mks_smooth)*np.exp(-mks_smooth*(x1_in - startx1))
  dthdx2 = dthgdx2 + np.exp(-mks_smooth*(x1_in - startx1))*(dthjdx2 - dthgdx2)

  # convert r and theta component of 4-vector #
  u1_out = drdx1*u1_in
  u2_out = dthdx1*u1_in+dthdx2*u2_in

  # phi doesn't change #
  u3_out = u3_in

  # time component doesn't change #
  u0_out = u0_in

  # return #
  return u0_out, u1_out, u2_out, u3_out 

#############################################
# convert 4-vector from KS to BL coordinate #
def uks2ubl(x1_in, a_in, u0_in, u1_in, u2_in, u3_in):
  r_in = np.exp(x1_in)
  delta=r_in**2-2*r_in+a_in**2
  u0_out = u0_in
  u1_out = u1_in
  u2_out = u2_in
  u3_out = u3_in
  u0_out = u0_in-2*r_in/delta*u1_in
  u3_out = u3_in-a_in/delta*u1_in
  return u0_out, u1_out, u2_out, u3_out 

################################################
# convert 4-vector from BL to LLNRF coordinate #
def ubl2lnrf(r_in, th_in, a_in, u0_in, u1_in, u2_in, u3_in):

  # assign #
  mu=np.cos(th_in)
  D_BL=r_in*r_in-2*r_in+a_in*a_in
  AR=(r_in*r_in+a_in*a_in)**2-a_in*a_in*D_BL*(1.-mu*mu)
  RHO_BL=r_in*r_in+a_in*a_in*mu*mu
  ENU=np.sqrt(D_BL*RHO_BL/AR)
  EMU1=np.sqrt(RHO_BL/D_BL)
  EMU2=np.sqrt(RHO_BL)
  EPSI=np.sqrt(1.-mu*mu)*np.sqrt(AR/RHO_BL)
  OM_BL=2*a_in*r_in/AR

  # convert #
  u0_out = ENU*u0_in
  u1_out = EMU1*u1_in
  u2_out = EMU2*u2_in
  u3_out = EPSI*(u3_in-OM_BL*u0_in)

  # zero out 4-vector in the horizon #
  u0_out[np.where(D_BL <= 0)] = 1
  u1_out[np.where(D_BL <= 0)] = 0
  u2_out[np.where(D_BL <= 0)] = 0
  u3_out[np.where(D_BL <= 0)] = 0

  # output #
  return u0_out, u1_out, u2_out, u3_out

######################################################

# terminal input #
files = sys.argv[1]
outfile = sys.argv[2]
rad = True
plot = False
rad_lim = 20

######################################################
# read #
dump = pyharm.load_dump(files)
dfile = h5py.File(files, "r")

# get x1, x2, x3 #
startx1 = dump['startx1']
startx2 = dump['startx2']
startx3 = dump['startx3']

# n1, n2, n3 #
n1 = dump['n1']
n2 = dump['n2']
n3 = dump['n3']

# dx1, dx2, dx3 #
dx1 = dump['dx1']
dx2 = dump['dx2']
dx3 = dump['dx3']

# FMKS parameters #
a_sim = dfile['header/geom/mmks/a'][()]
hslope = dfile['header/geom/mmks/hslope'][()]
mks_smooth = dfile['header/geom/mmks/mks_smooth'][()]
poly_alpha = dfile['header/geom/mmks/poly_alpha'][()]
poly_xt = dfile['header/geom/mmks/poly_xt'][()]

# load velocity #
ucon = dump['ucon']
ucon_ks = ucon.copy()
ucon_bl = ucon.copy()
ucon_lnrf = ucon.copy()

# load r and theta #
r = dump['r']
th = dump['th']
r = np.repeat(r, n3, axis=2)
th = np.repeat(th, n3, axis=2)

# load x and z #
if(plot):
  x = dump['x']
  z = dump['z']
  for i in range (0, n1):
    if(r[i,0,0] > rad_lim):
      i_cut = i + 1
      break

######################################################
# set up x1, x2, and x3 #
x1 = np.ndarray(shape=(n1, n2, n3), dtype=float)
x2 = np.ndarray(shape=(n1, n2, n3), dtype=float)
x3 = np.ndarray(shape=(n1, n2, n3), dtype=float)

# get x1, x2, x3 #
for i in range (0, n1):
  x1[i,:,:] = startx1 + dx1*i
for i in range (0, n2):
  x2[:,i,:] = startx2 + dx2*i
for i in range (0, n3):
  x3[:,:,i] = startx3 + dx3*i

######################################################
# plot and check #
if(plot):
  plt.contourf(x[0:i_cut,:,0], z[0:i_cut,:,0], ucon[1,0:i_cut,:,0]/ucon[0,0:i_cut,:,0], 100)
  plt.contourf(x[0:i_cut,:,n3//2], -z[0:i_cut,:,n3//2], ucon[1,0:i_cut,:,n3//2]/ucon[0,0:i_cut,:,n3//2], 100)
  plt.colorbar()
  plt.title(r'$U^{r}_{\rm mmks}/U^{t}_{\rm mmks}$')
  plt.savefig('uroverut_mmks.png')
  plt.clf()

  plt.contourf(x[0:i_cut,:,0], z[0:i_cut,:,0], ucon[2,0:i_cut,:,0]/ucon[0,0:i_cut,:,0], 100)
  plt.contourf(x[0:i_cut,:,n3//2], -z[0:i_cut,:,n3//2], ucon[2,0:i_cut,:,n3//2]/ucon[0,0:i_cut,:,n3//2], 100)
  plt.colorbar()
  plt.title(r'$U^{\theta}_{\rm mmks}/U^{t}_{\rm mmks}$')
  plt.savefig('uthoverut_mmks.png')
  plt.clf()

  plt.contourf(x[0:i_cut,:,0], z[0:i_cut,:,0], ucon[3,0:i_cut,:,0]/ucon[0,0:i_cut,:,0], 100)
  plt.contourf(x[0:i_cut,:,n3//2], -z[0:i_cut,:,n3//2], ucon[3,0:i_cut,:,n3//2]/ucon[0,0:i_cut,:,n3//2], 100)
  plt.colorbar()
  plt.title(r'$U^{\phi}_{\rm mmks}/U^{t}_{\rm mmks}$')
  plt.savefig('uphioverut_mmks.png')
  plt.clf()

######################################################

# convert 4-velocity to ks #
ucon_ks[0,:,:,:], ucon_ks[1,:,:,:] ,ucon_ks[2,:,:,:], ucon_ks[3,:,:,:] =\
ummks2uks(x1[:,:,:], x2[:,:,:], x3[:,:,:], ucon[0,:,:,:], ucon[1,:,:,:] ,ucon[2,:,:,:], ucon[3,:,:,:])

######################################################
# plot and check #
if(plot):
  plt.contourf(x[0:i_cut,:,0], z[0:i_cut,:,0], ucon_ks[1,0:i_cut,:,0]/ucon_ks[0,0:i_cut,:,0], 100)
  plt.contourf(x[0:i_cut,:,n3//2], -z[0:i_cut,:,n3//2], ucon_ks[1,0:i_cut,:,n3//2]/ucon_ks[0,0:i_cut,:,n3//2], 100)
  plt.colorbar()
  plt.title(r'$U^{r}_{\rm ks}/U^{t}_{\rm ks}$')
  plt.savefig('uroverut_ks.png')
  plt.clf()

  plt.contourf(x[0:i_cut,:,0], z[0:i_cut,:,0], ucon_ks[2,0:i_cut,:,0]/ucon_ks[0,0:i_cut,:,0], 100)
  plt.contourf(x[0:i_cut,:,n3//2], -z[0:i_cut,:,n3//2], ucon_ks[2,0:i_cut,:,n3//2]/ucon_ks[0,0:i_cut,:,n3//2], 100)
  plt.colorbar()
  plt.title(r'$U^{\theta}_{\rm ks}/U^{t}_{\rm ks}$')
  plt.savefig('uthoverut_ks.png')
  plt.clf()

  plt.contourf(x[0:i_cut,:,0], z[0:i_cut,:,0], ucon_ks[3,0:i_cut,:,0]/ucon_ks[0,0:i_cut,:,0], 100)
  plt.contourf(x[0:i_cut,:,n3//2], -z[0:i_cut,:,n3//2], ucon_ks[3,0:i_cut,:,n3//2]/ucon_ks[0,0:i_cut,:,n3//2], 100)
  plt.colorbar()
  plt.title(r'$U^{\phi}_{\rm ks}/U^{t}_{\rm ks}$')
  plt.savefig('uphioverut_ks.png')
  plt.clf()

######################################################

# convert 4-velocity to bl #
ucon_bl[0,:,:,:], ucon_bl[1,:,:,:] ,ucon_bl[2,:,:,:], ucon_bl[3,:,:,:] =\
uks2ubl(x1[:,:,:], a_sim, ucon_ks[0,:,:,:], ucon_ks[1,:,:,:] ,ucon_ks[2,:,:,:], ucon_ks[3,:,:,:])

######################################################
# plot and check #
if(plot):
  plt.contourf(x[0:i_cut,:,0], z[0:i_cut,:,0], ucon_bl[1,0:i_cut,:,0]/ucon_bl[0,0:i_cut,:,0], 100)
  plt.contourf(x[0:i_cut,:,n3//2], -z[0:i_cut,:,n3//2], ucon_bl[1,0:i_cut,:,n3//2]/ucon_bl[0,0:i_cut,:,n3//2], 100)
  plt.colorbar()
  plt.title(r'$U^{r}_{\rm bl}/U^{t}_{\rm bl}$')
  plt.savefig('uroverut_bl.png')
  plt.clf()

  plt.contourf(x[0:i_cut,:,0], z[0:i_cut,:,0], ucon_bl[2,0:i_cut,:,0]/ucon_bl[0,0:i_cut,:,0], 100)
  plt.contourf(x[0:i_cut,:,n3//2], -z[0:i_cut,:,n3//2], ucon_bl[2,0:i_cut,:,n3//2]/ucon_bl[0,0:i_cut,:,n3//2], 100)
  plt.colorbar()
  plt.title(r'$U^{\theta}_{\rm bl}/U^{t}_{\rm bl}$')
  plt.savefig('uthoverut_bl.png')
  plt.clf()

  plt.contourf(x[0:i_cut,:,0], z[0:i_cut,:,0], ucon_bl[3,0:i_cut,:,0]/ucon_bl[0,0:i_cut,:,0], 100)
  plt.contourf(x[0:i_cut,:,n3//2], -z[0:i_cut,:,n3//2], ucon_bl[3,0:i_cut,:,n3//2]/ucon_bl[0,0:i_cut,:,n3//2], 100)
  plt.colorbar()
  plt.title(r'$U^{\phi}_{\rm bl}/U^{t}_{\rm bl}$')
  plt.savefig('uphioverut_bl.png')
  plt.clf()

######################################################

# free array #
ucon = []
ucon_ks = []
x1 = []
x2 = []
x3 = []

######################################################

# set up LNRF 3-velocity #
lnrf_vel = np.ndarray(shape=(3, n1, n2, n3), dtype=float)

# convert 4-velocity to lnrf #
ucon_lnrf[0,:,:,:], ucon_lnrf[1,:,:,:] ,ucon_lnrf[2,:,:,:], ucon_lnrf[3,:,:,:] =\
ubl2lnrf(r[:,:,:], th[:,:,:], a_sim, ucon_bl[0,:,:,:], ucon_bl[1,:,:,:] ,ucon_bl[2,:,:,:], ucon_bl[3,:,:,:])

######################################################
# plot and check #
if(plot):
  plt.contourf(x[0:i_cut,:,0], z[0:i_cut,:,0], ucon_lnrf[1,0:i_cut,:,0]/ucon_lnrf[0,0:i_cut,:,0], 100)
  plt.contourf(x[0:i_cut,:,n3//2], -z[0:i_cut,:,n3//2], ucon_lnrf[1,0:i_cut,:,n3//2]/ucon_lnrf[0,0:i_cut,:,n3//2], 100)
  plt.colorbar()
  plt.title(r'$U^{r}_{\rm lnrf}/U^{t}_{\rm lnrf}$')
  plt.savefig('uroverut_lnrf.png')
  plt.clf()

  plt.contourf(x[0:i_cut,:,0], z[0:i_cut,:,0], ucon_lnrf[2,0:i_cut,:,0]/ucon_lnrf[0,0:i_cut,:,0], 100)
  plt.contourf(x[0:i_cut,:,n3//2], -z[0:i_cut,:,n3//2], ucon_lnrf[2,0:i_cut,:,n3//2]/ucon_lnrf[0,0:i_cut,:,n3//2], 100)
  plt.colorbar()
  plt.title(r'$U^{\theta}_{\rm lnrf}/U^{t}_{\rm lnrf}$')
  plt.savefig('uthoverut_lnrf.png')
  plt.clf()

  plt.contourf(x[0:i_cut,:,0], z[0:i_cut,:,0], ucon_lnrf[3,0:i_cut,:,0]/ucon_lnrf[0,0:i_cut,:,0], 100)
  plt.contourf(x[0:i_cut,:,n3//2], -z[0:i_cut,:,n3//2], ucon_lnrf[3,0:i_cut,:,n3//2]/ucon_lnrf[0,0:i_cut,:,n3//2], 100)
  plt.colorbar()
  plt.title(r'$U^{\phi}_{\rm lnrf}/U^{t}_{\rm lnrf}$')
  plt.savefig('uphioverut_lnrf.png')
  plt.clf()

######################################################

# assign #
lnrf_vel[0,:,:,:] = ucon_lnrf[1,:,:,:]/ucon_lnrf[0,:,:,:]
lnrf_vel[1,:,:,:] = ucon_lnrf[2,:,:,:]/ucon_lnrf[0,:,:,:]
lnrf_vel[2,:,:,:] = ucon_lnrf[3,:,:,:]/ucon_lnrf[0,:,:,:]

######################################################s

#free array#
ucon_bl = []
ucon_lnrf = []

# get vsquare #
gamma = np.ndarray(shape=(n1, n2, n3), dtype=float)
gamma = lnrf_vel[0,:,:,:]**2 + lnrf_vel[1,:,:,:]**2 + lnrf_vel[2,:,:,:]**2

######################################################
# plot and check #
if(plot):
  plt.contourf(x[0:i_cut,:,0], z[0:i_cut,:,0], gamma[0:i_cut,:,0], 100)
  plt.contourf(x[0:i_cut,:,n3//2], -z[0:i_cut,:,n3//2], gamma[0:i_cut,:,n3//2], 100)
  plt.colorbar()
  plt.title(r'$|v|^{2}_{\rm lnrf}$')
  plt.savefig('vsquare_lnrf.png')
  plt.clf()

######################################################

# get gamma #
gamma = 1/np.sqrt(1 - gamma)

######################################################
# plot and check #
if(plot):
  plt.contourf(x[0:i_cut,:,0], z[0:i_cut,:,0], gamma[0:i_cut,:,0], 100)
  plt.contourf(x[0:i_cut,:,n3//2], -z[0:i_cut,:,n3//2], gamma[0:i_cut,:,n3//2], 100)
  plt.colorbar()
  plt.title(r'$\Gamma_{\rm lnrf}$')
  plt.savefig('gamma_lnrf.png')
  plt.clf()

######################################################
# load coordinate #
phi = dump['phi']

# load variables #
rho = dump['rho']
gdet = dump['gdet']
sigma = dump['sigma']
gamma_mmks = dump['Gamma']

######################################################
# plot and check #
if(plot):
  plt.contourf(x[0:i_cut,:,0], z[0:i_cut,:,0], gamma_mmks[0:i_cut,:,0], 100)
  plt.contourf(x[0:i_cut,:,n3//2], -z[0:i_cut,:,n3//2], gamma_mmks[0:i_cut,:,n3//2], 100)
  plt.colorbar()
  plt.title(r'$\Gamma_{\rm mmks}$')
  plt.savefig('gamma_mmks.png')
  plt.clf()
#gamma = gamma_mmks # remove this line
gamma_mmks = []

######################################################

# volume #
vol = gdet*dx1*dx2*dx3

# repeat #
vol = np.repeat(vol, n3, axis=2)

#length#
length = vol**(1/3)

######################################################
# photon energy density #
if(rad):
  rmunu = dfile['Rmunu'][()]
  uphoton = np.ndarray(shape=(rmunu.shape[0], rmunu.shape[1], rmunu.shape[2]), dtype=float)
  ucon = dump['ucon']
  ucov = dump['ucov']

  # loop #
  uphoton[:,:,:] = 0
  for mu in range (0, rmunu.shape[3]):
    for nu in range (0, rmunu.shape[4]):
      uphoton[:,:,:]  = uphoton[:,:,:] + rmunu[:,:,:,mu,nu]*ucov[mu,:,:,:]*ucon[nu,:,:,:]

  # free #
  rmunu = []
  ucov = []
  ucon = []

######################################################
# get electron directional velocity #
vabs = np.sqrt(lnrf_vel[0,:,:,:]**2 + lnrf_vel[1,:,:,:]**2 + lnrf_vel[2,:,:,:]**2)
nhhat_r = lnrf_vel[0,:,:,:]/vabs
nhhat_th = lnrf_vel[1,:,:,:]/vabs
nhhat_phi = lnrf_vel[2,:,:,:]/vabs

#free#
vabs = []
lnrf_vel = []

######################################################
# intensity #
if(rad):
  intensity = uphoton
else:
  intensity = rho.copy()
  intensity[:,:,:] = 1

# scale down the variables #
rho = rho/rho.max()
intensity = intensity/intensity.max()
vol = vol/vol.max()
length = length/length.max()

######################################################

# output #
hf = h5py.File(outfile, 'w')

# reassign #
n1 = rho.shape[0]
n2 = rho.shape[1]
n3 = rho.shape[2]

# create data #
hf.create_dataset('n1',data=n1)
hf.create_dataset('n2',data=n2)
hf.create_dataset('n3',data=n3)
hf.create_dataset('r',data=r.T)
hf.create_dataset('th',data=th.T)
hf.create_dataset('phi',data=phi.T)
hf.create_dataset('length',data=length.T)
hf.create_dataset('vol',data=vol.T)
hf.create_dataset('gamma',data=gamma.T)
hf.create_dataset('rho',data=rho.T)
hf.create_dataset('intensity',data=intensity.T)
hf.create_dataset('nhat_r',data=nhhat_r.T)
hf.create_dataset('nhat_th',data=nhhat_th.T)
hf.create_dataset('nhat_phi',data=nhhat_phi.T)

# close #
hf.close()

######################################################
