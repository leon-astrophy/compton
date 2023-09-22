#####################################################
#
# Pre-process GRMHD snapshots to fed into xraypol
#
######################################################

#import#
import sys
import h5py
import pyharm
import numpy as np

######################################################

# terminal input #
files = sys.argv[1]
outfile = sys.argv[2]
radcut = float(sys.argv[3])
sigmacut = float(sys.argv[4])

######################################################
# read #
dump = pyharm.load_dump(files)

# assign #
n1 = dump['n1']
n2 = dump['n2']
n3 = dump['n3']

# assign #
dx1 = dump['dx1']
dx2 = dump['dx2']
dx3 = dump['dx3']

# load coordinate #
r = dump['r']
th = dump['th']
phi = dump['phi']

# load mass #
rho = dump['rho']
gdet = dump['gdet']
gamma = dump['Gamma']
sigma = dump['sigma']

######################################################
# sigma cut #
rho[np.where(sigma > sigmacut)] = 0

######################################################

# volume #
vol = gdet*dx1*dx2*dx3

# repeat #
r = np.repeat(r, n3, axis=2)
th = np.repeat(th, n3, axis=2)
vol = np.repeat(vol, n3, axis=2)

######################################################
# look for radius cut #
for i in range (0, r.shape[0]):
  if(r[i,0,0] > radcut):
    i_cut = i + 1
    break

######################################################
# reduce the size of the array #
r = r[0:i_cut,:,:]
th = th[0:i_cut,:,:]
phi = phi[0:i_cut,:,:]
vol = vol[0:i_cut,:,:]
rho = rho[0:i_cut,:,:]
gamma = gamma[0:i_cut,:,:]

######################################################

# dmass #
dmass = vol*rho
mass = np.sum(dmass.flatten())

# nebar #
nebar = rho/mass

######################################################

# output #
hf = h5py.File("pp-"+files, 'w')

# create data #
hf.create_dataset('n1',data=n1)
hf.create_dataset('n2',data=n2)
hf.create_dataset('n3',data=n3)
hf.create_dataset('r',data=r.T)
hf.create_dataset('th',data=th.T)
hf.create_dataset('phi',data=phi.T)
hf.create_dataset('vol',data=vol.T)
hf.create_dataset('gamma',data=gamma.T)
hf.create_dataset('nebar',data=nebar.T)

# close #
hf.close()

######################################################