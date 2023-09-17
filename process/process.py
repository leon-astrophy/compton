#####################################################
#
# Pre-process GRMHD snapshots
#
######################################################

#import#
import sys
import h5py
import pyharm
import numpy as np

######################################################

#read#
files = sys.argv[1]
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

# volume #
vol = gdet*dx1*dx2*dx3

######################################################

#rho[np.where(gamma < 1.1)] = 0

# repeat #
r = np.repeat(r, n3, axis=2)
th = np.repeat(th, n3, axis=2)
vol = np.repeat(vol, n3, axis=2)

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