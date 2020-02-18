#!/usr/bin/env python2.7
import numpy as np
from MDAnalysis import *
import MDAnalysis
from MDAnalysis.analysis import contacts
from MDAnalysis.analysis import distances
import matplotlib.colors as colors
from matplotlib import pyplot as plt

contact_cutoff = 6.0 #Angstrom

#SOP-SC topology
#ALA GLY CAs have the aminoacid name and no side chain
#All other amino acids have CA and amino acid name for the respective side chain
psf = 'Abeta40_SOP_sc_frags1.psf'
dcd = 'dump/prod1.dcd'

u = Universe(psf,dcd)

end_f   = u.trajectory.n_frames #Total number of frames in the trajectory
start_f = 0
skip    = 1


contact_mat_ens = []

u = Universe(psf,dcd)

for ts in u.trajectory[start_f:end_f:skip]:

    Calpha      = u.select_atoms("resname ALA GLY CA")
    Calpha_coor = Calpha.positions
    contact_mat = distances.contact_matrix(Calpha_coor, cutoff=contact_cutoff, returntype='numpy', box=None)
    contact_mat_ens.append(contact_mat)

mean_contact_mat = np.mean(np.asarray(contact_mat_ens),axis=0)
np.save('mean_contact_mat',mean_contact_mat)
plt.figure(0)
plt.pcolor(mean_contact_mat)
plt.colorbar()
plt.xlabel("Residue Index")
plt.ylabel("Residue Index")
plt.savefig("contact_map_cutoff_6AV3.png")

plt.figure(2)
mean_contact_mat[mean_contact_mat == 0.0] = 0.000000000000001 #Assign very small value to zeros to prevent log(0) error.
plt.pcolor(mean_contact_mat,norm=colors.LogNorm(vmin=mean_contact_mat.min(), vmax=mean_contact_mat.max()))
plt.colorbar()
plt.xlabel("Residue Index")
plt.ylabel("Residue Index")
plt.savefig("contact_map_cutoff_6A_logV3.png")