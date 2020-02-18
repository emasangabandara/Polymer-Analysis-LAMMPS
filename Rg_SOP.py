#!/usr/bin/env python2.7
import numpy as np
from MDAnalysis import *
import MDAnalysis
from matplotlib import pyplot as plt

def smoothHist(data) :
    y,binEdges=np.histogram(data,normed=True,bins=50)
    bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
    
    return bincenters , y

psf = 'Abeta42_SOP_sc_frags1.psf'
dcd = 'dump/prod1.dcd'

u = Universe(psf,dcd)

Rgyr = []

protein = u.select_atoms("name *")

for ts in u.trajectory:
    Rgyr.append((u.trajectory.time, protein.radius_of_gyration()))
	
Rgyr = np.array(Rgyr)
print 'mean-Rg', np.mean(Rgyr[:,1])

np.save('Rgyr',Rgyr)

plt.figure(0)
plt.plot(Rgyr[:,0], Rgyr[:,1], 'r--', lw=2, label=r"$<R_G>=$"+ str(np.mean(Rgyr[:,1])))
plt.xlabel("Time (steps)")
plt.ylabel(r"Radius of gyration $R_G$ ($\AA$)")
plt.tight_layout()
plt.legend(loc='best')
plt.savefig("Rg.png")

plt.figure(1)
plt.plot(smoothHist(Rgyr[:,1])[0],smoothHist(Rgyr[:,1])[1],'-',color='black', alpha=1.0, linewidth=1.5,label=r"$R_G-distribution$")
plt.xlabel(r"Radius of gyration $R_G$ ($\AA$)")
plt.ylabel('Probability Density')
plt.legend(loc='best')
plt.savefig("Rg-dist.png")
