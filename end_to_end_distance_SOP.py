#!/usr/bin/env python2.7
import numpy as np
from MDAnalysis import *
import MDAnalysis
from scipy.spatial import distance
from matplotlib import pyplot as plt
#@TODO Fix label texts
def smoothHist(data) :
    y,binEdges=np.histogram(data,normed=True,bins=50)
    bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
    
    return bincenters , y

psf = 'Abeta40_SOP_sc_frags1.psf'
dcd = 'dump/prod1.dcd'

u = Universe(psf,dcd)

Retoe = []

protein = u.select_atoms("name *")

for ts in u.trajectory:
    res_first = protein[0].position
    res_last  = protein[-1].position    
    etoe = distance.euclidean(res_first,res_last)
    print etoe 
    Retoe.append((u.trajectory.time, etoe))
	
Retoe = np.array(Retoe)
print 'mean-etoe', np.mean(Retoe[:,1])

np.save('etoe',Retoe)

plt.figure(0)
plt.plot(Retoe[:,0], Retoe[:,1], 'r--', lw=2, label=r"$<R_{end-to-end}>=$"+ str(np.mean(Retoe[:,1])))
plt.xlabel("Time (steps)")
plt.ylabel(r"$R_{end-to-end}$ ($\AA$)")
plt.tight_layout()
plt.legend(loc='best')
plt.savefig("Retoe.png")

plt.figure(1)
plt.plot(smoothHist(Retoe[:,1])[0],smoothHist(Retoe[:,1])[1],'-',color='black', alpha=1.0, linewidth=1.5,label=r"$R_{end-to-end}$")
plt.xlabel(r"$R_{end-to-end}$ ($\AA$)")
plt.ylabel('Probability Density')
plt.legend(loc='best')
plt.savefig("etoe-dist.png")