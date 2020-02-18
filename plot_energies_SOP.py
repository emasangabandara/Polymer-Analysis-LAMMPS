#!/usr/bin/env python2.7
import numpy as np
from matplotlib import pyplot as plt

def smoothHist(data) :
    y,binEdges=np.histogram(data,normed=True,bins=50)
    bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
    
    return bincenters , y
    
    
data = np.loadtxt('logs/data.dat')
Step, Temp, E_pair, E_vdwl, E_coul, E_bond, PotEng, KinEng = np.hsplit(data,8)
#Step, Temp, Press, E_pair, E_vdwl, E_coul, E_bond, PotEng, KinEng = np.hsplit(data,9)

#Plot time series
plt.figure(0)
plt.plot(Step,Temp,color='k')
plt.axhline(y=np.mean(Temp), color='r')
plt.xlabel('steps')
plt.ylabel('Temperature (K)')
plt.tight_layout()
plt.savefig('temp.ts.png')

plt.figure(1)
plt.plot(Step,E_pair,color='k',label=r'$E_{pair}$')
plt.plot(Step,E_vdwl,color='r',label=r'$E_{vdwl}$')
plt.plot(Step,E_coul,color='b',label=r'$E_{coul}$')
plt.plot(Step,E_bond,color='g',label=r'$E_{bond}$')
plt.plot(Step,PotEng,color='c',label=r'$Pot-Eng$')
plt.plot(Step,KinEng,color='m',label=r'$Kin-Eng$')
plt.xlabel('steps')
plt.ylabel('Energies (Kcal/mol')
plt.tight_layout()
plt.legend(loc='best')
plt.savefig('eng.ts.png')

#Plot distributions series
plt.figure(2)
plt.plot(smoothHist(Temp)[0],smoothHist(Temp)[1],'-',color='black', alpha=1.0, linewidth=1.5,label=r"$Temperature$")
plt.xlabel(r"$Temperature$ ($K$)")
plt.ylabel('Probability Density')
plt.legend(loc='best')
plt.savefig("temp-dist.png")

plt.figure(3)
Enb = E_pair + E_vdwl + E_coul
plt.plot(smoothHist(Enb)[0],smoothHist(Enb)[1],'-',color='black', alpha=1.0, linewidth=1.5,label=r"$E_{non-bonded} = E_{pair} + E_{vdwl} + E_{coul}$")
plt.plot(smoothHist(E_bond)[0],smoothHist(E_bond)[1],'-',color='red', alpha=1.0, linewidth=1.5,label=r"$E_{bond}$")
plt.xlabel(r"$Energies$ ($Kcal/mol$)")
plt.ylabel('Probability Density')
plt.legend(loc='best')
plt.savefig("Eenergy-dist.png")