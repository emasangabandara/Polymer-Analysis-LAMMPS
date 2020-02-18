#!/usr/bin/env python2.7
import numpy as np
from MDAnalysis import *
import MDAnalysis
from MDAnalysis.analysis.rdf import InterRDF
from matplotlib import pyplot as plt

psf = 'Abeta42_SOP_sc_frags1.psf'
dcd = 'dump/prod1.dcd'

u = Universe(psf,dcd)

protein = u.select_atoms("resname *")

rdf = InterRDF(protein, protein, nbins=100)
rdf.run()

plt.plot(rdf.bins[1:], rdf.rdf[1:])
plt.ylabel(r"$g(r)$")
plt.xlabel(r"$r$ ($\AA$)")
plt.savefig("RDF.png")
