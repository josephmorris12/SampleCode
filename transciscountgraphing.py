#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 12:35:43 2023

@author: joseph
"""

# This script plots the actual deccorelated average amounts of cis and trans counts there are during a simulation
# It does this using average errorbar plots. 

import getprefixandfilenames as gpaf
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import mdtraj as md
import scipy.stats
import pymbar
import heparinchaincalcs as hcc
import math
import numpy as np

# Set up subplots
fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(7, 7), squeeze=False)
# plt.subplots_adjust(left=#, right=#, bottom=#, top=#, wspace=#, hspace=#)

### Make x-axis labels and ticks
xaxislabelsfloat = gpaf.getprefixlistfloats(2)
xaxislabelsstr = gpaf.getprefixliststrings(2)
axes[0, 0].set_xlim(left=-0.8, right=-4.1) # !!! Change this to change the x-axis start and end points
axes[0, 0].set_xticks(ticks=xaxislabelsfloat, labels=xaxislabelsstr, rotation=0, fontsize=7)


### Make error plots using direct quantities of cis and trans
pdblist = gpaf.getpdblist()
dcdlist = gpaf.getdcdlist()
ymeanscis = []
yerrscis = []
ymeanstrans = []
yerrstrans = []

for dcd,pdb in zip(dcdlist, pdblist):
    # Initialize variables
    numciseachframemat = []
    numtranseachframemat = []
    
    # Load the trajectory
    traj = md.load(dcd, top=pdb)
    dihedralsmat = hcc.dihedrals(traj)
    dihedralsmatdegrees = dihedralsmat * (180/math.pi)
    
    # Count the number of cis and trans dihedrals in each frame
    for matrow in dihedralsmatdegrees: # Iterate over rows of matrix
        # Initialize counters
        numciscountereachframe = 0
        numtranscountereachframe = 0
        for di in matrow: # Iterate over columns of each row
            if di >= -45 and di <= 45:
                numciscountereachframe += 1
            if di >= 135 or di <= -135:
                numtranscountereachframe += 1
        # Add the counters to the matrix
        numciseachframemat.append(numciscountereachframe)
        numtranseachframemat.append(numtranscountereachframe)
    
    # Decorrelate the data for cis
    cisinddec = pymbar.timeseries.subsample_correlated_data(numciseachframemat) # Get indices
    cisdecorrelateddata = [numciseachframemat[i] for i in cisinddec]
    cisdecorrelateddatamat = np.array([cisdecorrelateddata], float)
    cismeandec = cisdecorrelateddatamat.mean() # Calculate mean of uncorrelated data
    cisstedec = scipy.stats.sem(cisdecorrelateddata) # Calculate standard error of mean of uncorrelated data
    
    # Append meandec and stedec to matrices for cis
    ymeanscis.append(cismeandec)
    yerrscis.append(cisstedec)
    
    # Decorrelate the data for trans
    transinddec = pymbar.timeseries.subsample_correlated_data(numtranseachframemat) # Get indices
    transdecorrelateddata = [numtranseachframemat[i] for i in transinddec]
    transdecorrelateddatamat = np.array([transdecorrelateddata], float)
    transmeandec = transdecorrelateddatamat.mean() # Calculate mean of uncorrelated data
    transstedec = scipy.stats.sem(transdecorrelateddata) # Calculate standard error of mean of uncorrelated data
    
    # Append meandec and stedec to matrices for trans
    ymeanstrans.append(transmeandec)
    yerrstrans.append(transstedec)

# Plot direct number plots
axes[0, 0].errorbar(xaxislabelsfloat, ymeanscis, yerrscis, linestyle="dashed", capsize=5, fmt="o", mfc="w", color="red")
axes[0, 0].errorbar(xaxislabelsfloat, ymeanstrans, yerrstrans, linestyle="dashed", capsize=5, fmt="o", mfc="w", color="blue")


### Label plots
# Make key
cislabel = mpatches.Patch(color="red", label="Cis")
translabel = mpatches.Patch(color="blue", label="Trans")
fig.legend(handles=[cislabel, translabel], )



# Labeling titles for single plot
xaxistitle = "A1 of Dihedral_H_H_H_H"
axes[0, 0].set_xlabel(xaxistitle)
yaxistitle = "Average Number of Cis or Trans Dihedrals"
axes[0, 0].set_ylabel(yaxistitle)
header = "Cis vs Trans Variation due A1 if Dihedral_H_H_H_H"
axes[0, 0].set_title(header, loc="center")
bottomtxt = "Command: ./run.py title 500 350:350 “30”. "
bottomtxt += "The x-axis labels refer to what A1 of Dihedral_H_H_H_H was modified to. "
bottomtxt += "Uncorrelated data is used for the mean and standard error calculations." #!!!
plt.figtext(x=0.03, y=0.03, s=bottomtxt, wrap=True, fontsize=7)


### Output
plt.savefig("cistrans", dpi = 900, bbox_inches="tight")

