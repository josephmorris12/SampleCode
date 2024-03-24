#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  9 01:57:50 2023

@author: joseph
"""

# This script will make two plots. One with the deccorelated averages of the number of
# parallel and antiparallel phf6 beta bridges

import matplotlib.pyplot as plt
import mdtraj as md
import scipy.stats
import pymbar
import numpy as np
import gzip
import pickle
import getprefixandfilenames as gpaf
import statistics

# Set up plots
fig, axes = plt.subplots(nrows=1 , ncols=1 , figsize=(6,6), squeeze=False)

# Initialize list of deccorelated means and standard errors for the beta bridges
ymeanlistfractionparallel = []
ystdlistparallel = []

# Get list of dat file to iterate over
# Change the replicate to do different replicate simulations
replicatelist = ["r0","r1","r2"]
for replicate in replicatelist:
    datlist = gpaf.getdatlist(replicate)
    
    # Load the dat files one at a time and find the standard error and mean of each one
    for datname in datlist: 
        with gzip.open(datname) as data_file:
            data = pickle.load(data_file)
            
            # Initialize the lists and variables that will hold the number of parallel and antiparallel beta bridges in each frame
            i_bb = data["i_bb"]
            n_bb = data["n_bb"]
            bbs = data["bbs"]
    
            n_parallel = []
            n_antiparallel = []
            
            # Iterate over all of the frames of one of the dat files to find the number of parallel and anti parallel beta bridges
            for i_bb_fr, n_bb_fr in zip(i_bb, n_bb):
                bbs_fr = bbs[i_bb_fr:i_bb_fr + n_bb_fr]
                lengths = bbs_fr["i_tail"] - bbs_fr["i_head"] + 1
    
                n_parallel.append(np.sum(lengths[bbs_fr["kind"] == 1]))
                n_antiparallel.append(np.sum(lengths[bbs_fr["kind"] == 2]))
                
                # Find the fraction of parallel beta bridges for each frame = parallel/(parallel+antiparallel)
                fractionparallel = []
                for (par,anti) in zip(n_parallel,n_antiparallel):
                    if (par+anti) == 0:
                        fractionparallel.append(0.0)
                    else:
                        fractionparallel.append(par/(par+anti))
                        
            # Find the mean and standard deviation and append it to lists to be graphed
            ymeanlistfractionparallel.append(statistics.mean(fractionparallel))
            ystdlistparallel.append(statistics.stdev(fractionparallel))
           
            
# Plot the data in an error bar plot
xlabelslist = replicatelist # Use floats to get it to fit in with the labels given by xlim and xticks
#axes[0, 0].set_xlim(left=-0.9, right=-3.6)
axes[0, 0].errorbar(xlabelslist, ymeanlistfractionparallel, ystdlistparallel, linestyle="dashed", capsize=3, fmt="o", mfc="w", color="black" )

# Label the axes for the parallel plot
axes[0, 0].set_ylabel("Fraction of Parallel Beta Bridges")
axes[0, 0].set_xlabel("Replicate Simulation")
axes[0, 0].set_title("a)", loc="left") 

# Bottom text
bottomtxt = "a) Refers to the average fraction of parallel beta bridges during the simulation. "
bottomtxt += "The average fraction was calculated by average fraction = (mean of sum of all parallel by frame)/(mean of sum of all parallel by frame + mean of sum of all antiparallel by frame). "
bottomtxt += "The error bars are standard deviations. "
plt.figtext(x=0.01, y=-0.025, s=bottomtxt, wrap=True, fontsize=7) 

# Save file
plt.savefig("fractionparallelbetabridges", dpi = 900, bbox_inches="tight")




