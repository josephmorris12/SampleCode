#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  8 22:13:24 2023

@author: joseph
"""

# This script will use the dat files in the analysis directory 
# to make plots for the replicates of a simulation. The simulation will be decided 
# by the prefixes in the analysis directory. The plots will plot frames against the number of phf6 
# in the largest cluster in a line graph. This should always work as it uses 
# the custom analyze_traj.py, so you don't need to mess with slicing trajectories
# to get rid of chains.


import gzip
import matplotlib.pyplot as plt
import pickle
import getprefixandfilenames as gpaf
import numpy as np
from pathlib import Path

## Initialize variables
pos = (0,0)
mypath = Path(__file__).parent
prefixlist = gpaf.getprefixliststrings(2) # Make sure proper number of decimals here for glob to work!!!!!!!!!
replist = gpaf.replicatesstringslist() # Get corresponding replicate strings for dat files


## Plot
for prefix in prefixlist: # Iterate over prefixes to do plots
    fig, axes = plt.subplots(nrows=1 , ncols=1, squeeze=False)
    datlist = sorted(mypath.glob(f"{prefix}*.dat")) # Check for dat files based on prefix and sort based on replicate number
    for (dat,rep) in zip(datlist,replist): ## Iterate over dat files and replicate numbers to do plot for one prefix at a time
        with gzip.open(dat) as data_file: # Open each dat file with its respective replicate string and do the plot
            data = pickle.load(data_file)
            clustersizesarray = data["sc_clust_sizes"]
            largestclusterlist = []
            for framerow in clustersizesarray:
                largestclusterlist.append(np.max(framerow))
            
            numberframes = len(largestclusterlist)
            xaxisframeslist = np.arange(1,numberframes+1)
            axes[pos].plot(xaxisframeslist, largestclusterlist, label=rep)
    plt.legend(bbox_to_anchor=(1,1.1))
    
    # Make ticks and labels for x-axis, set limits for x-axis, plot scales for x-axis
    startxticks = 0
    stopxticks = numberframes
    numofxticks = 8
    roundto = 0
    ticksx = list(np.linspace(start=startxticks, stop=stopxticks, num=numofxticks))
    labelsx = [str(round(xti, roundto)) for xti in ticksx]
    # Set x-axis labels
    axes[pos].set_xlim(left=startxticks, right=stopxticks)
    axes[pos].set_xticks(ticks=ticksx, labels=labelsx) # ticks is location, and labels is actual string labels
    
    bottomtxt = "Aggregation of phf6 over time."
    plt.figtext(x=0.4, y=0, s=bottomtxt, wrap=True, fontsize=7) 
    axes[pos].set_xlabel("Frame")
    axes[pos].set_ylabel("Number of phf6 in the largest cluster")
    
    plt.savefig(f"{prefix}phf6aggregation.png", dpi = 900, bbox_inches="tight")
