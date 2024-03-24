#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 18:31:01 2023

@author: joseph
"""

import matplotlib.pyplot as plt
import calcssingleevalutaeachframehep as hp
import numpy as np
import scipy.stats
import pymbar
import mdtraj as md
import math 
import getprefixandfilenames as gpaf


def heparinerrorplot(pp, fu, he, xal, yal, conv):
    """
    This function will go through all of the files in listfilenames and make an
    error plot at the desired plot position. Use for the plots that use distance.

    Parameters
    ----------
    prefixmainfilenames : TYPE List of strings
        DESCRIPTION. List of all the prefix+main filenames of interest in the 
        analysis directory.
    fu : TYPE String
        DESCRIPTION. Aka function. Name of a function that will be used in calculating the 
        property of interest for the plot.
    pp : TYPE Tuple of ints
        DESCRIPTION. Aka plot position. Contains the row and column position for 
        the plot of interest.
    he : TYPE String
        DESCRIPTION. Aka header. The particular header you want for this plot.
    xal : TYPE String
        DESCRIPTION. Aka x-axis label. The particular label for the x-axis that
        you want for this plot.
    yal : TYPE String
        DESCRIPTION. Aka y-axis label. The particular label for the y-axis that
        you want for this plot.

    Returns
    -------
    None.
    Does a plot at the given plotposition using the given function. It will also
    put the header, do x-axis labeling, put the y-axis title, and put the x-axis 
    title as well.

    """
    ### Unpack plotposition
    pprow = pp[0]
    ppcol = pp[1]
    
    ### Make x-axis labels and ticks
    xaxislabels = gpaf.getprefixliststrings(2)
    xaxislabelsfloat = gpaf.getprefixlistfloats(2)
    axes[pprow, ppcol].set_xlim(left=-0.8, right=-4.1) # !!! Change this to change the x-axis start and end points
    axes[pprow, ppcol].set_xticks(ticks=xaxislabelsfloat, labels=xaxislabels, rotation=270, fontsize=4)
    
    # Get pdb and dcd filenames
    pdbfilenames = gpaf.getpdblist()
    dcdfilenames = gpaf.getdcdlist()
    
    ### Iterate over all files, do calculations with the specified function, 
    ### and put averages and error on the plot of interest to make the error plot.
 
    ymeans = np.array([], float)
    yerr = np.array([], float)
    for (pdb,dcd) in zip(pdbfilenames,dcdfilenames):
        traj = md.load(dcd, top=pdb)
        fulistcalcs = fu(traj)
        
        # Check for the proper conversion factor
        if conv == "A":
            fulistcalcsamp = fulistcalcs * 10.0 # convert nm to A
        else:
            fulistcalcsamp = fulistcalcs * (180/math.pi) # Convert radians to degrees
        
        
        inddec = pymbar.timeseries.subsample_correlated_data(fulistcalcsamp) # Find indices to uncorrelate data #!!!
        decorfulistcalcs = fulistcalcsamp[inddec]
        meandecor = decorfulistcalcs.mean(axis=0) # Calculate mean of uncorrelated data
        stedecor = scipy.stats.sem(decorfulistcalcs) # Calculate standard error of mean of uncorrelated data
        ymeans = np.hstack((ymeans, meandecor))
        yerr = np.hstack((yerr, stedecor))
    axes[pprow, ppcol].errorbar(xaxislabelsfloat, ymeans, yerr, linestyle="dashed", capsize=5, fmt="o", mfc="w", color="black" )
    
    ### Place titles and headers for the plot
    # Labeling titles
    axes[pprow, ppcol].set_xlabel(xal, fontsize=6)
    axes[pprow, ppcol].set_ylabel(yal, fontsize=8)
    axes[pprow, ppcol].set_title(he, loc="left")



if __name__ == '__main__':
    # Set up subplots
    fig, axes = plt.subplots(nrows=3 , ncols=2, figsize=(10.2,13.2), squeeze=False) #!!!
    # plt.subplots_adjust(left=, right=, bottom=, top=, wspace=, hspace=)
    
    ###### Make sure everything is in the proper order!!!!!!!
    
    ### Set up lists to iterate over for graphing
    # Set up positions for plots
    plotposlistoftuples = [(0,0), (1,0), (2,0), (0,1), (1,1), (2,1)] #!!!
    # Set up list of functions we want to use
    funclist = [hp.bondlengths, hp.bondangle, hp.dihedrals, hp.persistancelength, hp.radiusofgyration, hp.endtoenddistance]
    # Set up the list of headers for the individual graphs
    headerlist = ["a)", "b)", "c)", "d)", "e)", "f)"] #!!!
    # Set up the list of x-axis labels
    xaxislabels = ["A1 of Dihedral_H_H_H_H"] * 6 #!!!
    # Set up the list of y-axis labels
    yaxislabels = ["Bond Length (A)",
                   "Bond Angle (Degrees)",
                   "Dihedral Angle (Degrees)",
                   "Persistance Length (A)",
                   "Radius of Gyration (A)",
                   "End to End Distance (A)"] #!!!
    conversionlist= ["A", "Degrees", "Degrees", "A", "A", "A", "A"]

    
    ### Use for loop to make graph. 
    # Each loop should plot the data from all of the given files onto one of the plots using one of the respective functions.
    for (pp,fu,he,xal,yal,conversion) in zip(plotposlistoftuples, funclist, headerlist, xaxislabels, yaxislabels, conversionlist):
        heparinerrorplot(pp,fu,he,xal,yal,conversion)
    
    
    # Bottom text
    bottomtxt = "Command: ./run.py title 500 350:350 “30”. "
    bottomtxt += "The x-axis labels refer to A1 of dihedral_H_H_H_H which was changed in order to make trans more favorable by making the delta E between cis and trans greater. "
    bottomtxt += "a) Average bond lengths for the entire simulation along with standard error bars of one heparin chain is graphed. The average is found for each frame for every trajectory. The decorrelation and the average for the overall simulation for that trajectory is based on these averages. " 
    bottomtxt += "b) Average bond angles for the entire simulation along with standard error bars of one heparin chain is graphed. The average is found for each frame for every trajectory. The decorrelation and the average for the overall simulation for that trajectory is based on these averages. "
    bottomtxt += "c) Average dihedrals for the entire simulationalong with standard error bars of one heparin chain is graphed. The average is found for each frame for every trajectory. The decorrelation and the average for the overall simulation for that trajectory is based on these averages. "
    bottomtxt += "d) Average persistance length for the entire simulation along with standard error bars of one heparin chain is graphed. " 
    bottomtxt += "e) Average radius of gyration for the entire simulation along with standard error bars of one heparin chain is graphed. "
    bottomtxt += "f) Average end to end distance for the entire simulation along with standard error bars of one heparin chain is graphed. "
    bottomtxt += "Uncorrelated data is used for the mean and standard error calculations." #!!!
    plt.figtext(x=0.01, y=0.03, s=bottomtxt, wrap=True, fontsize=6.5) #!!!
    
    # Output
    plt.savefig("sixplots", dpi = 900, bbox_inches="tight")


