#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  4 21:11:17 2023

@author: joseph
"""

# The purpose of this file is to automatically use the concatenate_traj.py script to 
# concatenate files together to get the full concatenated dcd and pdb files for analysis. 
# Then the analyze_traj.py script will be used to convert these into the dat file 
# that can also be used for analysis. This should be placed inside the analysis folder. 
# All of the fully concatenated dcd, pdb, and dat files will be put into the graphing folder. 
# All md folders and their contents should be named prefix-main-filetype with folders having no filetype.

import os
import subprocess
from pathlib import Path
import argparse


parser = argparse.ArgumentParser(
    description = """The purpose of this file is to automatically use the concatenate_traj.py script to 
    concatenate files together to get the full concatenated dcd and pdb files for analysis. 
    Then the analyze_traj.py script will be used to convert these into the dat file 
    that can also be used for analysis. This should be placed inside the main directroy
    inside the graphing folder. All of the fully concatenated dcd, pdb, and 
    dat files will be put into the graphing folder. All md folders and their contents
    should be named prefix-main-filetype with folders having no filetype."""
)
parser.add_argument("--doanalysis", "-d", default = "True",
                    help = """True or False input. Decides if you want to use
                    analyze_traj.py on the concatenated files. The default is
                    True.""")
args = parser.parse_args()
boolanalysis = args.doanalysis


### Get mainname of the folders and files based on the analysis-mainname folder
analysispath = Path(__file__).parent
analysispathname = analysispath.name
main = analysispathname.split("~")[1]
main = "*" + main + "*"


### Get names of all of the trajectory folder paths
scrpathobj = Path(__file__)
mainobj = scrpathobj.parent.parent
mdtrajpathsobj = mainobj.glob(main)
mdnamelist = []
# Convert path objects to strings
for pt in mdtrajpathsobj:
    mdnamelist.append(str(pt.name))
# Get rid of the analysis filename
for ind,nm in enumerate(mdnamelist):
    if "analysis" in nm:
        mdnamelist.pop(ind)
mdnamelist.sort()


### Get paths of all of the trajectory folder paths
graphingdir = os.getcwd()
mdpathslist = []
mdtrajpathsobj = mainobj.glob(main)
for pat in mdtrajpathsobj:
    mdpathslist.append(str(pat))
# Get rid of the analysis path
for inx,nu in enumerate(mdpathslist):
    if "analysis" in nu:
        mdpathslist.pop(inx)
mdpathslist.sort()


### Make list of replicates suffixes to iterate over
# Find the number of replicates
firstmdpath = Path(mdpathslist[0]) # Make path object from the first path (the same number of replicates will be in each folder)
firstdcdfirstfolderfilename = mdnamelist[0] + "_000.r?.dcd" # Use the first dat file to find the number of replicates
firstpathglob = firstmdpath.glob(firstdcdfirstfolderfilename) # Make list of specific file that has r#
numreplicates = len(list(firstpathglob)) # Take list of specific file that has r# so now how many replicates there are

# Make replicates suffix list to interate over
intlist = list(range(numreplicates))
repsuffixlist = ["r" + str(myint) for myint in intlist]


### Use concatenate_traj.py on all of the md folders and mv the important dcd, pdb, and dat outputs to the graphing folder
for pa, mdname in zip(mdpathslist, mdnamelist): # Iterate over the folder names and their paths
    # Change directory to a md folder
    os.chdir(pa)
    subprocess.run(["concatenate_traj.py", f"{mdname}"])
    if boolanalysis != "True": # Only use analysis_traj.py if you want to   
        subprocess.run(["mv", "--", f"{mdname}.cat.map.top.pdb", graphingdir])
        for repsuf in repsuffixlist: # Iterate over replicate suffix list to move all of them to analysis folder
            subprocess.run(["mv", "--", f"{mdname}.cat.map.{repsuf}.dcd", graphingdir])  
    else:
        for repsuf in repsuffixlist: # Make dat files by iterate over replicate suffixes
            subprocess.run(["analyze_traj.py", f"{mdname}.cat.map.{repsuf}.dcd", f"{mdname}.cat.map.top.pdb", f"{mdname}.{repsuf}.dat"])
        for repsuf in repsuffixlist: # Move the files that have replicates suffixes to analysis directory
            subprocess.run(["mv", "--", f"{mdname}.cat.map.{repsuf}.dcd", f"{mdname}.{repsuf}.dat", graphingdir])
        subprocess.run(["mv", "--", f"{mdname}.cat.map.top.pdb", graphingdir]) # Move pdb to analysis directory (does not have a replicate suffix)

   
