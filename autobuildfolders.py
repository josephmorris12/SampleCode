#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 26 15:18:23 2023

@author: joseph
"""

# This script should be placed inside the autosetupresearchfiles folder with all of 
# the other files it will copy and modify. This script will make a new main directory and place
# all of the supplied trajectory files inside of their own respective trajectory folders within that 
# directory. This will make it easy to set up the research files for pod. Copy this script
# to a new file every time you want to use it and make edits to it.

##### Always make sure to keep the proper orders of names, etc.

import os
import shutil
import changestring as cs
from pathlib import Path
import numpy as np
import subprocess



###############################################################################
# This is not linked to anything in the script. Use it to find inputs for cgnewstrings and runnewstrings
# The input of cgnewstrings and runnewstrings should be edited to be a list of strings
# with the exact format you want in the file. If needed just rewrite this part of the code,
# and use lists and lists of lists as needed to modify fileprefixlist, cg/runoriginalstrings, and cg/runnewstrings.

### 1. Find strings to change to by using multipliers
originalfloats = [5.2384242763002895] #!!! List of original floats you want changed
multipliers = [0.9,1.1,1.3,1.5,1.7,1.8] #!!! List of multipliers for those floats
cs.multipliers(originalfloats, multipliers)
raise Exception("Only doing the calculations with multipliers.") 

### 2. Find strings to change to by using a hard coded list of floats directly
originalfloat = [2.5] #!!! List of single original float you want changed
newfloatslist = np.linspace(start=80, stop=170, step=5) #!!! List of new floats you want the original floats to be changed to
cs.changeonefloateachfilelists(originalfloat, newfloatslist)
raise Exception("Converting list of floats into list of lists.") 

# This will give the list of lists that you iterate over while using cg.editfilenewstrings.
# Each list inside the list of lists corresponds to the new values you want for a 
# particular md file. It will also give you a list of multipliers as strings.
###############################################################################



# Make names for prefix_name
# gives the prefixes of the titles in order
filetitleonly = "mytest" #!!! Title for new directory and main part of title for new files
filetitle = "~" + filetitleonly # Add '~' so it is easy to split file titles later
fileprefixlist = ["0.9", "1.1", "1.3", "1.5", "1.7", "1.8"] #!!!
titlelist = [(i + filetitle) for i in fileprefixlist] # List of folder and file titles

# Path to analysispythonscripts which holds this file and the others needed inside
mymaindir = str(Path(__file__).parent)
# Path to folder where want new working directory to put files in, the hardcoded path
# input here is just the location where you want the new main directory
dirwherewantnewsimfiles = "/home/joseph/chaininvestigation" #!!!
myoutputpath = dirwherewantnewsimfiles + "/" + filetitleonly 

# Make new directory to put folders in
os.mkdir(myoutputpath)
os.chdir(myoutputpath)
# Make new folders in the new directory
newfiledir = [] # List of absolute paths of folders in new directory
for var in titlelist:
    newdir = os.path.join(myoutputpath, var)
    newfiledir.append(newdir)
for dirname in newfiledir:
    os.mkdir(dirname)
os.mkdir("analysis" + filetitle)
os.chdir(mymaindir)


### 1. Put cg.profile into each new directory
for x in newfiledir:
    subprocess.run(["cp", "cg.profile", x]) # Copy cg.profile which does not have write permissions (will make deleting folders hard)


### 2. Put run.py into each new directory
# MAKE SURE USING CORRECT run.py!!!
# Put runoriginalstrings as [] or runnewstrings as [[]] when just want to copy "run.py" with no changes
runname = "nophf61heparinrun.py" #!!!
runoriginalstrings = [] #!!! List of original strings you want to change
runnewstrings = [[]] #!!! List of lists; each "row" corresponds to strings that you want to change originalstrings to

# Case where want to copy with no changes
if (runoriginalstrings == []) or (runnewstrings == [[]]):
    print(f"No strings were set. The {runname} files will be output without modifications.")
    for y in newfiledir:
        shutil.copy2(runname, y)
# Case where want to copy with changes
for pas,newstrs in zip(newfiledir,runnewstrings):
    cs.editfilenewstrings(runname, runoriginalstrings, newstrs, outputpath=pas)


### 3. Adjust and put default_slurm.sh into each new directory
slurmname = "default_slurm.sh" 
originalstrings = ["myfilename", "mycommand", "168"] 
maxlengthofsim = "168" #!!! The max length of the simulation in hours. Short max simulation times get nodes faster, but the simulation may cut off early if you make the max time to short
for ind,file in enumerate(newfiledir):
    newstrings = [titlelist[ind], f"./{runname} -- {titlelist[ind]} 500 350:350 '30'", maxlengthofsim] #!!! Put command here
    cs.editfilenewstrings(slurmname, originalstrings, newstrings, outputpath=file)


### 4. Put cg.ff into each directory
# Use when want to change cg.ff
# Put cgoriginalstrings as [] or cgnewstrings as [[]] when just want to copy "cg.ff" with no changes
cgname = "cg.ff" 
cgoriginalstrings = ['5.2384242763002895'] #!!! List of original strings you want to change
cgnewstrings = [['4.714581848670261'], 
                ['5.762266703930319'], 
                ['6.8099515591903765'], 
                ['7.857636414450434'], 
                ['8.905321269710491'], 
                ['9.429163697340522']] #!!! List of lists each "row" corresponds to strings that you want to change originalstrings to

# Case where want to copy with no changes
if (cgoriginalstrings == []) or (cgnewstrings == [[]]):
    print(f"No strings were set. The {cgname} files will be output without modifications.")
    for z in newfiledir:
        #shutil.copy2(cgname, z)
        subprocess.run(["cp", cgname, z]) # Copy cg.ff which does not have write permissions 
        subprocess.run(["chmod", "+w", z + f'/{cgname}']) # Give the copy of cg.ff write permissions so it is easier to delete folders
# Case where want to copy with changes (copy through editfilenewstrings() gives write permission automatically)
for pa,newstr in zip(newfiledir,cgnewstrings):
    cs.editfilenewstrings(cgname, cgoriginalstrings, newstr, outputpath=pa)


### 5. Fill the analysis folder with catanaauto.py and sbatchauto.py
analysispath = myoutputpath + "/analysis" + filetitle
shutil.copy2("catanaauto.py", analysispath)
shutil.copy2("sbatchauto.py", analysispath)
shutil.copy2("getprefixandfilenames.py", analysispath)
shutil.copy2("centerscript.py", analysispath)







