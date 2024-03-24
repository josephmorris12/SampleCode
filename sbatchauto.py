#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 25 01:18:33 2023

@author: joseph
"""

#This python script should be placed in the autosetupresearchfiles 
#directory, and it should be copied into every prefix~analysis folder. This 
#script will go through all of the trajectory folders in the main directory, 
#and it will sbatch default_slurm for that trajectory.

import os
import subprocess
from pathlib import Path

# Get mainname of the folders and files based on the analysis~mainname folder
analysispath = Path(__file__).parent
analysispathname = analysispath.name
mainname = analysispathname.split("~")[1]

# Make path objects for all of the trajectory folders in the main directory
mainpath = Path(__file__).parent.parent
mainpathstr = str(mainpath)
mainpathglob = sorted(mainpath.glob(f"*{mainname}"))
noanalysismainpath = []
for mp in mainpathglob:
    if ("analysis" in mp.name) != True:
        noanalysismainpath.append(mp)

# Iterate over noanalysismainpath to sbatch default_slurm in each trajectory folder
for trajfolder in noanalysismainpath:
    trajfolderpathstr = str(trajfolder.absolute())
    os.chdir(trajfolderpathstr)
    subprocess.run(["sbatch", "default_slurm.sh"])
    





