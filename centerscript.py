#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 24 18:10:53 2023

@author: joseph
"""

# This script will take all of the dcd and pdb files in the directory it is currently in
# and output those same dcd files with the heparin chain in the center in chimerax. This will
# only work if there is a single heparin chain. Use this after using catanauto.py.

import mdtraj
import sys
import center_periodic
import getprefixandfilenames as gpaf 
from pathlib import Path


mydir = Path(__file__).parent
pdblist = gpaf.getpdblist()
replist = gpaf.replicatesstringslist()


for rep in replist: # Iterate over replicates to center them
    repdcdlist = gpaf.getdcdlistbyrep(rep) # For each replicate get the list of dcd files
    for (dcd,pdb) in zip(repdcdlist,pdblist): # Iterate over the list of dcd files and their corresponding pdb files
        # Center the dcd
        traj = mdtraj.load(dcd, top=pdb)
        center_periodic.center(traj, center_periodic.SelectionCriterion([-1]), query="name CA or name XX")
        traj.xyz += traj.unitcell_lengths[:, None] / 2
        
        # Make new name for centered dcd
        dcdpath = mydir / dcd
        dcdpathstr = str(dcdpath).replace(f".cat.map.{rep}.dcd", f".cat.map.centered.{rep}.dcd") 
        traj.save_dcd(dcdpathstr)





