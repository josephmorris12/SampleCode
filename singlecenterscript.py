#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 10:54:30 2023

@author: joseph
"""

# This script makes a centered dcd file for a single heparin chain simulation.
# Input the dcd, pdb, and new dcd name through command line
# Generally use centerscript.py as it does this for all of the simulations in a directory.

import mdtraj
import sys
import center_periodic

dcd, pdb, out_dcd = sys.argv[1:]
traj = mdtraj.load(dcd, top=pdb)
center_periodic.center(traj, center_periodic.SelectionCriterion([-1]), query="name CA or name XX")
traj.xyz += traj.unitcell_lengths[:, None] / 2
traj.save_dcd(out_dcd)
