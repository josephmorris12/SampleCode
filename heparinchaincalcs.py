#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 20 18:49:36 2023

@author: joseph
"""

# This script will hold various functions that do calculations for a single chain
# of heparin. The calculations will always give the actual output of the normal mdtraj function they
# utilize. Some of the functions in here are shared with calcssingleevaluaeachframe.py.

import numpy as np
import mdtraj as md

def persistancelength(traj):
    """
    This function will calculate the peristance length for every frame of the 
    simulation.
    
    Parameters
    ----------
    traj : TYPE Traj
        DESCRIPTION. Mdtraj of the simulation.

    Returns
    -------
    1D array with the radius of gyration for every frame.
    
    """
    ### Calc end to end distance squared in every frame (h^2)
    # Put start and end index into matrix to calculate end to end distance
    natoms = traj.xyz.shape[1]
    endtoendind = np.array([0, natoms-1], int)
    endtoendind = np.expand_dims(endtoendind, axis=0)

    # Calc end to end distance squared in every frame
    endtoenddist = md.compute_distances(traj, endtoendind)
    endtoenddistsquared = np.square(endtoenddist)

    ### Calc average length between all bond pairs in every frame (l)
    # Make matrix of atoms pairs to calc distance inbetween
    atomlist = list(range(1,natoms-1))
    atomarray = np.array(atomlist)
    atomarray = np.repeat(atomarray, 2)
    atomarray = np.append(0, atomarray)
    atomarray = np.append(atomarray, natoms-1)
    atomarray = atomarray.reshape((natoms-1, 2))

    # Calc distances between each pair of atoms
    distarray = md.compute_distances(traj, atomarray)

    # Calc average bond distance between each pair of atoms over each frame
    avgbonddistarray = distarray.mean(axis=1)

    ### Calc number of links (n)
    nlinks = natoms - 1

    ### Calc persistance length = (h^2)/(nl2) + l/2 in every frame
    nl2 = nlinks*avgbonddistarray*2.0
    ldivided2 = avgbonddistarray * 0.5
    persislen = np.divide(endtoenddistsquared.T, nl2) #+ (avgbonddistarray * 0.5)
    persislen = np.add(persislen, ldivided2).T
    persislen = persislen.flatten()
    return persislen


def radiusofgyration(traj):
    """
    This function calcs the radius of gyration for every frame for a heparin
    molecule.

    Parameters
    ----------
    traj : TYPE Traj
        DESCRIPTION. Mdtraj of the simulation.

    Returns
    -------
    1D array with the radius of gyration for every frame.
    
    """
    # Make mass array
    numheparin = int(traj.xyz.shape[1] / 2)
    massarray = np.array([310.1967, 248.12505], float)
    massarray = np.tile(massarray, int(numheparin))
    
    # Compute radius of gyration for each frame
    rgarray = md.compute_rg(traj, masses=massarray)
    rgarray = rgarray.flatten()
    return rgarray


def endtoenddistance(traj):
    """
    This function will plot calc the end to end distance in each frame.

    Parameters
    ----------
    traj : TYPE Traj
        DESCRIPTION. Mdtraj of the simulation.

    Returns
    -------
    None.
    1D array with all the end to end distance for every frame.

    """
    # Put start and end index into matrix to calculate end to end distance
    natoms = traj.xyz.shape[1]
    endtoendind = np.array([0, natoms-1], int)
    endtoendind = np.expand_dims(endtoendind, axis=0)

    # Calc end to end distance in every frame
    distarray = md.compute_distances(traj, endtoendind)
    distarray = distarray.flatten()
    return distarray


def dihedrals(traj):
    """
    This function calcs the bond angles of all the bonds in each frame for a trajectory.
    The dihedral angles are computed in quadruplets:
    (0,1,2,3), (1,2,3,4), (2,3,4,5), ...

    Parameters
    ----------
    traj : TYPE traj
        DESCRIPTION. TYPE Mdtraj of the simulation.

    Returns
    -------
    None.
    2D array with all the dihedral angles of all the dihedral quadruplets for every frame.
    """
    # Make quadruplets matrix to calculate dihedrals
    natoms = traj.xyz.shape[1]
    col1 = np.arange(0, natoms-3)
    col2 = np.arange(1, natoms-2)
    col3 = np.arange(2, natoms-1)
    col4 = np.arange(3, natoms)
    quadmatrix = np.c_[col1,col2,col3,col4]

    # Calc dihedrals over the quartets of atoms
    dihedirals = md.compute_dihedrals(traj, quadmatrix)
    return dihedirals


def bondangle(traj):
    """
    This function calcs the bond angles of all the bonds in each frame for a trajectory.
    The distances between atoms are calculated in triplets:
    (0,1,2), (1,2,3), (2,3,4), ...
    
    Parameters
    ----------
    traj : TYPE traj
        DESCRIPTION. Mdtraj of the simulation.

    Returns
    -------
    None.
    2D array with all the bond angles of all the angle triplets for every frame.
    
    """
    # Make matrix that holds atoms of simulation in triplets
    natoms = traj.xyz.shape[1]
    col1 = np.arange(0, natoms-2)
    col2 = np.arange(1, natoms-1)
    col3 = np.arange(2, natoms)
    tripletsmatrix = np.c_[col1,col2,col3]

    # Calc angles over the triples of atoms
    anglesarray = md.compute_angles(traj, tripletsmatrix)
    return anglesarray


def bondlengths(traj):
    """
    This function will calc the distances of all the chemical bond lengths in 
    each frame for a trajectory. The distances between atoms are calculated pairwise:
    (0,1), (1,2), (2,3), ...
    
    Parameters
    ----------
    traj : TYPE traj
        DESCRIPTION. Mdtraj of the simulation.

    Returns
    -------
    None.
    2D array with all the distances of all the bond pairs for every frame.
    
    """
    # Make matrix of atoms pairs to calc distances inbetween them
    natoms = traj.xyz.shape[1]
    atomlist = list(np.arange(1,natoms-1))
    atomarray = np.array(atomlist)
    atomarray = np.repeat(atomarray, 2)
    atomarray = np.append(0, atomarray)
    atomarray = np.append(atomarray, natoms-1)
    atomarray = atomarray.reshape((natoms-1, 2))

    # Calc distances between each pair of atoms
    distarray = md.compute_distances(traj, atomarray)
    return distarray
    
