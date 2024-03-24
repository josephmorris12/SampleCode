#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 08:58:29 2023

@author: joseph
"""

# This module should be placed in autosetupresearchfiles, and it should be placed
# in the analysis~mainname folder of every simulation directory. This file will contain
# functions that find the mainname, prefixes, pdb names, dcd names, and dat names
# of the files within this modules working directory.

# All files/folders are of the form prefix~mainname.suffix
# Suffix will always have multiple "."

# This script can also deal with multiple replicates

from pathlib import Path


def getmainname():
    """
    This function will find the mainname of the folders and files.

    Returns
    -------
    A string of the mainname of the folders and files.

    """
    # Get mainname of the folders and files based on the analysis~mainname folder
    analysispath = Path(__file__).parent
    analysispathname = analysispath.name
    mainname = analysispathname.split("~")[1]
    return mainname


def getpdblist():
    """
    This function will find all of the pdb filenames in your working directory.

    Returns
    -------
    A sorted list of strings with all of the pdb filenames from your working directory.
    They are sorted from least to greatest.

    """
    workingpath = Path(__file__).parent
    allpdbpaths = sorted(workingpath.glob("*.pdb"))
    allstrpdbs = [f.name for f in allpdbpaths]
    return allstrpdbs


def getdcdlistbyrep(replicate="r0"):
    """
    This function will find all of the dcd filenames in your working directory 
    for a specified replicate (accept those that are "centered").
    
    Parameters
    ----------
    numdecimals : TYPE str
        DESCRIPTION. Replicate suffix of the replicate you are interested in. The default 
        is r0.

    Returns
    -------
    A sorted list of strings with all of the dcd filenames for the replicate of 
    interest from your working directory. They are sorted from least to greatest.

    """
    workingpath = Path(__file__).parent
    alldcdpaths = sorted(workingpath.glob(f"*map.{replicate}.dcd"))
    allstrdcd = [f.name for f in alldcdpaths]
    return allstrdcd


def getdcdlistbyprefix(prefix):
    """
    This function will find all of the dcd filenames in your working directory 
    for a specified mainname (accept those that are "centered").
    
    Parameters
    ----------
    numdecimals : TYPE str,int,float
        DESCRIPTION. Prefix of the simulation you are interested in.

    Returns
    -------
    A sorted list of strings with all of the dcd filenames for the simulation of 
    interest from your working directory. They are sorted from least to greatest.

    """
    prefix = str(prefix)
    workingpath = Path(__file__).parent
    alldcdpaths = sorted(workingpath.glob(f"{prefix}*.dcd"))
    allstrdcd = [f.name for f in alldcdpaths]
    allstrdcdnocentered = []
    for mn in allstrdcd:
        if 'centered' in mn:
            continue
        else:
            allstrdcdnocentered.append(mn)
    return allstrdcdnocentered


def getdatlistbyrep(replicate="r0"):
    """
    This function will find all of the dat filenames in your working directory
    for a specified replicate.
    
    Parameters
    ----------
    numdecimals : TYPE str
        DESCRIPTION. Replicate suffix of the replicate you are interested in. The default 
        is r0.

    Returns
    -------
    A sorted list of strings with all of the dat filenames from your working directory
    for the replicate of interest. They are sorted from least to greatest.

    """
    workingpath = Path(__file__).parent
    alldatpaths = sorted(workingpath.glob(f"*{replicate}.dat"))
    allstrdat = [f.name for f in alldatpaths]
    return allstrdat


def getdatlistbyprefix(prefix):
    """
    This function will find all of the dat filenames in your working directory 
    for a specified mainname (accept those that are "centered").
    
    Parameters
    ----------
    numdecimals : TYPE str,int,float
        DESCRIPTION. Prefix of the simulation you are interested in.

    Returns
    -------
    A sorted list of strings with all of the dat filenames for the simulation of 
    interest from your working directory. They are sorted from least to greatest.

    """
    prefix = str(prefix)
    workingpath = Path(__file__).parent
    alldatpaths = sorted(workingpath.glob(f"{prefix}*.dat"))
    allstrdat = [f.name for f in alldatpaths]
    allstrdatnocentered = []
    for mn in allstrdat:
        if 'centered' in mn:
            continue
        else:
            allstrdatnocentered.append(mn)
    return allstrdatnocentered


def truncate_float(number, length):
    """
    Truncates a number to a given number of decimals. Used in getprefixlistfloats.

    Parameters
    ----------
    number : TYPE String,float,int
        DESCRIPTION. Number want truncated.
    length : TYPE Int
        DESCRIPTION. The number of decimals you want. If you pass -1, then
        all of the entries will keep the number same number of decimals and will
        just be converted to floats.

    Returns
    -------
    number : TYPE Int,float
        DESCRIPTION. Number with truncated number of decimals.

    """
    if length == -1: # Case where want to keep all decimals
        number = float(number)
    elif '.' in number: # The case when the number is a float
        number = float(number)
        number = number * pow(10, length)
        number = int(number)
        number = float(number)
        number /= pow(10, length)
    else: # The case when the number is an integer (don't want an extra decimal place)
        number = int(number)
    return number


def getprefixlistfloatstruncate(numdecimals):
    """
    This function will find all of the prefixes in your working directory.
    It will truncate the prefixes to the specified numdecimals.

    Parameters
    ----------
    numdecimals : TYPE Int
        DESCRIPTION. Number of decimals you want your prefixes rounded to.

    Returns
    -------
    prefixlistfloatsprecise : TYPE List of floats
        DESCRIPTION. Sorted list of prefixes as floats. Each prefix will be rounded
        to the specified number of decimals.

    """
    # Get a list of strings of pdbs
    pdblist = getpdblist() 
    prefixlist = [file.split("~")[0] for file in pdblist]
    
    # Change the number of decimals
    ###prefixlistfloats = [float(fil) for fil in prefixlist] 
    prefixlistfloatsprecise = [truncate_float(fi, numdecimals) for fi in prefixlist]
    return prefixlistfloatsprecise


def getprefixliststrings():
    """
    This function will find all of the prefixes in your working directory.
    It will not round the prefixes to a specified numdecimals.

    Parameters
    ----------
    numdecimals : TYPE Int
        DESCRIPTION. Number of decimals you want your prefixes rounded to.

    Returns
    -------
    prefixliststrings : TYPE List of strings
        DESCRIPTION. Sorted list of prefixes as strings. Each prefix will be rounded
        to the specified number of decimals.

    """
    # Get the prefix list as a sorted list of floats and convert the to strings
    pdblist = getpdblist() 
    prefixlist = [file.split("~")[0] for file in pdblist]
    prefixliststrings = [str(f) for f in prefixlist]
    return prefixliststrings
    

def replicatesstringslist():
    """
    This function will make a list of replicate suffixes based on the number
    of replicates in the current working directory.

    Returns
    -------
    repsuflist : TYPE list of strings
        DESCRIPTION. List of replicate suffixes based on the number
        of replicates in the current working directory.
        ex. ["r0", "r1"]

    """
    # Count the number of replicates
    workingpath = Path(__file__).parent
    for repnum in range(100): # Arbitrarily choose 100 as did not want to use forever loop and will never have 100 replicates
        if len(list((workingpath.glob(f"*r{repnum}*")))) == 0: # Check to see if there is any replicates with that number in the directory
            numofrep = repnum
            break
    
    # Use the counted number of replicate to make a list of strings of replicate suffixes
    repsuflist = ["r"+str(num) for num in range(numofrep)] 
    return repsuflist



if __name__ == "__main__":
    # Test getmainname()
    print(getprefixliststrings())
    print()

    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
