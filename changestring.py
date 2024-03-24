#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 25 16:15:00 2023

@author: joseph
"""
# This script changes as many unique strings in a file as you want

import os
from pathlib import Path
import numpy as np

def multipliers(originalfloats, multipliers):
    """
    This function will take in a list of floats and a list of corresponding
    multipliers for those floats. Then, it will print out a list with all of original floats
    multiplied by the corresponding multipliers. Each row of the outputted array
    will be the list of originalfloats all multiplied by the same multiplier.

    Parameters
    ----------
    originalfloats : TYPE list of floats
        DESCRIPTION. List of original floats want multiplied.
    multipliers : TYPE list of floats
        DESCRIPTION. List of corresponding multipliers for the entire list of original floats.

    Returns
    -------
    None.
    Prints list of lists with each column being the original floats list multiplied by the corresponding
    multiplier for that particular column. Also, print originalfloats and multipliers of list of strings.
    The list of lists is also printed as only a single list. In this case only the first string in each 
    list inside of the list of lists is printed.

    """
    # Do calculations and make into list of lists
    originalfloatsarray = np.array(originalfloats, float)
    newarray = np.empty((0,len(originalfloats)))
    for mult in multipliers:
        newarray = np.vstack((newarray, (originalfloatsarray*mult))) 
    strlist = []
    for row in newarray:
        strlist.append(list(map(str, row)))

    # Print the list of lists so it is only a list
    # Useful for when want the fileprefixlist to be the actual values the file was changed to
    # If multiple values are changed do not use this output
    print("ONLY USE THIS FUNCTION IF CHANGING BY A MULTIPLIER!!!")
    print()
    print("when want actual values as prefixes (use when only changing one value) fileprefixlist = ")
    for inx,ro in enumerate(strlist):
        if inx == 0:
            print(f"['{ro[0]}'", end=", ")
        elif inx == len(strlist) - 1:
            print(f"'{ro[0]}']")
        else:
            print(f"'{ro[0]}'", end=", ")
    print()
    # Adjust multipliers so that it is a list of strings and print this out
    strmultipliers = list(map(str, multipliers))
    print(f"when want multipliers as prefixes (use when changing multiple values by same multiplier) fileprefixlist = \n{strmultipliers}")
    print()

    # Adjust originalstrings so that it is a list of strings and print this out
    stroriginalfloats = list(map(str, originalfloats))
    print(f"cg/runoriginalstrings = \n{stroriginalfloats}")
    print()

    # Print the list of lists in a format so it can be easily copied and read
    for ind,myrow in enumerate(strlist):
        if ind == 0:
            print(f"cg/runnewstrings = [{myrow},")
        elif ind == len(strlist) - 1:
            print(f"                {myrow}]")
        else:
            print(f"                {myrow},")
    print()
    


def changeonefloateachfilelists(originalfloat, newfloatslist):
    """
    This function will take in a list with a single original float and a list of new floats.
    It will print the single original float as a list with one string. It will print 
    the newfloats list as a list of lists and a list of strings. Use this when
    want to find strings to change to by using a hard coded list of floats directly.
    Parameters
    ----------
    originalfloats : TYPE List of float
        DESCRIPTION. Single float in a list.
    newfloatslist : TYPE List of floats
        DESCRIPTION. List of floats.

    Returns
    -------
    None.
    Prints the inputs for fileprefixlist, cg/runoriginalstring and cg/runnewstrings.
    
    """
    print("ONLY USE THIS CODE WHEN CHANGING ONE STRING IN EACH cg.ff DIRECTLY FROM A LIST THAT YOU HARD CODE IN!!!")
    
    # Make list of string for originalfloat
    stroriginalfloat = list(map(str, originalfloat))
    print()
    
    # Make list of lists for cg/runnewstrings
    strnewfloatslist = list(map(str, newfloatslist))
    listoflistnewstrings = []
    for stri in strnewfloatslist:
        listoflistnewstrings.append([stri])
        
    # Print the list of lists as a list of strings so get the fileprefixlist
    print("fileprefixlist = ")
    for inx,ro in enumerate(listoflistnewstrings):
        if inx == 0:
            print(f"['{ro[0]}'", end=", ")
        elif inx == len(listoflistnewstrings) - 1:
            print(f"'{ro[0]}']")
        else:
            print(f"'{ro[0]}'", end=", ")
    print()
    
    # Print the original float list as a string list
    print(f"cg/runoriginalstring = \n{stroriginalfloat}")
    print()
    
    # Print the list of lists in a format so it can be easily copied and read
    for ind,myrow in enumerate(listoflistnewstrings):
        if ind == 0:
            print(f"cg/runnewstrings = [{myrow},")
        elif ind == len(listoflistnewstrings) - 1:
            print(f"                {myrow}]")
        else:
            print(f"                {myrow},")
    print()



def editfilenewstrings(filename, originalstrings, newstrings, outputpath=Path.cwd(), filenamenew=None):
    """
    This function will take in a list of original strings, and replace them with 
    the strings in the newstrings list of lists. Then, it will insert these into a copy 
    of the given file. All occurances of the original string will be replaced so
    do not put in the same argument twice.

    Parameters
    ----------
    filename : TYPE str
        DESCRIPTION. The name of the file to be edited.
    originalstrings : TYPE list of strings
        DESCRIPTION. List of floats to be replaced in file. They must all be unique. 
    newstrings : TYPE list of strings
        DESCRIPTION. List of strings to be used to replace corresponding originalstrings.
    outputpath : TYPE str
        DESCRIPTION. The desired outputpath. The default is the current working directory.
        Omit the actual name of the file. This will be added elsewhere.
    filenamenew : TYPE str
        DESCRIPTION. The desired filename for the output. If you do not want to change the 
        files actual name then omit this.
    Returns
    -------
    None.
    Outputs copy of the specified file that is editted accordingly.

    """
    # Open original file and import it as text
    with open(filename, 'r') as myfile:
        text = myfile.read()
        # Make sure the string your are replacing only occurs once or at least once
        for ori in originalstrings:
            oricount = str(text.count(str(ori)))
            if int(oricount) > 1:
                print(f"A string occurs more then once. The output file will replace all of those {oricount} occurances.")
            if int(oricount) == 0:
                raise Exception("A string is not in the file. Check your command for spelling mistakes.")

    # Change the imported text and output it to the desired path
    # Will overwrite files if same name and path
    if filenamenew == None:
        filenamenew = filename
    mydir = os.getcwd()
    os.chdir(outputpath)
    f = open(filenamenew, "w")
    for ind, old in enumerate(originalstrings): 
        text = text.replace(old, newstrings[ind])
    f.write(text)
    f.close()
    os.chdir(mydir)



