#!/usr/bin/env python
"""
antenna_config_utils.py
=======================

A collection of functions for beam patterns of OSKAR stations.

"""

__version__ = "0.1"
__author__  = "Danny Price"

import sys
import numpy as np
import pyfits as pf

def openBeamPattern(filename):
    """ Opens an OSKAR beam pattern file using PyFITS
    
    Parameters
    ----------
    filename: str
        path to file
    """
    try:
        beam = pf.open(filename)
        return beam  
    except:
        print "Error: cannot open %s"%filename
        raise
      
    
    



    

