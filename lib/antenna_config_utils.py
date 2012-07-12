#!/usr/bin/env python
"""
antenna_config_utils.py
=======================

A collection of functions for generating antenna configurations for OSKAR stations.

"""

__version__ = "0.1"
__author__  = "Danny Price"

import sys
import numpy as np

def applyTaper(ant_coords, radius):
    """ Apply a taper to convert a square antenna grid to be circular
    
    Parameters
    ----------
    ant_coords: numpy.array
        an Nx2 numpy array of (x,y) antenna coordinates
    radius: float
        a radial distance under which antennas are not removed from grid
    """
    
    x, y = ant_coords[:,0], ant_coords[:,1]
    z = np.sqrt(x**2 + y**2)

    return ant_coords[z<radius]

def generateGrid(ant_xnum, ant_ynum, ant_spacing):
    """ Generate a rectangular grid of antennas 
    
    returns a numpy array of (x,y) antenna coordinates
    
    Parameters
    ----------
    ant_xnum: int
        number of antennas across in x-direction
    ant_ynum: int
        number of antennas across in y-direction
    ant_spacing: float
        spacing between antennas, in metres.
    """
    
    list_x, list_y = [], []
    
    # This is a slow implementation, probably a better way to do this than a nested loop
    for x in range(0, ant_xnum):
        for y in range(0, ant_ynum):
            pos_x = ant_spacing * x
            pos_y = ant_spacing * y
            list_x.append(pos_x)
            list_y.append(pos_y)
    
    # Centre on (0,0) - this works for both odd and even
    x = np.array(list_x) - float(ant_xnum)/2  * ant_spacing + 0.5 * ant_spacing
    y = np.array(list_y) - float(ant_xnum)/2  * ant_spacing + 0.5 * ant_spacing
        
    ant_coords = np.column_stack((x, y))
    
    return ant_coords

def generateStations(station_file):
    """ Convert an OSKAR station file into multiple stations of single antennas 
    
    TODO: make this better - output directory path choice, generate new output dir
    
    Parameters
    ----------
    station_file: str
        path to antenna station configuration file (normally config.txt)
    """
    
    # Generate a list of antennas from the station file
    try:
        antennas = openAntConfig(station_file)
    except:
        print 'Error: Could not open station file'
        raise
    
    # Loop through the antenna list, creating a new folder and station for each antenna
    try:
        for i in range(0,len(antennas)):
            station_id = 'station%03i'%(i+1)
            print 'creating %s/config.txt'%station_id
            os.mkdir(station_id)
            f = open(os.path.join(station_id, 'config.txt'), 'w')
            f.write('0,0,0\n\n')
            f.close()
    
    except:
        print 'Error: could not generate stations'
        raise

def openAntConfig(filename):
    """ Open an antenna configuration file
    
    Parameters
    ----------
    filename: str
        path to file
    """   
    
    return np.genfromtxt(filename, delimiter=',')

def saveAntConfig(ant_coords, filename, mode='w'):
    """ Save an antenna configuration to file
    
    Parameters
    ----------
    ant_coords: numpy.array
        array of antenna coordinates
    filename: str
        path to file
    mode: str
        mode to open file. Defaults to 'w'
    """   
    
    f = open(filename, mode)
    for row in ant_coords:
        f.write("%s,%s\n"%(row[0], row[1]))
    f.write("\n")   # OSKAR needs trailing newline                
    f.close()
    print "Antenna configuration written to %s"%filename


    

