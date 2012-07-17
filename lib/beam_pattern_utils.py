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
import scipy.ndimage as ndimage
import scipy.ndimage.filters as filters

import matplotlib
import pylab as plt   

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


def findMaxima(beam_data, threshold=0.5, size=10):
    """ Find the local maxima of an array
    
    This function uses filters.maximum_filter and minimum_filter from scipy.ndimage. Based off
    http://stackoverflow.com/questions/9111711/get-coordinates-of-local-maxima-in-2d-array-above-certain-value
    
    Parameters
    ----------
    beam_data: np.array (2D)
        beam data in which to search for maxima
    threshold: float
        threshold value over which (max-min) yields a detection
    size: scalar or tuple
        neighbourhood size gives the shape that is taken from the input array, 
        at every element position, to define the input to the filter function.
    
    """
    
    data = np.nan_to_num(beam_data)
    
    # Apply scipy max filter

    data_max = filters.maximum_filter(data, size)
    maxima = (data == data_max)
    data_min = filters.minimum_filter(data, size)
    diff = ((data_max - data_min) > threshold)
    maxima[diff == 0] = 0

    labeled, num_objects = ndimage.label(maxima)
    slices = ndimage.find_objects(labeled)
    
    x, y = [], []
    for dy,dx in slices:
        x_center = (dx.start + dx.stop - 1)/2
        x.append(x_center)
        y_center = (dy.start + dy.stop - 1)/2    
        y.append(y_center) 
    
    return np.column_stack((x,y))

def findMainLobe(beam_data):
    """ Find the main lobe of a beam 
    
    This function assumes the main lobe is the maxima closest to the centre of the beam.
    Returns the coordinates of the main lobe (xpixel, ypixel)
    
    Parameters
    ----------
    beam_data: np.array (2D)
        beam pattern data array
    
    """
    xdim, ydim = beam_data.shape
    
    # Find all local maxima
    max_loc = findMaxima(beam_data, size=10, threshold=0.5)
    
    # Compute the distance from centre for each local maxima
    dist = np.sum(np.abs(max_loc - np.array([xdim/2, ydim/2])), axis=1)
    
    # Return the local maxima closest to the centre
    return max_loc[np.argmin(dist)]
    
    

def findFWHM(beam_data):
    """ Find the full width half maximum of a beam pattern
    
    returns the fwhm over x and y axes, as a tuple (fwhm_x, fwhm_y)
    
    Parameters
    ----------
    beam_data: np.array (2D)
        beam pattern data array
    """
    
    xdim, ydim = beam_data.shape
    ml_coords = findMainLobe(beam_data)
    ml_max = beam_data[ml_coords[0], ml_coords[1]]
    ml_3db = ml_max / 2
    
    w = 1
    while ml_3db <= beam_data[ml_coords[0]+w, ml_coords[1]]:
        w += 1
        if w + ml_coords[0] -1 >= xdim:
            break
    x1 = w
    
    w = -1
    while ml_3db <= beam_data[ml_coords[0]-w, ml_coords[1]]:
        w -= 1
        if w + ml_coords[0] -1 >= xdim:
            break
    x2 = w
    
    w = 1
    while ml_3db <= beam_data[ml_coords[0], ml_coords[1]+w]:
        w +=1
        if w + ml_coords[0] -1 >= xdim:
            break
    y1 = w

    w = -1
    while ml_3db <= beam_data[ml_coords[0], ml_coords[1]-w]:
        w -=1
        if w + ml_coords[0] -1 >= xdim:
            break
    y2 = w    
    
    fwhm_x = np.average([np.abs(x1), np.abs(x2)])
    fwhm_y = np.average([np.abs(y1), np.abs(y2)])
    
    return (fwhm_x, fwhm_y)

    
    
def run_tests():
    """ Some test routines for stats """

    print "Running tests/examples..."
    
    # Load test data
    print "\nTesting openBeamPattern()"
    beam = openBeamPattern('beampattern.fits')
    xdim, ydim = beam[0].data.shape[3:]
    beam_data  = beam[0].data[0,0,0]
    print "xdim: %s, ydim: %s"%(xdim, ydim)

    print "\nTesting findMaxima()"
    max_loc = findMaxima(beam_data, size=10, threshold=0.4)    
    
    print "\nTesting findMainLobe()"
    ml_coords = findMainLobe(beam_data)
    print "Main lobe coordinates: %s"%ml_coords

    print "\nTesting findFWHM()"
    fwhm = findFWHM(beam_data)
    print "FWHM in x-dir: %2.2f, FWHM in y-dir %2.2f"%fwhm

    print "\nPlotting results..."
    plt.imshow(beam_data)
    plt.colorbar()
    
    # Overlay location of maxima
    plt.plot(max_loc[:,0],max_loc[:,1], 'ro')
    plt.plot(ml_coords[0], ml_coords[1], marker='x', markersize=16, markeredgewidth=1.5, color='#ffffff')
    
    # Draw box around FWHM
    plt.axvline(ml_coords[0]-fwhm[0], color='#FFFFFF', linestyle='dashed')
    plt.axvline(ml_coords[0]+fwhm[0], color='#FFFFFF', linestyle='dashed')
    plt.axhline(ml_coords[1]-fwhm[1], color='#FFFFFF', linestyle='dashed')
    plt.axhline(ml_coords[1]+fwhm[1], color='#FFFFFF', linestyle='dashed')
    
    plt.xlim(0,beam_data.shape[0])
    plt.ylim(0,beam_data.shape[1])
    plt.show()



if __name__ =='__main__':
    run_tests()

        
  
      
    
    



    

