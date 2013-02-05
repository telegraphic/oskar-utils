#!/usr/bin/env python
"""
beam_pattern_viewer.py
====================

This script starts a graphical user interface for viewing OSKAR beam pattern.

Requirements
------------
Qt4
PySide 1.1 
numpy, matplotlib (1.1)
pyFITS

"""

__version__ = "0.1"
__author__  = "Danny Price"

# Imports
import sys
from optparse import OptionParser

try:
    from PySide import QtCore, QtGui
except:
    print "Error: cannot load PySide. Please check your install."
    exit()
    
try:    
    import numpy as np
except:
    print "Error: cannot load Numpy. Please check your install."
    exit()

try:    
    import pyfits as pf
except:
    print "Error: cannot load PyFITS. Please check your install."
    exit()

import matplotlib
if matplotlib.__version__ == '0.99.3':
    print "Error: your matplotlib version is too old to run this. Please upgrade."
    exit()
else:
    matplotlib.use('Qt4Agg')
    matplotlib.rcParams['backend.qt4']='PySide'
    from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
    from matplotlib.figure import Figure
try:
    import pylab as plt
except:
    print "Error: cannot load Pylab. Check your matplotlib install."
    exit()

from lib.beam_pattern_utils import *

    

class OskarGui(QtGui.QWidget):
    """ OSKAR GUI class
    
    A Qt4 Widget that uses matplotlib to display antenna configurations
    
    Parameters
    ----------
    filename: str
        Name of file to open. Defaults to blank, in which case no file is opened.
    """
    def __init__(self, filename=''):
        super(OskarGui, self).__init__()
        
        # Initialize user interface
        self.filename = filename
        self.initUI(width=800, height=800)
        
       

    def initUI(self, width=1024, height=768):
        """ Initialize the User Interface 
        
        Parameters
        ----------
        width: int
            width of the UI, in pixels. Defaults to 1024px
        height: int
            height of the UI, in pixels. Defaults to 768px
        """
        
        self.main_frame = QtGui.QWidget()    
        #self.gen_gui = generateGui()
        
        # Create buttons/widgets
        self.but_stats    = QtGui.QPushButton("Statistics")
        #self.but_save     = QtGui.QPushButton("Save As")
        self.but_open     = QtGui.QPushButton("Open")
        self.lab_info     = QtGui.QLabel(" ")
        
        
        self.but_stats.clicked.connect(self.onButStats)
        #self.but_save.clicked.connect(self.onButSave)
        self.but_open.clicked.connect(self.onButOpen)
        
        self.beam_data = np.zeros([5,5])
        
        # Create plots
        self.sb_fig, self.sb_ax = self.createBeamPlot()
        
        
        # generate the canvas to display the plot
        self.sb_canvas = FigureCanvas(self.sb_fig)
        self.mpl_toolbar = NavigationToolbar(self.sb_canvas, self.main_frame)
        
        # Widget layout
        layout = QtGui.QVBoxLayout()
        
        
        layout.addWidget(self.sb_canvas)
        layout.addWidget(self.mpl_toolbar)
        
        bbox = QtGui.QHBoxLayout()
        bbox.addWidget(self.lab_info)
        bbox.addStretch(1)
        bbox.addWidget(self.but_open)
        bbox.addWidget(self.but_stats)
        #bbox.addWidget(self.but_save)
        layout.addLayout(bbox)
        
        
        self.setLayout(layout)    
        
        self.setGeometry(300, 300, width, height)
        self.setWindowTitle('OSKAR beam pattern tool')    
        self.show()
        
        # Load file if command line argument is passed
        if self.filename != '':
            try:
                self.beam_data = openBeamPattern(self.filename)
                self.updatePlot()
            except:
                print "Error: cannot open %s"%self.filename
        
        
    def createBeamPlot(self):
          """ Creates a single pylab plot for displaying a beam pattern """

          fig = plt.figure(figsize=(8,6),dpi=80)
          fig.set_facecolor('#ededed')
          
          # Format plot
          ax = plt.subplot(111)

          ax.set_xlabel("Theta (deg)")
          ax.set_ylabel("Phi (deg)")
          ax.set_xlim(0,1)
          ax.set_ylim(0,1)
          plt.xticks(np.linspace(0,1,5), [-90, -45, 0, 45, 90])
          plt.yticks(np.linspace(0,1,5), [-90, -45, 0, 45, 90])
          
        
          fig.canvas.draw()
      
          return fig, ax

    def updatePlot(self):
        """ Updates the antenna config plot"""
                
        self.sb_ax.clear()
        self.sb_fig.clear()
        
        xdim, ydim = self.beam_data[0].data.shape[3:]
        beam_data = self.beam_data[0].data[0,0,0]
        
        
        self.sb_plot = plt.imshow(beam_data, interpolation='None')
        
        plt.xlabel(self.beam_data[0].header.get('CTYPE1'))
        plt.ylabel(self.beam_data[0].header.get('CTYPE2'))
        plt.xticks(np.linspace(0,xdim,5), [-90, -45, 0, 45, 90])
        plt.yticks(np.linspace(0,xdim,5), [-90, -45, 0, 45, 90])
        plt.colorbar()
        
        self.sb_fig.canvas.draw()
        
        ra  = float(self.beam_data[0].header.get('OBSRA'))
        dec = float(self.beam_data[0].header.get('OBSDEC'))
        freq = float(self.beam_data[0].header.get('CRVAL5')) / 1e6
        
        self.lab_info.setText("RA: %2.2f, DEC: %2.2f, FREQ: %2.2fMHz"%(ra, dec, freq))       

    def onButSave(self):
        """ Button action: Save to file 
        TODO: NOT CURRENTLY IMPLEMENTED
        """
        self.file_dialog    = QtGui.QFileDialog()
        fileparts = self.file_dialog.getSaveFileName(caption="Save beam pattern", filter="FITS files (*.fits)")
        filename = fileparts[0]
        
        saveAntConfig(self.ant_coords, filename)        
                      
    def onButOpen(self):
        """ Button action: Open station file """
        self.file_dialog    = QtGui.QFileDialog()
        fileparts = self.file_dialog.getOpenFileName(caption="Select OSKAR beam pattern file", filter="FITS files (*.fits)")
        filename = fileparts[0]
        
        self.beam_data = openBeamPattern(filename)
        self.updatePlot()

    def onButStats(self):
        """ Button action: Open statistics """
        beam_data  = self.beam_data[0].data[0,0,0]
        xdim, ydim = self.beam_data[0].data.shape[3:]
        
        #print "\nTesting findMaxima()"
        max_loc = findMaxima(beam_data, size=10, threshold=0.4)
        #print "%i local maxima found"%len(max_loc)    

        ml_coords = findMainLobe(beam_data)
        (fov_ra, fov_dec) = findFieldOfView(self.beam_data)
        (fwhm_x, fwhm_y) = findFWHM(beam_data)
        fwhm_ra  = fwhm_x / xdim * fov_ra
        fwhm_dec = fwhm_y / ydim * fov_dec
        
        print "Main lobe coords: %s"%ml_coords
        print "FoV  in RA: %2.2fdeg, DEC: %2.2f"%(fov_ra, fov_dec)
        print "FWHM in RA: %2.2fdeg, DEC: %2.2f"%(fwhm_ra, fwhm_dec)
        
        # Overlay location of maxima
        stat_plot = self.sb_fig.add_subplot(111)
        stat_plot.plot(max_loc[:,0],max_loc[:,1], 'ro')
        stat_plot.plot(ml_coords[0], ml_coords[1], marker='x', markersize=16, markeredgewidth=1.5, color='#ffffff')
    
        # Draw box around FWHM
        stat_plot.axvline(ml_coords[0]-fwhm_x/2, color='#FFFFFF', linestyle='dashed')
        stat_plot.axvline(ml_coords[0]+fwhm_x/2, color='#FFFFFF', linestyle='dashed')
        stat_plot.axhline(ml_coords[1]-fwhm_y/2, color='#FFFFFF', linestyle='dashed')
        stat_plot.axhline(ml_coords[1]+fwhm_y/2, color='#FFFFFF', linestyle='dashed')
    
        stat_plot.set_xlim(0,beam_data.shape[0])
        stat_plot.set_ylim(0,beam_data.shape[1])
        self.sb_fig.canvas.draw()


def main():
    
    # Basic option parsing 
    p = OptionParser()
    p.set_usage('beam_pattern_viewer.py [filename] [options]')
    p.set_description(__doc__)
    (options, args) = p.parse_args()

    print "Starting OSKAR beam pattern tool..."
    global main_gui
    app = QtGui.QApplication(sys.argv)
    
    try:
        filename = args[0]
        main_gui = OskarGui(filename)
    except:
        main_gui = OskarGui()
    app.exec_()
    sys.exit()    

if __name__ == '__main__':
    main()