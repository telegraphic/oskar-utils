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
    """
    def __init__(self):
        super(OskarGui, self).__init__()
        
        # Initialize user interface
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
        self.but_save     = QtGui.QPushButton("Save As")
        self.but_open     = QtGui.QPushButton("Open")
        self.lab_info     = QtGui.QLabel(" ")
        
        
        self.but_stats.clicked.connect(self.onButStats)
        self.but_save.clicked.connect(self.onButSave)
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
        bbox.addWidget(self.but_save)
        layout.addLayout(bbox)
        
        
        self.setLayout(layout)    
        
        self.setGeometry(300, 300, width, height)
        self.setWindowTitle('OSKAR beam pattern tool')    
        self.show()
        
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
        
        
        plt.imshow(beam_data, interpolation='None')
        
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
                
    def onGenGuiClosed(self):
        """ Close action: Generate Station"""
        self.updatePlot()

    def onButSave(self):
        """ Button action: Save to file """
        self.file_dialog    = QtGui.QFileDialog()
        fileparts = self.file_dialog.getSaveFileName(caption="Save antenna configuration", filter="Text files (*.txt *.dat *.csv)")
        filename = fileparts[0]
        
        saveAntConfig(self.ant_coords, filename)        
                      
    def onButOpen(self):
        """ Button action: Open station file """
        self.file_dialog    = QtGui.QFileDialog()
        fileparts = self.file_dialog.getOpenFileName(caption="Select OSKAR antenna configuration file", filter="FITS files (*.fits)")
        filename = fileparts[0]
        
        self.beam_data = openBeamPattern(filename)
        self.updatePlot()

    def onButStats(self):
        """ Button action: Open statistics """
        pass


def main():
    global main_gui
    
    print "Starting OSKAR beam pattern tool..."
    app = QtGui.QApplication(sys.argv)
    main_gui = OskarGui()
    app.exec_()
    sys.exit()    

if __name__ == '__main__':
    main()