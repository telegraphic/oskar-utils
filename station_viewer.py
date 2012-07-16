#!/usr/bin/env python
"""
station_viewer.py
====================

This script starts a graphical user interface for viewing and generating OSKAR stations.

Requirements
------------
Qt4
PySide 1.1 
numpy, matplotlib (1.1)

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

from lib.antenna_config_utils import *

class OskarGui(QtGui.QWidget):
    """ OSKAR GUI class
    
    A Qt4 Widget that uses matplotlib to display antenna configurations
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
        
        self.gen_gui = generateGui()
        
        # Create buttons/widgets
        self.but_gen      = QtGui.QPushButton("Generate")
        self.but_save     = QtGui.QPushButton("Save As")
        self.but_open     = QtGui.QPushButton("Open")
        self.lab_nant     = QtGui.QLabel("No. antennas: 0")
        
        
        self.but_gen.clicked.connect(self.onButGen)
        self.but_save.clicked.connect(self.onButSave)
        self.but_open.clicked.connect(self.onButOpen)
        
        self.ant_coords = []
        
        # Create plots
        self.sb_fig, self.sb_ax, self.sb_title = self.createStationPlot()
        
        
        # generate the canvas to display the plot
        self.sb_canvas = FigureCanvas(self.sb_fig)
        self.mpl_toolbar = NavigationToolbar(self.sb_canvas, self.main_frame)
        
        # Widget layout
        layout = QtGui.QVBoxLayout()
        
        
        layout.addWidget(self.sb_canvas)
        layout.addWidget(self.mpl_toolbar)
        
        bbox = QtGui.QHBoxLayout()
        bbox.addWidget(self.lab_nant)
        bbox.addStretch(1)
        bbox.addWidget(self.but_open)
        bbox.addWidget(self.but_gen)
        bbox.addWidget(self.but_save)
        layout.addLayout(bbox)
        
        
        self.setLayout(layout)    
        
        self.setGeometry(300, 300, width, height)
        self.setWindowTitle('OSKAR station tool')    
        self.show()
        
        # Load file if command line argument is passed
        if self.filename != '':
            try:
                self.ant_coords = openAntConfig(self.filename)
                self.updatePlot()
            except:
                print "Error: cannot open %s"%self.filename
        
        
    def createStationPlot(self):
          """ Creates a single pylab plot for antenna layout """

          fig = plt.figure(figsize=(8,6),dpi=80)
          fig.set_facecolor('#ededed')
          title = fig.suptitle("Antenna Positions")
          title.set_fontsize(14)
      
          # Format plot
          ax = plt.subplot(111)
          ax.set_ylim(-100,100)
          ax.set_xlim(-100,100)
          ax.set_xlabel("X Position (m)")
          ax.set_ylabel("Y Position (m)")
          ax.grid()  
        
          fig.canvas.draw()
      
          return fig, ax, title

    def updatePlot(self):
        """ Updates the antenna config plot"""
        
        x,y = self.ant_coords[:,0], self.ant_coords[:,1]
        nants = len(x)
        
        self.sb_ax.clear()
        self.sb_ax.plot(x, y, marker='.', markersize=80.0/np.sqrt(nants), lw=0, color='black')
        self.sb_ax.set_xlim(np.min(x) - 1, np.max(x) + 1)      
        self.sb_ax.set_ylim(np.min(y) - 1, np.max(y) + 1)
        self.sb_ax.grid()
        self.sb_ax.set_xlabel("X Position (m)")
        self.sb_ax.set_ylabel("Y Position (m)")
        self.lab_nant.setText("No. antennas: %i"%len(x))
        self.sb_fig.canvas.draw()       
        
    def onButGen(self):
        """ Button action: Generate Station"""
        self.gen_gui.show()
                
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
        fileparts = self.file_dialog.getOpenFileName(caption="Select OSKAR antenna configuration file", filter="Text files (*.txt *.dat *.csv)")
        filename = fileparts[0]
        
        self.ant_coords = openAntConfig(filename)
        self.updatePlot()


class generateGui(QtGui.QWidget):
    """ Antenna configuration generator
    
    A widget that generates some basic antenna configurations. 
    """
    
    global main_gui     # Sure there is a better way, but this works
    
    def __init__(self):
        super(generateGui, self).__init__()
        
        # Initialize user interface
        self.initUI(width=400, height=600)
        
        
    def initUI(self, width, height):
        """ Initialize the generate form """
        
        self.ant_coords = []
    
        
        # Create widgets
        self.ant_spacing = QtGui.QLineEdit("0.5")
        self.ant_xnum    = QtGui.QLineEdit("4")
        self.ant_ynum    = QtGui.QLineEdit("4")
        self.but_gen     = QtGui.QPushButton("Generate")
        self.check_taper = QtGui.QCheckBox("Apply taper")
                
        self.but_gen.clicked.connect(self.onButGen)
        
    
        # Widget layout
        layout = QtGui.QVBoxLayout()
        cbox = QtGui.QHBoxLayout()
        cbox.addWidget(QtGui.QLabel("Antenna Spacing (m)"))
        cbox.addWidget(self.ant_spacing)
        layout.addLayout(cbox)
                
        cbox = QtGui.QHBoxLayout()
        cbox.addWidget(QtGui.QLabel("Dimensions"))
        cbox.addWidget(self.ant_xnum)
        cbox.addWidget(QtGui.QLabel("x"))
        cbox.addWidget(self.ant_ynum)
        layout.addLayout(cbox)
        
        layout.addWidget(self.check_taper)
        
        cbox = QtGui.QHBoxLayout()
        cbox.addStretch(1)
        cbox.addWidget(self.but_gen)
        layout.addLayout(cbox)

        self.setLayout(layout)    
        
        self.setGeometry(300, 300, width, height)
        self.setWindowTitle('OSKAR Station Generator')    
        #self.show()
        
    def onButGen(self):
        """ Button action: Generate Station"""
        
        ant_xnum = int(self.ant_xnum.text())
        ant_ynum = int(self.ant_ynum.text())
        ant_spacing = float(self.ant_spacing.text())        
        
        ant_coords = generateGrid(ant_xnum, ant_ynum, ant_spacing)
        if self.check_taper.isChecked():
            ant_coords = applyTaper(ant_coords, float(ant_xnum)/2 * ant_spacing) 
                                 
        main_gui.ant_coords = ant_coords
        main_gui.updatePlot()
        self.close()
        

def main():
    # Basic option parsing 
    p = OptionParser()
    p.set_usage('beam_pattern_viewer.py [filename] [options]')
    p.set_description(__doc__)
    (options, args) = p.parse_args()
    
    print "Starting OSKAR antenna config tool..."
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