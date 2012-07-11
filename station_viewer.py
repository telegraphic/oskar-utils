#!/usr/bin/env python
"""
station_viewer.py
====================

This script starts a graphical user interface for viewing and generating OSKAR stations.

Requirements
------------
PySide for Qt4 bindings, numpy, matplotlib.

"""

__version__ = "0.1"
__author__  = "Danny Price"

# Imports
import sys
from optparse import OptionParser

from PySide import QtCore, QtGui

import numpy as np

import matplotlib
matplotlib.use('Qt4Agg')
matplotlib.rcParams['backend.qt4']='PySide'
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.gridspec as gridspec
import pylab as plt

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
        sb_canvas = FigureCanvas(self.sb_fig)
        
        # Widget layout
        layout = QtGui.QVBoxLayout()
        
        layout.addWidget(sb_canvas)
        
        bbox = QtGui.QHBoxLayout()
        bbox.addWidget(self.lab_nant)
        bbox.addStretch(1)
        bbox.addWidget(self.but_open)
        bbox.addWidget(self.but_gen)
        bbox.addWidget(self.but_save)
        layout.addLayout(bbox)

        self.setLayout(layout)    
        
        self.setGeometry(300, 300, width, height)
        self.setWindowTitle('OSKAR Station Generator')    
        self.show()
        

    ################################
    ##    ANTENNA CONFIG PLOTS    ##
    ################################
        
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



    ##########################
    ##    BUTTON ACTIONS    ##
    ##########################
        
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
        
        f = open(filename, 'w')
        for row in self.ant_coords:
            f.write("%s,%s\n"%(row[0], row[1]))
        f.write("\n")   # OSKAR needs trailing newline                
        f.close()
        print "Antenna configuration written to %s"%filename
        
                      
    def onButOpen(self):
        """ Button action: Open station file """
        self.file_dialog    = QtGui.QFileDialog()
        fileparts = self.file_dialog.getOpenFileName(caption="Select OSKAR antenna configuration file", filter="Text files (*.txt *.dat *.csv)")
        filename = fileparts[0]
        
        self.ant_coords = np.genfromtxt(filename, delimiter=',')
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
        
        # This loop is going to be slow - redo with numpy at some stage
        list_x, list_y = [], []
        ant_xnum = int(self.ant_xnum.text())
        ant_ynum = int(self.ant_ynum.text())
        ant_spacing = float(self.ant_spacing.text())
        
        
        for x in range(0, ant_xnum):
            for y in range(0, ant_ynum):
                pos_x = ant_spacing * x
                pos_y = ant_spacing * y
                list_x.append(pos_x)
                list_y.append(pos_y)
        
        main_gui.ant_coords = np.column_stack((list_x, list_y))
        main_gui.updatePlot()
        self.close()
        

def main():
    global main_gui
    
    print "Starting OSKAR Station Generator..."
    app = QtGui.QApplication(sys.argv)
    main_gui = OskarGui()
    app.exec_()
    sys.exit()    

if __name__ == '__main__':
    main()
