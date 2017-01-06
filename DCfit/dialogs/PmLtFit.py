'''
    v0.6.1  
    This is a package  PmLtFit is a master of protein dose curves fit 
    packages.

    Copyright (C) 2015 Baltic Institute of Advanced Technology

    Authors: Sarunas Azna, Piotras Cimmperman, Vytautas Rafanavicius

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''

import psaFit as bp
import qti
import numpy as np
import math
from PyQt4 import QtCore, QtGui

class Window(QtGui.QDialog):
    def __init__(self):
        QtGui.QDialog.__init__(self, qti.app)
        self.setWindowModality(QtCore.Qt.NonModal)
        self.setWindowTitle('Pm (Lt) fit')
        self.initUI()
        self.tablename = ''
        fitwiz = qti.app.table('fitWizard')
        fitwiz.showMaximized()
        fitwiz.showNormal()
        qti.app.resultsLog().append('Dose curves GUI started \n ')

    def initUI(self):

        self.mode = 'N'
        ## model choosing block creation
        rbox = QtGui.QVBoxLayout(self)
        topright = QtGui.QFrame(self)
        topright.setFrameShape(QtGui.QFrame.StyledPanel)
        combo = QtGui.QComboBox(self)
        combo.addItem("N")
        combo.addItem("NU")
        combo.activated[str].connect(self.onActivated)
        self.lbl2 = QtGui.QLabel('Binding model', self)
        rbox.addWidget(self.lbl2)
        rbox.addWidget(combo)
        topright.setLayout(rbox)
        ##action block creation
        vbox = QtGui.QVBoxLayout(self)
        topleft = QtGui.QFrame(self)
        topleft.setFrameShape(QtGui.QFrame.StyledPanel)
        lbl1 = QtGui.QLabel('Action', self)
        tpb1 = QtGui.QPushButton('Fit', self)
        tpb1.clicked.connect(self.pushed)
        '''
        tpb2 = QtGui.QPushButton('Logarithmic differentiate', self)
        tpb2.clicked.connect(self.pushed)
        tpb3 = QtGui.QPushButton('Differentiation parameters', self)
        tpb3.clicked.connect(self.pushed)
        '''
        vbox.addWidget(lbl1)
        vbox.addWidget(tpb1)
        '''
        vbox.addWidget(tpb2)
        vbox.addWidget(tpb3)
        '''
        topleft.setLayout(vbox)
        ## delete block creation
        bbox = QtGui.QVBoxLayout(self)
        bottom = QtGui.QFrame(self)
        bottom.setFrameShape(QtGui.QFrame.StyledPanel)
        bb1 = QtGui.QPushButton('Delete Fit Graphs', self)
        bb1.clicked.connect(bp.deleteFitGraphs)
        bbox.addWidget(bb1)
        # bb2 = QtGui.QPushButton('Delete Fit Tables', self)
        # bb2.clicked.connect(bp.deleteFitTables)
        # bbox.addWidget(bb2)
        bb3 = QtGui.QPushButton('Delete Fit Graphs and Tables', self)
        bb3.clicked.connect(bp.deleteFitGraphsTables)
        bbox.addWidget(bb3)
        bottom.setLayout(bbox)
        ## splitting widget
        hbox = QtGui.QHBoxLayout(self)
        splitter1 = QtGui.QSplitter(QtCore.Qt.Horizontal)
        splitter1.addWidget(topleft)
        splitter1.addWidget(topright)
        splitter2 = QtGui.QSplitter(QtCore.Qt.Vertical)
        splitter2.addWidget(splitter1)
        splitter2.addWidget(bottom)
        hbox.addWidget(splitter2)
        self.setLayout(hbox)
        QtGui.QApplication.setStyle(QtGui.QStyleFactory.create('Cleanlooks'))
        self.setGeometry(300, 300, 300, 200)
        self.show()

    def onActivated(self, text):
        self.mode = text
    def pushed(self):
        self.m = self.sender()
        if self.m.text() == 'Fit':
            if not bp.objectMessanger('fitWizard', 'Table'):
                return
            t = qti.app.currentTable()
            if isinstance(t, qti.Table):
                name = t.objectName()
                if name == 'fitWizard':
                    t = qti.app.table(self.tablename)
            if  not isinstance(t, qti.Table):
                t = qti.app.table(self.tablename)
                if not isinstance(t, qti.Table):
                    msgBox = QtGui.QMessageBox()
                    msgBox.setText("Please select data table")
                    ret = msgBox.exec_()
                    return
            try:
                name = t.objectName()
            except Exception, e:
                msgBox = QtGui.QMessageBox()
                msgBox.setText(str(e))
                ret = msgBox.exec_()
            self.tablename = name
            if self.mode == 'N':
                qti.app.resultsLog().append('N-model fit started \n ')
                bp.fitPmLt(name, "N")
            elif self.mode == 'NU':
                qti.app.resultsLog().append('NU-model fit started \n ')
                bp.fitPmLt(name, 'NU')
        
win = Window()
win.move(600, 400)

win.show()
