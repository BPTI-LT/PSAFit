'''
    v0.6.1  
    This is a package  PmLtGen is a master of protein dose curves generation
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
        self.setWindowTitle('Dose curves analysis')
        self.initUI()

    def initUI(self):

        self.mode = 'N'
        self.amplitude = 0
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
        tpb1 = QtGui.QPushButton('Generate', self)
        tpb1.clicked.connect(self.pushed)
        tpb2 = QtGui.QPushButton('Logarithmic differentiate', self)
        tpb2.clicked.connect(self.pushed)
        tpb3 = QtGui.QPushButton('Differentiation parameters', self)
        tpb3.clicked.connect(self.pushed)
        vbox.addWidget(lbl1)
        vbox.addWidget(tpb1)
        vbox.addWidget(tpb2)
        vbox.addWidget(tpb3)
        topleft.setLayout(vbox)
        # noise amplitude place creation
        mrbox = QtGui.QVBoxLayout(self)
        newright = QtGui.QFrame(self)
        newright.setFrameShape(QtGui.QFrame.StyledPanel)
        qle3 = QtGui.QLineEdit(self)
        qle3.move(100, 100)
        qle3.textChanged[str].connect(self.noiseAmplitude)
        lbl4 = QtGui.QLabel('<i>Noise amplitude</i>, MPa', self)
        mrbox.addWidget(lbl4)
        mrbox.addWidget(qle3)
        newright.setLayout(mrbox)
        ## delete block creation
        bbox = QtGui.QVBoxLayout(self)
        bottom = QtGui.QFrame(self)
        bottom.setFrameShape(QtGui.QFrame.StyledPanel)
        bb1 = QtGui.QPushButton('Delete Graphs', self)
        bb1.clicked.connect(bp.deleteGraphs)
        bbox.addWidget(bb1)
        bb2 = QtGui.QPushButton('Delete Tables', self)
        bb2.clicked.connect(bp.deleteTables)
        bbox.addWidget(bb2)
        bb3 = QtGui.QPushButton('Delete Graphs and Tables', self)
        bb3.clicked.connect(bp.deleteGraphsTables)
        bbox.addWidget(bb3)
        bottom.setLayout(bbox)
        ## splitting widget
        hbox = QtGui.QHBoxLayout(self)

        splitter3 = QtGui.QSplitter(QtCore.Qt.Horizontal)
        splitter3.addWidget(newright)
        splitter3.addWidget(topright)

        splitter1 = QtGui.QSplitter(QtCore.Qt.Horizontal)
        splitter1.addWidget(topleft)
        splitter1.addWidget(splitter3)
        splitter2 = QtGui.QSplitter(QtCore.Qt.Vertical)
        splitter2.addWidget(splitter1)
        splitter2.addWidget(bottom)
        hbox.addWidget(splitter2)
        self.setLayout(hbox)
        QtGui.QApplication.setStyle(QtGui.QStyleFactory.create('Cleanlooks'))
        self.setGeometry(300, 300, 300, 200)
        self.setWindowTitle('Dose curves analysis')
        self.show()

    def noiseAmplitude(self, text):
       # qti.app.resultsLog().append(text)
        self.amplitude = float(text)

    def onActivated(self, text):
        self.mode = text
    def pushed(self):
        self.m = self.sender()
        if self.m.text() == 'Generate':
            if not bp.objectMessanger('initParams', 'Table'):
                return
            elif not bp.objectMessanger('SimParam', 'Table'):
                return
            if self.mode == 'N':
                bp.generPmLt(self.amplitude, 'N')
            elif self.mode == 'NU':
                bp.generPmLt(self.amplitude, 'NU')
        elif self.m.text() == 'Logarithmic differentiate':
            if self.mode == 'N':
                bp.logDerivat('N')
            elif self.mode == 'NU':
                bp.logDerivat('NU')
        elif self.m.text() == 'Differentiation parameters':
            if not bp.objectMessanger('initParams', 'Table'):
                return
            elif not bp.objectMessanger('SimParam', 'Table'):
                return
            if self.mode == 'N':
                bp.derivatParams('N')
            elif self.mode == 'NU':
                bp.derivatParams('NU')

win = Window()
win.move(600, 400)

win.show()
