"""
    v0.6.1  
    This is a package  fitPmTm  which holds function used in protein melting 
    pressure-temperature phase diagram fitting program.

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

"""
import numpy as np
from math import exp, log, sqrt
import psaFit as bp

import qti
from PyQt4 import QtCore, QtGui
import traceback
from lmfit import minimize, Parameters, Parameter, report_fit, fit_report
import time
import sys



## Function for fitting Pm(Tm) curves (pressure-temperature phase diagram).
#  @param fitWizard Qti Table Object : template table
#  @param currentTableCols Qti Table Object : selected columns of the current table
#  @return Qti Table Object : results table(data curves)
#  @return Qti Table Object : results params table(parameters values)
#  @return Qti Graph Object : results graph(plot curves)
def newfitPmTm(tablename = ''):
    try:
        ## initializing parameters
        startTime = time.time()
        tempnum = '_1'
        fit_vars = ["DuG", "DuV", "DuBeta", "DuAlpha", "DuCp", "DuH"]
        t = qti.app.table(tablename)
        if  not isinstance(t, qti.Table):
            t = qti.app.currentTable()
            if not isinstance(t, qti.Table):
                msgBox = QtGui.QMessageBox()
                msgBox.setText(str(e))
                ret = msgBox.exec_()
        wizard = qti.app.table("fitWizard")
        useGlobal = False # logic parameter for using globality in the module
        answer = bp.columnToArray(tablename)
        x = answer[0]
        yy =  answer[1]
        colindexes = answer[2]
        
        if yy.shape[0] == 1:
            useGlobal = False
        if yy.shape[0]>14:
            rows = yy.shape[0]
        else:
            rows = 15
        fitWiz = bp.fitWizardParameters(wizard, useGlobal, fit_vars)
        prams = fitWiz[0]
        other_params = fitWiz[1]
        if useGlobal:
            globalD = fitWiz[2]
        else:
            globalD = {}
        
       
            #if fit_params['DuBeta'+ '_%i' % (yi + 1)].vary == True:
               # fit_params.add('eplips_%i' % (yi + 1), value = 100, min = 1e-09, vary = True)
               # fit_params['DuBeta'+ '_%i' % (yi + 1)].expr = '-(eplips_%i  + DuAlpha_%i **2)*' % (yi + 1, yi + 1)+ str(other_params['T']) +'/DuCp_%i' % ( yi + 1)
                
            #elif fit_params['DuAlpha'+ '_%i' % (yi + 1)].vary == True:
                #fit_params.add('eplips_%i' % (yi + 1), value = 2, min = 1e-09, vary = True)
                #fit_params['DuAlpha'+ '_%i' % (yi + 1)].expr = 'sqrt(-eplips_%i  - DuBeta_%i/' % (yi + 1, yi + 1)+ str(other_params['T']) +'*DuCp_%i)' % ( yi + 1)
            #fit_params.add('deltaDuCp_%i' % (yi+1), value = -0.3218791946, max=(0-1e-09), vary=False)

            #fit_params.add('DuCp_%i' % (yi+1), expr='(deltaDuCp_%i -DuAlpha_%i **2)*T_%i /DuBeta_%i' % (yi + 1, yi + 1, yi + 1, yi + 1))
            #bp.ellipsoidEnthalpy(fit_params, x, yi)
        def objective(params, x, data, funcname, other_params):
                ndata, nx = data.shape
                resid = 0.0*data[:]
                for i in range(ndata):
                    for k in range(nx):
                        dat = data[i, k]
                        model = funcname(params, x[k], dat, other_params= other_params)
                        resid[i, k] = model
                return resid.flatten()
                    
        ## initiating fit parameters
        for iter in range(yy.shape[0]):
            fit_params = Parameters()
            y = np.array([yy[iter]])
            currselcols = [colindexes[0] , [colindexes[1][iter]]]
            yNames = [t.colName(colindexes[1][iter] + 1)]
            for name in prams:
                for yi, yj in enumerate(y):
                    fit_params.add( name + '_%i' % (yi + 1)    , value = prams[name].value, min = prams[name].min, max = prams[name].max, vary = prams[name].vary)
                       
            ## ititiating and assigning global fit parameters
            if useGlobal:
                for name in globalD:
                    if globalD[name]:
                        fit_params.add( name , value = prams[name].value, min = prams[name].min, max = prams[name].max, vary = prams[name].vary)
                        for yi, yj in enumerate(y):
                            fit_params[name + '_%i' % (yi + 1)].expr = name
                            
            fit_value = []
            for u  in fit_vars:
                fit_value.append(u)
                fit_value.append(fit_params[u+tempnum].value)
                fit_value.append(fit_params[u+tempnum].min)
                fit_value.append(fit_params[u+tempnum].max)
                fit_value.append(fit_params[u+tempnum].vary)
            funcname = bp.Gibsenergy
            
            
            
            ## FIT
            fitStartTime = time.time()
            minimizer = minimize(objective, fit_params, args=(x, y, funcname, other_params))

            rezStartTime = time.time()
            funcnameR = bp.ellipsoidPmTm
            fitRezTName = bp.fitResultElip(x, funcnameR, minimizer.params,  other_params, t, currselcols, fitTableName = "TfitResElip_1")
            a = bp.Rsquare(qti.app.table(fitRezTName))
            p_values = bp.runs_test(qti.app.table(fitRezTName))
            chi_squares = bp.chi_squares(qti.app.table(fitRezTName), \
                minimizer.nvarys)
            chi_sqr = chi_squares[0]
            red_chi = chi_squares[1]

            dict = {}
            dict['yNames'] = yNames
            dict['intable'] = t.objectName()
            dict['outtable'] = bp.existTableName('TfitParElip_1')
            dict['fittable'] = qti.app.table(fitRezTName).objectName()
            dict['fitvars'] =  fit_vars
            dict['globalD'] = globalD
            dict['initialvals'] = fit_value
            dict['otherparams'] = other_params
            dict['p_values'] = p_values
            dict['r2'] = a
            dict['chi_sqr'] = chi_sqr
            dict['red_chi'] = red_chi
            outputTable = bp.outputTable(minimizer, dict)
            pdc ={ 'y_name':"<i>P<sub>m</sub></i>, MPa",
                    'x_name':"<i>T<sub>m</sub></i>, K",
                    'change_x': True,
                    'x_scale' : 0,
                    'x_num_format' : 4,
                    'DotsLine' : 1,
                    'x_min' : min(x),
                    'x_max' : max(x),
                    'graphName' : 'GfitElip_1' }
            bp.Geditor(dict=pdc, fromTable=True, table=qti.app.table(fitRezTName))

            g = qti.app.currentGraph()
            fitwiz = qti.app.table('fitWizard')
            fitwiz.showMaximized()
            fitwiz.showNormal()
            g.showMaximized()
            g.showNormal()

            ##Write rezults to log
            qti.app.updateLog(' \n ')
            qti.app.updateLog("Date: " +time.strftime("%Y-%m-%d %H:%M:%S") + ' \n ')

            qti.app.updateLog("Input table name: " + str(t.objectName()) + ' \n ')
            qti.app.updateLog("Output table name: " + str(outputTable.objectName()) + ' \n ')
            qti.app.updateLog("Fit result table name: " + fitRezTName + ' \n ')
            qti.app.updateLog("Fit result graph name: " + g.objectName() + ' \n ')
            qti.app.updateLog("function: " + funcname.__name__ + ' \n ')
            for yi, i in enumerate(a):
                qti.app.updateLog( 'Rsquare_%i' % (yi+1) + ": " + str(i) + ' \n ')
            for yi, i in enumerate(p_values):
                qti.app.updateLog( 'P_%i' % (yi+1) + ": " + str(i) + ' \n ')
            for yi, i in enumerate(chi_sqr):
                qti.app.updateLog( 'Chisqr_%i' % (yi+1) + ": " + str(i) + ' \n ')
            for yi, i in enumerate(red_chi):
                qti.app.updateLog( 'redChi_%i' % (yi+1) + ": " + str(i) + ' \n ')
            qti.app.updateLog(fit_report(minimizer) + ' \n ')
            endTime = time.time()
            #msgBox = QtGui.QMessageBox()
            #msgBox.setText("Done")
            qti.app.updateLog('Time: ' + str("%.4g"%( endTime - startTime) )+'s \n ' +
                                                            'Fit time: ' + str(  "%.4g"%(rezStartTime - fitStartTime))+'s \n ' +
                                                            'Data initiating time: ' + str("%.4g"%(fitStartTime - startTime))+'s \n ' +
                                                            'Results displaying time: ' + str("%.4g"%(endTime - rezStartTime))+'s \n ' +
                                                            '--------------------------------------------------------------------------------------------------------------------  \n'
                                                            )
            #ret = msgBox.exec_()

    except Exception, e:
        msgBox = QtGui.QMessageBox()
        msgBox.setText("Error!!!")
        if str(e) == "'NoneType' object has no attribute 'selectedColumns'" or str(e) == "integer division or modulo by zero":
            msgBox.setInformativeText("You must to select x and y columns from the table")
        else:
            msgBox.setInformativeText(str(e))
        msgBox.setDetailedText(str(traceback.format_exc()))
        ret = msgBox.exec_()
        qti.app.resultsLog().append(str(traceback.format_exc()))
