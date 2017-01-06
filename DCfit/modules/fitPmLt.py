"""
    v0.6.1  
    This is a package  fitPmLt  which holds function used in protein dose
    curves fitting programs.

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

## Function for fitting Pm(Lt) curves ([N]ative protein binding to ligand model).
#  @param fitWizard Qti Table Object : template table
#  @param currentTableCols Qti Table Object : selected columns of the current table
#  @return Qti Table Object : results table(data curves)
#  @return Qti Table Object : results params table(parameters values)
#  @return Qti Graph Object : results graph(plot curves)
def fitPmLt(tablename, model='N'):
    try:
        ## initializing parameters
        startTime = time.time()
        tempnum = '_1'
        useGlobal = True # logic parameter for using globality in the module
        t = qti.app.table(tablename)
        if  not isinstance(t, qti.Table):
            t = qti.app.currentTable()
            if not isinstance(t, qti.Table):
                msgBox = QtGui.QMessageBox()
                msgBox.setText("Selected window is not of type 'Table' y")
                ret = msgBox.exec_()
        wizard = qti.app.table("fitWizard")
        if model == 'N':    
            fit_vars = ["DuG", "DuV", "DuBeta", "Kbn", "DbnV","DbnBeta", "Mt"]
            funcname = bp.LtUpper
            funcnameR = bp.LtUpper
            def objective(params, x, data, funcname, other_params):
                ndata, nx = data.shape
                resid = 0.0*data[:]
                for i in range(ndata):

                    for k in range(nx):
                        hll = bp.brent(funcname, x[k], params, i, other_params)
                        resid[i, k] = data[i, k] - hll
                        
                return resid.flatten()
        else:
            fit_vars = ["DuG", "DuV", "DuBeta", "Kbn","DbnV","DbnBeta", "Kbu","DbuV","DbuBeta","Mt"]
            funcname = bp.PmLtFull
            funcnameR = bp.LtFull
            def objective(params, x, data, funcname, other_params):
                ndata, nx = data.shape
                resid = 0.0*data[:]
                for i in range(ndata):
                    for k in range(nx):
                        resid[i, k] = data[i, k] - funcname(x[k], params, i, other_params)
                return resid.flatten()
        
        
        answer = bp.columnToArray(tablename)
        x = answer[0]
        yy =  answer[1]
        colindexes = answer[2]
        colnames = answer[3]
        if yy.shape[0] == 1:
            useGlobal = False
        if yy.shape[0]>14:
            rows = yy.shape[0]
        else:
            rows = 15
        fitWiz = bp.fitWizardParameters(wizard, useGlobal, fit_vars)
        prams = fitWiz[0]
        other_params = fitWiz[1]
        globalD = fitWiz[2]
        if model == 'N':
            del other_params['Kbu']
            del other_params['DbuV']
            del other_params['DbuBeta']
        if len(globalD) == 0:
            useGlobal = False
       
        if useGlobal:
            y = yy
            fit_params = Parameters()
            ## initiating fit parameters
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
            
            
            
            ## FIT
            fitStartTime = time.time()
            minimizer = minimize(objective, fit_params, args=(x, y, funcname, other_params))     
         
            rezStartTime = time.time()
            fitRezTName = bp.fitResult(x, funcnameR, minimizer.params, other_params, t, colindexes, brent = True, fitTableName = "TfitRes"+ model +'_1')

            a = bp.Rsquare(qti.app.table(fitRezTName))
            p_values = bp.runs_test(qti.app.table(fitRezTName))
            chi_squares = bp.chi_squares(qti.app.table(fitRezTName),  minimizer.nvarys)
            chi_sqr = chi_squares[0]
            red_chi = chi_squares[1]

            #for yi, i in enumerate(a):
                    #fit_params.add( 'Rsquare_%i' % (yi+1)    , value = i)
            dict = {}
            sYcols = t.selectedYColumns()
            yNames = []
            for i in sYcols:
                index = t.colIndex(i)
                yNames.append (t.colName(index + 1))
            dict['yNames'] = yNames
            dict['intable'] = t.objectName()
            dict['outtable'] = bp.existTableName('TfitParam'+ model +'_1')
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
            #minimizer.params, fit_value, fit_vars, y, fitRezName=fitRezTName, useGlobal = useGlobal, globalD = globalD, tName = 'TfitParam'+ model +'_1' )
              

            pdc ={ 'y_name':"<i>P<sub>m</sub></i>, MPa",
                    'x_name':"<i>L<sub>t</sub></i>, M",
                    'change_x': True,
                    'x_min' : x[0],
                    'x_max' : x[-1],
                    'x_scale' : 1,
                    'x_num_format' : 5,
                    'DotsLine' : 1,
                    'graphName' : 'Gfit' + model +'_1' }
            bp.Geditor(dict=pdc, fromTable=True, table=qti.app.table(fitRezTName))

            g = qti.app.currentGraph()
            fitwiz = qti.app.table('fitWizard')
            fitwiz.showMaximized()
            fitwiz.showNormal()
            g.showMaximized()
            g.showNormal()
            ##Write rezults to log
            log = ''
            log += '\n'
            log += "Date: " + time.strftime("%Y-%m-%d %H:%M:%S") + '\n'
            log += "Model: " + model + '\n'
            log += '[[Input data]]'
            log += "  Input table name: " + str(t.objectName()) + '\n'
            log += "  X data name: " + str(colnames[0][0])+'\n'
            for i, val in enumerate(colnames[1]):
                log += "  Y" +str(i + 1) +" data name: " + val +'\n'
            log += '[[Output data]]'
            log += "  Output table name: " + str(outputTable.objectName()) + '\n'
            log += "  Fit statistics table name: " + fitRezTName + '\n'
            log += "  Fit result graph name: " + g.objectName() + '\n'
            log += '[[Fit individual statistics]]'
            for yi, i in enumerate(a):
                log +=  '  R^2_%i' % (yi+1) + ": " + str(i) + ' \n '
            for yi, i in enumerate(p_values):
                log +=  '  P_%i' % (yi+1) + ": " + str(i) + ' \n '
            for yi, i in enumerate(chi_sqr):
                log +=  '  chi-square_%i' % (yi+1) + ": " + str(i) + ' \n '
            for yi, i in enumerate(red_chi):
                log +=  '  reduced chi-square_%i' % (yi+1) + ": " + str(i) + ' \n '
            log += fit_report(minimizer) + '\n'
            endTime = time.time()
            log += '[[Execution Time]]' + '\n'
            log += '  Total: ' + str("%.4g"%( endTime - startTime) )+'s \n ' 
            log += '  Fit time: ' + str(  "%.4g"%(rezStartTime - fitStartTime))+'s \n' 
            log += '  Data initiating time: ' + str("%.4g"%(fitStartTime - startTime))+'s \n' 
            log += '  Results displaying time2: ' + str("%.4g"%(endTime - rezStartTime))+'s \n' 
            log += '------------------------------------------------------------------------ \n'
            qti.app.resultsLog().append(log)
                            
        else:
            for iter in range(yy.shape[0]):
                startTime = time.time()
                fit_params = Parameters()
                y = np.array([yy[iter]])
                currselcols = [colindexes[0] , [colindexes[1][iter]]]
                yNames = [t.colName(colindexes[1][iter] + 1)]
                ## initiating fit parameters
                for name in prams:
                    for yi, yj in enumerate(y):
                        fit_params.add( name + '_%i' % (yi + 1)    , value = prams[name].value, min = prams[name].min, max = prams[name].max, vary = prams[name].vary)
                       
               
                fit_value = []
                for u  in fit_vars:
                    fit_value.append(u)
                    fit_value.append(fit_params[u+tempnum].value)
                    fit_value.append(fit_params[u+tempnum].min)
                    fit_value.append(fit_params[u+tempnum].max)
                    fit_value.append(fit_params[u+tempnum].vary)
                
                
                
                ## FIT
                fitStartTime = time.time()
                minimizer = minimize(objective, fit_params, args=(x, y, funcname, other_params))     
             
                rezStartTime = time.time()
                fitRezTName = bp.fitResult(x, funcnameR, minimizer.params, other_params, t, currselcols, brent = True, fitTableName = "TfitRes"+ model +'_1')

                a = bp.Rsquare(qti.app.table(fitRezTName))
                p_values = bp.runs_test(qti.app.table(fitRezTName))
                chi_squares = bp.chi_squares(qti.app.table(fitRezTName),  minimizer.nvarys)
                chi_sqr = chi_squares[0]
                red_chi = chi_squares[1]


                dict = {}
                dict['yNames'] = yNames
                dict['intable'] = t.objectName()
                dict['outtable'] = bp.existTableName('TfitParam'+ model +'_1')
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
                #minimizer.params, fit_value, fit_vars, y, fitRezName=fitRezTName, useGlobal = useGlobal, globalD = globalD, tName = 'TfitParam'+ model +'_1' )
                  
                pdc ={ 'y_name':"<i>P<sub>m</sub></i>, MPa",
                        'x_name':"<i>L<sub>t</sub></i>, M",
                        'change_x': True,
                        'x_min' : x[0],
                        'x_max' : x[-1],
                        'x_scale' : 1,
                        'x_num_format' : 5,
                        'DotsLine' : 1,
                        'graphName' : 'Gfit' + model +'_1' }
                bp.Geditor(dict=pdc, fromTable=True, table=qti.app.table(fitRezTName))

                g = qti.app.currentGraph()
                fitwiz = qti.app.table('fitWizard')
                fitwiz.showMaximized()
                fitwiz.showNormal()
                g.showMaximized()
                g.showNormal()
                
                
                ##Write rezults to log
                
                qti.app.updateLog('\n')
                qti.app.updateLog("Date: " +time.strftime("%Y-%m-%d %H:%M:%S") + '\n')
                qti.app.updateLog("Model: " + model + '\n')
                qti.app.updateLog('[[Input data]]')
                qti.app.updateLog("  Input table name: " + str(t.objectName()) + '\n')
                qti.app.updateLog("  X data name: " + str(colnames[0][0])+'\n')
                for i, val in enumerate(colnames[1]):
                    qti.app.updateLog("  Y" +str(i + 1) +" data name: " + val +'\n')
                qti.app.updateLog('[[Output data]]') 
                qti.app.updateLog("  Output table name: " + str(outputTable.objectName()) + '\n')
                qti.app.updateLog("  Fit result table name: " + fitRezTName + '\n')
                qti.app.updateLog("  Fit result graph name: " + g.objectName() + '\n')
                qti.app.updateLog('[[Fit individual statistics]]') 
                for yi, i in enumerate(a):
                    qti.app.updateLog( 'R^2_%i' % (yi+1) + ": " + str(i) + ' \n ')
                for yi, i in enumerate(p_values):
                    qti.app.updateLog( 'P_%i' % (yi+1) + ": " + str(i) + ' \n ')
                for yi, i in enumerate(chi_sqr):
                    qti.app.updateLog( 'chi-square_%i' % (yi+1) + ": " + str(i) + ' \n ')
                for yi, i in enumerate(red_chi):
                    qti.app.updateLog( 'reduced chi-square_%i' % (yi+1) + ": " + str(i) + ' \n ')
                qti.app.updateLog(fit_report(minimizer) + '\n')
                endTime = time.time()
                qti.app.updateLog('[[Execution Time]]' + '\n') 
                qti.app.updateLog('  Total: ' + str("%.4g"%( endTime - startTime) )+'s \n ' +
                            '  Fit time: ' + str(  "%.4g"%(rezStartTime - fitStartTime))+'s \n' +
                            '  Data initiating time: ' + str("%.4g"%(fitStartTime - startTime))+'s \n' +
                            '  Results displaying time3: ' + str("%.4g"%(endTime - rezStartTime))+'s \n' +
                            '------------------------------------------------------------------------ \n'
                            )
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

