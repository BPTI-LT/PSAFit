"""
    v0.6.1  
    This is a package  basic  which holds function used in all other packages.

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
import psaFit as bp
import random as rnd
import time
import qti
from lmfit import Parameters
import logging
from PyQt4 import QtCore, QtGui
from math import exp, log, sqrt
import scipy.stats as st
import collections
from math import log10, floor

'''
Code below is used for QtiPlot tables and graphs editing
'''


def round_sig(x, sig=3):
    if abs(x) < 0.0000000000001:
        return x
    if abs(x) < 1:
        return round(x, sig-int(floor(log10(abs(x))))-1)
    else:
        return round(x, sig)

## Function for getting QtiPlot tables column names in one string array.
#  @param inputTable Qti Object table : table for column names reading
#  @return String Array : list with column names
def colNames(inputTable):
    answer = []
    for i in range(inputTable.numCols()):
        answer.append(inputTable.colName(i+1))
    return answer

## Function for converting numpy array data points to Qtiplot table data points.
#  @param output_table Qti Tabele Object : QtiPlot table for data
#  @param x Numpy 1D array : absicissae data points
#  @param data Numpy 2D array : ordinates data points
def arrayToTable(output_table,x,data):
    for i in range(0,len(x)):
        output_table.setCellData(1, i+1, x[i])
    for k in range(len(data)):
        for i in range(len(data[k])):
            output_table.setCellData(k+2, i+1, data[k][i])


## Function for converting QtiPlot SimParam table to dictionary
#  @param input_table=qti.app.table("SimParam")  Qti Table Object : QtiPlot SimParam table
#  @return Dict : dictionary with keys - column names, and values Numpy 1D arrays with all numerical values in those columns
def simParamDict(inputTable = qti.app.table("SimParam")):
    answer = {}
    for m in range(inputTable.numCols()):
        key = str(inputTable.colLabel(m)) # note that m = 0 (like in arrays), but colIndex = 1
        value = []
        for k in range(inputTable.numRows()):
            if not emptyCell(inputTable,m + 1,k + 1):
                value.append(inputTable.cell(m + 1,k + 1))
        value = np.array(value)
        answer.update({key:value})
    return answer



## Function for converting QtiPlot current table selected columns data to numpy array data points (x y data points).
#  @return 2D Numpy Array : x (return[0]) and y (return[1]) values
def columnToArray(tablename = ''):
    x = []
    dat = []
    xcolindexes = []
    ycolindexes = []
    xcolnames = []
    ycolnames = []
    t = qti.app.table(tablename)
    if  not isinstance(t, qti.Table):
        t = qti.app.currentTable()
        if not isinstance(t, qti.Table):
            msgBox = QtGui.QMessageBox()
            msgBox.setText('Select data table please')
            ret = msgBox.exec_()
            return
    colvalues = []
    xx = []
    yy = []
    cols = t.selectedColumns()
    if len(cols) == 0:
        msgBox = QtGui.QMessageBox()
        msgBox.setText('There are no selected columns in current table')
        ret = msgBox.exec_()
        return
    for k in cols:
        index = t.colIndex(k)
        des = t.colPlotDesignation(index)
        if des == 1:
            xx.append(k)
            xcolindexes.append(index)
            xcolnames.append(t.colName(index+1))
        elif des == 2:
            yy.append(k)
            ycolindexes.append(index)
            ycolnames.append(t.colName(index+1))
    xy = xx
    if len(xx) == 0:
        msgBox = QtGui.QMessageBox()
        msgBox.setText('There is no column of type X in selected columns.')
        ret = msgBox.exec_()
        return
    xy.extend(yy)
    if len(xy) == 0:
        msgBox = QtGui.QMessageBox()
        msgBox.setText('Data in table does not have x column. Please select correct table.')
        ret = msgBox.exec_()
        return
    for j in range(t.numRows()):
        empty = 0
        for k in xy:
            if emptyCell(t, k, j+1):
                empty = empty + 1
        if empty%len(xy) != 0:
            qti.app.resultsLog().append('Row '+str(j+1) + ' is filled incorrectly. The results might be damaged')
    answer = 1
    for m in xy:
        index = t.colIndex(m)
        data_y = []
        for k in range(1, t.numRows()+1):
            if not emptyCell(t,m,k):
                if t.colPlotDesignation(index) == 1:
                    x.append(t.cell(m,k))
                elif t.colPlotDesignation(index) == 2:
                    data_y.append(t.cell(m,k))
        if len(data_y) >1:
            dat.append(data_y)
    x =     np.array(x)
    data = np.array(dat)
    answer = [x, data, [xcolindexes, ycolindexes], [xcolnames, ycolnames]]
    return answer

## Function for converting QtiPlot table data points to numpy array data points.
#  @param input_table Qti Table Object : input table
#  @return Numpy 2D Array : x (return[0]) and y (return[1]) values
def tableToArray(input_table):
    x = []
    dat = []

    for m in range(1, input_table.numCols()+1):
        data_y = []
        for k in range(1, input_table.numRows()+1):
            if not emptyCell(input_table,m,k):
                if m == 1:
                    x.append(input_table.cell(m,k))
                else:
                    data_y.append(input_table.cell(m,k))
            if emptyCell(input_table,m,k): #it prints which cell is empty
                message = 'there is missing data point in cell('+str(m)+','+str(k)+')'
                print message
        if len(data_y) >1:
            dat.append(data_y)
    x =     np.array(x)
    data = np.array(dat)
    answer = [x, data]
    return answer

## Function for checking if particular QtiPlot table cell is empty.
#  @param table Qti Table Object : input table
#  @param column_number Int : column number
#  @param row_number  Int : row number
def emptyCell( table, column_number, row_number):
    should_break = False
    m=int(table.cell(column_number,row_number))
    try:
        n=unicode(m)
    except ValueError:
        should_break=True
    return should_break


## Function for generating PmLt results QtiPlot table.
#  @param x 1D Numpy Array : Lt values
#  @param x 2D Numpy Array : Pm values
#  @param param="none" String : name of varied parameter
#  @param outTableName="output_table_1" String : name of output table
#  @param SimParam=qti.app.table("SimParam") Qti Table Object : table with simulated parameters
#  @param InitParams=qti.app.table("initParams") Qti Table Object : table with initial parameters
#  @return Qti Table Object : resul table of generated Pm(Lt) curves
def generatorResultTable(x, data, dict, param = "none", outTableName = "output_table_1" , SimParam = qti.app.table("SimParam"), InitParams = qti.app.table("initParams") ):

    outTableName = existTableName(outTableName)
    rows = len(x)
    if rows <= 29:
        rows = 29
    cols = len(data) + 5
    t1 = qti.app.newTable(outTableName, rows, cols)

    if param != "none":
        sims = simParamDict(SimParam)
        sim = sims[param]
        t1.setColName(len(sim) + 3, "initPar", False)
        t1.setColName(len(sim) + 4, "value", False)
        for yi, y in enumerate(sim):
            a = "%.4g"% y
            a = param + "=" + a
            t1.setColName(yi + 2, a, False)
    t1.setColName(1, "Lt", False)
    arrayToTable(t1,x,data)
    t1.setColTextFormat(len(data) + 3)
    kk = []
    vv = []
    for key, value in dict.iteritems():
        kk.append(key)
        vv.append(value)

    for i in range(len(kk)):
        t1.setCellData(len(data) + 3, i + 1, kk[i])
        t1.setCellData(len(data) + 4, i + 1, vv[i])


    mm = [param,
            time.strftime("%Y-%m-%d %H:%M:%S") ]
    nn = ['varied parameter',
            'date/time']
    for i in range(len(nn)):
        t1.setColTextFormat(len(data) + 2)
        t1.setText(len(data)+2, len(kk) + i + 2, nn[i] )
        t1.setText(len(data)+3, len(kk) + i + 2, str(mm[i]))
    return t1


## Function for parameters initializing from table to dictionary.
#  @param init_table Qti Table Object :  parameters initializing table
#  @return Dictionary : dictionary for parameters keys and values.
def initParamsDict(init_table):
    dict = collections.OrderedDict()
    for i in range(1, init_table.numRows()):
        if not emptyCell(init_table, 1, i):
            if not emptyCell(init_table, 2, i):
                dict.update({init_table.cellData(1,i):init_table.cell(2,i)})
    return dict

## Function for parameters initializing from table to Parameters Class.
#  @param fit_wiz_table Qtiplot table :  parameters initializing table
#  @param useGlobal=False Boolean : parameter for reading globality value
#  @return Parameters Class : parameters for fitting
#  @return name None|Str : fitting parameter name
#  @return value Float : fitting parameter value
#  @return min Float : fitting parameter minimum value
#  @return max Float : fitting parameter mazimum value
#  @return vary Boolean : fitting parameter vary on/off value
#  @return Dictionary : keys - parameters names, values value of globality(True or False); additional return case
#  @return LmFit Parameters Object : name, value, min, max, vary values of all parameters; default return case then useGlobal=False
#  @return List : default return case placed in 0 place, additional return case place in 1 place; then useGlobal=True
def fitWizardParameters(fit_wiz_table, useGlobal = False, fitvars = {}):
    fit_params = Parameters()
    dict = collections.OrderedDict() # local dictionary which is used for creating FitParams - lmfit Parameters object
    other_params = collections.OrderedDict() # local dictionary for non fit params
    globalD = {} # local dictionary which is used for storing information about parameters globality
    douseGlobal = False
    for i in range(1, fit_wiz_table.numRows()):
        if not emptyCell(fit_wiz_table, 1, i):
            name = fit_wiz_table.cellData(1,i)
            if name in fitvars:
                dict = {'value':fit_wiz_table.cell(2,i),
                                'min':fit_wiz_table.cell(4,i),
                                'max':fit_wiz_table.cell(5,i),
                                'vary':bool(fit_wiz_table.cell(3,i)),}
                if emptyCell(fit_wiz_table, 2, i):
                    dict.update({'value':None})
                if emptyCell(fit_wiz_table, 4, i):
                    dict.update({'min':float("-inf")})
                if emptyCell(fit_wiz_table, 5, i):
                    dict.update({'max':float("inf")})
                if emptyCell(fit_wiz_table, 3, i):
                    dict.update({'vary':False})
                fit_params.add( name, value = dict['value'], min = dict['min'], max = dict['max'], vary = dict['vary'])
                if useGlobal:
                    if emptyCell(fit_wiz_table, 6, i) or not bool(fit_wiz_table.cell(3,i)):
                        globalD.update({name:False})
                    else:
                        globalD.update({name: bool(fit_wiz_table.cell(6,i))})
                        douseGlobal = True
            else:
                other_params[name] = fit_wiz_table.cell(2,i)
    if not douseGlobal:
        globalD = {}
    return [fit_params, other_params, globalD]
  


## Function for generating original table names.
#  @param table_name Str : name for checking
#  @return Str : original table name
def existTableName(table_name):
    i = 1
    while i < 1000:
        to = -1
        if i > 10:
            to = -2
        if i > 100:
            to = -3
        table_name = table_name[:to] + str(i)
        try:
            x = qti.app.table(table_name)
            x.cell(1,1)
        except:
            break
        i=i+1
    return table_name

## Function for creating fit result table.
#  @param x Array : all x data points
#  @param funcname Function : mathematical function for mathematical formula ex. Gaussian
#  @param params LmFit Parameters Object : parameters values
#  @param experm_table Str : experimental data table name
#  @param selcolindexes Array : 2D array of indexes of selected columns
#  @param brent=False Boolean : logic parameter for using Brent's method
#  @param fitTableName="fit_result_1" Str : output table name
#  @return Str : name of fit result table
def fitResult(x, funcname, params, other_params, experm_table, selcolindexes, brent = False, fitTableName = "fit_result_1"):
    fitTableName = existTableName(fitTableName)
    cols = len(selcolindexes[1])
    xcol = selcolindexes[0][0]
    t1 = qti.app.newTable(fitTableName, len(x), 2 * cols + 1)
    t2 = experm_table
    sc = selcolindexes[1]

    for i in range(len(x)):
        t1.setCellData(1, i+1, x[i])
        for k in range(cols):
            if not brent:
                try:
                    y = funcname(params, x[i],data = [], i = k, other_params = other_params)
                except:
                    y = funcname(params, x[i], other_params)
                t1.setCellData((k+1) * 2 + 1, i+1, y)
            else:
                if funcname == bp.LtFull:
                    Pm_range_calc = bp.findPmRange(bp.Pm_low, bp.LtFull, params, k, other_params)
                    xa = Pm_range_calc['low']
                    xb = Pm_range_calc['high']
                    other_params.update({'xa':xa})
                    other_params.update({'xb':xb})
                y = bp.brent(funcname, x[i], params, k, other_params)
                t1.setCellData((k+1) * 2 + 1, i+1, y)
    for i in range(len(x)):
        for yk, k in enumerate(sc):
            y = t2.cell(k + 1, i+1)
            t1.setCellData((yk + 1) * 2 , i+1, y)
    ##  setting column values
    t1.setColName(1, t2.colName(xcol + 1))
     
    for yi, i in enumerate(sc):
        value = t2.colName(i + 1)
        t1.setColName((yi + 1 ) * 2 + 1, 'fit%i' % (yi+1))
        t1.setColName(yi + 2, value)
    
    
    
    for yi, i in enumerate(sc):
        value = t2.colName(i + 1)
        if not brent:
            t1.setColName(yi + 1, value)
        else:
            t1.setColName((yi  + 1) * 2 + 1, 'fit%i' % (yi + 1))
            t1.setColName( yi + 2, value)
    return fitTableName

## Function for QtiPlot graphs editing.
#  @param fromTable=False Boolean : if True - editing current table, if False - editing current Graph
#  @param dict Dictionary : graphs editing arguments
#  @param dict['graphName']="graph_1" String : graph name
#  @param dict['title']=False String : graph title
#  @param dict['x_name']="<i>x</i>, UNIT" String : x axis name
#  @param dict['y_name']="<i>y</i>, UNIT" String : y axis name
#  @param dict['Line_on']=True Boolean : logic parameter for turning on/of curves lines
#  @param dict['Dots_on']=True Boolean : logic parameter for turning on/of curves data points dots
#  @param dict['DotsLine']=0 Int : parameter for line/dots variation. 0 - False, 1 - dots line variantion, 2 line dots variation
#  @param dict['changeVerticalAxis']=0 Int : parameter for curves axis variation. 0 - all curves on left axis, 1 - left rigth..., 2 - right left, 3 - all right
#  @param dict['Grid_on']=False Boolean : logic parameter for turning on/of grid
#  @param dict['change_x']=False Boolean : logic parameter for turning on/of bottom axis formating (x_min, x_max, x_scale, x_num_format)
#  @param dict['x_min']=0.0 Float : bottom axis start
#  @param dict['x_max']=1000.0 Float : bottom axis end
#  @param dict['x_scale']=0 Int : bottom axis scale- [0, 1, 2, 3, 4, 5, 6] = [Linear, Log10, ln, Log2, Reciprocal(1/T), Probability, Logit]
#  @param dict['x_num_format']=1 Int : bottom axis numerical format- [0, 1, 2, 3, 4, 5] = [automatic, 10000.0, 1e4, 1x10^4, 10k, 1 #dot 10^4]
#  @param dict['change_y']=False Boolean : logic parameter for turning on/of left axis formating (y_min, y_max, y_scale, y_num_format)
#  @param dict['y_min']=0.0 Float : left axis start
#  @param dict['y_max']=1000.0 Float : left axis end
#  @param dict['y_scale']=0 Int : left axis scale- [0, 1, 2, 3, 4, 5, 6] = [Linear, Log10, ln, Log2, Reciprocal(1/T), Probability, Logit]
#  @param dict['y_num_format']=1 Int : left axis numerical format- [0, 1, 2, 3, 4, 5] = [automatic, 10000.0, 1e4, 1x10^4, 10k, 1 #dot 10^4]
#  @param dict['change_canvas']=False Boolean : logic parameter for turning on/of canvas size changing
#  @param dict['width']=400 Int : graph width by points
#  @param dict['height']=400 Int : graph height by points
#  @param dict['twoYAxis']=400 Boolean : parameter for enableing right axis(y2)
#  @return Qti Graph Object : formated graph
def Geditor(dict={'editor':True}, fromTable=False, table=qti.app.currentTable()):
    dict.setdefault('graphName', 'graph_1')
    dict.setdefault('title', False)
    dict.setdefault('x_name', "<i>x</i>, UNIT")
    dict.setdefault('y_name', "<i>y</i>, UNIT")
    dict.setdefault('y2_name', "")
    dict.setdefault('Line_on', True)
    dict.setdefault('Dots_on', True)
    dict.setdefault('DotsLine', 0) # 0 - False, 1 - dots line variantion, 2 line dots variation
    dict.setdefault('changeVerticalAxis', 0) # 0 - all curves on left axis, 1 - left rigth..., 2 - right left, 3 - all right
    dict.setdefault('Grid_on', False)
    dict.setdefault('change_x', False)
    dict.setdefault('x_min', 0.0)
    dict.setdefault('x_max', 1000.0)
    dict.setdefault('x_scale', 0) # [0, 1, 2, 3, 4, 5, 6] = [Linear, Log10, ln, Log2, Reciprocal(1/T), Probability, Logit]
    dict.setdefault('x_num_format', 1)  #[0, 1, 2, 3, 4, 5] = [automatic, 10000.0, 1e4, 1x10^4, 10k, 1 #dot 10^4]
    dict.setdefault('change_y', False)
    dict.setdefault('y_min', 0.0)
    dict.setdefault('y_max', 1000.0)
    dict.setdefault('y_scale', 0) # [0, 1, 2, 3, 4, 5, 6] = [Linear, Log10, ln, Log2, Reciprocal(1/T), Probability, Logit]
    dict.setdefault('y_num_format', 1)  #[0, 1, 2, 3, 4, 5] = [automatic, 10000.0, 1e4, 1x10^4, 10k, 1 #dot 10^4]
    dict.setdefault('change_canvas', False)
    dict.setdefault('width', 400)
    dict.setdefault('height', 300)
    dict.setdefault('twoYAxis', False)
    Line_on = 0
    Dots_on = 0

    if fromTable == True:
        t = table
        graphName = existGraphName(dict['graphName'])
        gg = qti.app.newGraph(graphName, 1, 1, 1)
        g = gg.activeLayer()
        for i in range(2, t.numCols()+1):
            c = g.addCurve(t, t.colName(i), 2)
            try:
                if emptyCell( t, i+1, 1):
                    break
            except:
                continue
    else:
        gg = qti.app.currentGraph()
        g = gg.activeLayer()
    ## set curves to right axes
    for i in range(g.numCurves()):
        curve = g.curve(i)
        x_axis = curve.xAxis()
        if dict['changeVerticalAxis'] == 1 and i%2 == 1:
            x_axis = 1
        elif dict['changeVerticalAxis'] == 2 and i%2 == 0:
            x_axis = 1
        elif dict['changeVerticalAxis'] == 3:
            x_axis = 1
        elif dict['changeVerticalAxis'] == 0:
            break
    if dict['Line_on'] == True:
        Line_on = 1
    if dict['Dots_on'] == True:
        Dots_on = 1
    DL = dict['DotsLine']

    x = g.numCurves()
    for i in range(x):
        if DL == 0:
            ii = 2
            color = i
        elif DL == 1:
            ii = i % 2
            color = i/2
        elif DL == 2:
            ii = (i + 1) % 2
            color = i/2
        g.setCurveLineColor(i,  (i%17))
        #g.setCurveLineColor(i,  (color%17))
        #g.setCurveLineColor(i,  0)
        #g.setBrush(i, (color%17))


        g.setCurveLineStyle(i, (i%5+1)*int(Line_on)*(ii-ii/2))
        g.setCurveLineWidth(i, 3)
        c = g.curve(i)
        #c.setSkipSymbolsCount(2)
        try:
            s=c.symbol
        except:
            break
        s = c.symbol()
        s.setSize(10)

        #s.setBrush(QtGui.QBrush(1))
        #s.setPen(QtGui.QPen(0, 3))
        s.setStyle(i%15*int(Dots_on) * abs(ii - 1) + int(Dots_on) * abs(ii -1) -1)
    g.removeTitle()
    if dict['title'] != False:
        g.setTitle(dict['title'])
    g.setAxisTitle(0, dict['y_name'])
    g.setAxisTitle(2, dict['x_name'])
    g.setAxisTitle(1, dict['y2_name'])
    for i in range(4):
        g.setAxisTitleFont(i, QtGui.QFont("DejaVu Sans", 11))
    if dict['change_x'] == True and dict['change_y'] == True:
        g.setAxisNumericFormat(2, dict['x_num_format'], 0 )
        g.setScale(2, dict['x_min'], dict['x_max'], 0.0, 10, 5, dict['x_scale'], 0)
        g.setAxisNumericFormat(0, dict['y_num_format'], 0 )
        g.setScale(0, dict['y_min'], dict['y_max'], 0.0, 10, 5, dict['y_scale'], 0)
    elif dict['change_y'] == True and dict['change_x'] == False:
        g.setAutoScale()
        g.setAxisNumericFormat(0, dict['y_num_format'], 0 )
        g.setScale(0, dict['y_min'], dict['y_max'], 0.0, 10, 5, dict['y_scale'], 0)
    elif dict['change_y'] == False and dict['change_x'] == True:
        g.setAutoScale()
        g.setAxisNumericFormat(2, dict['x_num_format'], 0 )
        g.setScale(2, dict['x_min'], dict['x_max'], 0.0, 10, 5, dict['x_scale'], 0)
    elif dict['change_y'] == False and dict['change_x'] == False:
        g.setAutoScale()

    if dict['Grid_on']:
        g.showGrid()
    for k in range(0,4):
        g.setMajorTicksType(k, 3)
        g.setMinorTicksType(k, 3)
        g.setMajorTicksType(k, 3)
        g.setMinorTicksType(k, 3)
        if k%2!=0:
            g.enableAxisLabels(k, False)
        if dict['twoYAxis']:
            g.enableAxisLabels(1, True)
            g.setScale(1, dict['y_min'], dict['y_max'], 0.0, 10, 5, dict['y_scale'], 0)
            g.setAutoScale()
    if dict['change_canvas'] == True:
        g.setCanvasSize(dict['width'], dict['height'])
    g.setAntialiasing(True, True)
    g.removeLegend()
    g.newLegend()
    return gg


## Function for creating LmFit fit output table..
#  @param minimizer Object : LmFit minimize() function result object
#  @param dict Dict : dictionary with varios parameters values after fit
#  @return Qti Table Object : output table
def outputTable(minimizer, dict):
# fitParams, fitStartParams, rr, data, tName = 'outputPmFit_1', fitRezName = 'undefined',useGlobal = False, globalD = {}):

    ndata = len(dict['yNames'])
    #minimiserresult parameters
    mm = collections.OrderedDict()
    mm['Lmfit results'] = ''
    mm['nfev'] = minimizer.nfev
    mm['success'] = minimizer.success
    mm['errorbars'] = minimizer.errorbars
    mm['message'] = minimizer.message
    mm['ier'] = minimizer.ier
    mm['lmdif_message'] = minimizer.lmdif_message
    mm['nvarys'] = minimizer.nvarys
    mm['ndata'] = minimizer.ndata
    mm['nfree'] = minimizer.nfree
    mm['-----'] = ''
    mm['Fit statistics:'] = ''
    mm['chisqr'] = round_sig(minimizer.chisqr,3)
    mm['redchi'] = round_sig(minimizer.redchi,3)
    mm['aic'] = round_sig(minimizer.aic,3)
    mm['bic'] = round_sig(minimizer.bic, 3)
         
    if ndata == 1:
        mm['R^2'] = round_sig(dict['r2'][0], 3)
        mm['P'] = round_sig(dict['p_values'][0], 3)

        
    mm['------'] = ''    
    mm['Other information'] = ''
    mm['input table'] = dict['intable']
    mm['output table'] = dict['outtable']
    mm['fit result table'] = dict['fittable']
    mm['date/time'] = time.strftime("%Y-%m-%d %H:%M:%S") 
    kk = ['lmfit var',
            'lmfit val',
            '',
            'fit param',
            'value',
            'min',
            'max',
            'vary',
            'global']
    
    gener = {}
    genval = 2
    
    rr = dict['fitvars']
    fitvarsnumber = len(rr)
    tName = dict['outtable']
    fitStartParams = dict['initialvals']
    other_params = dict['otherparams']
    rows = 50
    globalD = dict['globalD']
    if len(globalD) > 0:
        useGlobal = True
    else:
        useGlobal = False
    
    #writes parameters used to generate curves
    '''inc = inputTable.colNames()
    if inputTable.objectName() == 'Example' and 'initPar' in inc and 'value' in inc:
        genval = 3
        for row in range(1, inputTable.numRows() + 1):
            if not emptyCell(inputTable, 'value', row):
                key = inputTable.cellData('initPar', row)
                value = inputTable.cellData('value', row)
                if key in rr:
                    gener.update({key:value})
                    '''
    names = dict['yNames']
    outputTable = qti.app.newTable(tName, rows, 10 + genval * ndata)
    ##parameters names in first column
    outputTable.setColTextFormat(1)
    outputTable.setColName(1, 'param')
    for yi, y in enumerate(rr):
        outputTable.setText(1, yi+1, y)
    if ndata > 1:
        outputTable.setText(1, fitvarsnumber + 2, 'Individual statistics of fit')
        outputTable.setText(1, fitvarsnumber + 3, 'R^2')
        outputTable.setText(1, fitvarsnumber + 4, 'P')
        outputTable.setText(1, fitvarsnumber + 5, 'chisqr')
        outputTable.setText(1, fitvarsnumber + 6, 'redchi')
        offset = fitvarsnumber + 6
    else:
        offset = fitvarsnumber

        
    '''
    if genval == 3:
        for cols in range(ndata):
            outputTable.setColName((cols) * genval + 2, names[cols])
            
            for yi, y in enumerate(rr):
                try:
                    value = gener[y]
                    outputTable.setCellData((cols) * genval + 2, yi + 1, value)
                except KeyError:
                    continue
                    '''
    ## fitted parameters values
    rez = minimizer.params.valuesdict()
    for cols in range(ndata):
        if genval == 3:
            outputTable.setColName((cols) * genval + genval, 'fit_' + str(cols + 1))
        else:
            outputTable.setColName((cols) * genval + 2, names[cols])
        for yi, y in enumerate(rr):
            y = y + '_' + str(cols + 1)
            try:
                value = rez[y]
                outputTable.setCellData((cols) * genval + genval, yi + 1, round_sig(value))
            except KeyError:
                continue
        if ndata > 1:
            outputTable.setCellData((cols) * genval + genval, fitvarsnumber + 3, round_sig(dict['r2'][cols], 4))
            outputTable.setCellData((cols) * genval + genval, fitvarsnumber + 4, round_sig(dict['p_values'][cols],3))
            outputTable.setCellData((cols) * genval + genval, fitvarsnumber + 5, round_sig(dict['chi_sqr'] [cols], 3))
            outputTable.setCellData((cols) * genval + genval, fitvarsnumber + 6, round_sig(dict['red_chi'][cols], 3))

    ## other parameters
    outputTable.setText(1, offset + 2, 'Other parameters')
    for yi, y in enumerate(other_params.keys()):
        outputTable.setText(1, offset + 3 + yi, y)
        outputTable.setCellData(2, offset + 3 + yi, other_params[y])
    
    
    ## adding minimize parameters
    outputTable.setColTextFormat(ndata * genval + 2)
    outputTable.setColTextFormat(ndata *genval + 3)
    for yi, y in enumerate(mm.keys()):
        outputTable.setCellData((ndata) * genval + 2, yi + 1, y)
        outputTable.setCellData((ndata) * genval + 3, yi + 1, mm[y])

    ## adding fitting parameters
    #fitStartParams
    outputTable.setColTextFormat(ndata *genval + 5)
    for yi, y in enumerate(fitStartParams):
        outputTable.setCellData((ndata) * genval + 5 + yi % 5, yi/5 + 1, y)
    if useGlobal:
        for yi, name in enumerate(rr):
            try:
                outputTable.setCellData((ndata) * genval + 5 + 5, yi + 1, globalD[name] )
            except KeyError:
                outputTable.setCellData((ndata) * genval + 5 + 5, yi + 1, False)
    ## naming columns
    for yi, y in  enumerate(kk):
        outputTable.setColName((ndata) * genval + 2 + yi, y)
    return qti.app.currentTable()

## Function for doing logarithmic(10) derivative dy/d(log(x)) on selected QtiPlot table columns.
#  @param tName="logDeriv_1" String : name of output table
#  @paramtable table=qti.app.currentTable() Qti Table Object : input table
#  @return Qti Table Object : output table
def logDerivTable(tName="logDeriv_1", table=qti.app.currentTable()):
    answer = columnToArray()
    x = answer[0]
    data = answer[1]
    tName = existTableName(tName)
    ndata, nx = data.shape
    outputTable = qti.app.newTable(tName, nx, ndata * 2 + 1)
    for ki, k in enumerate(x):
        outputTable.setCellData(1, ki + 1, k)
    der = bp.logDerivative(x, data)
    for ki, k in enumerate(der):
        for yi, dy in enumerate(k):
            outputTable.setCellData(ki * 2 + 2, yi + 1, der[ki][yi] )
            outputTable.setCellData(ki * 2 + 3, yi + 1, data[ki][yi])
    # naming columns
    scols = table.selectedColumns()
    names = []
    for yi, i in enumerate(scols):
        index = table.colIndex(i)
        k = table.colName(index+1)
        if yi > 0:
            names.append('(' + k + ")'")
        names.append(k)
    for cols in range(ndata * 2 + 1):
        outputTable.setColName(cols + 1, names[cols])
    return outputTable

## Output table of derivative characteristics generator.
#  @param paramDict Dictionary : keys and values of derivative characteristics
#  @param param="default" String : simulated parameter name
#  @param initParams={} Dictionary : keys and values of initial parameters
#  @param simParams=qti.app.table("SimParam") Qti Table Object : simulated parameters table
#  @param tName="derivParams_1" String : output table name
#  @return Qti Table Object : output table
def outTabDeriv(paramDict, param = "default", initParams = {}, simParams = qti.app.table("SimParam"), tName = "derivParams_1" ):
    tName = existTableName(tName)
    kk = [param,
            'xHmax',
            'hMax',
            'sHalf',
            'hEnd',
            'MxHmax',
            'MhMax',
            'MsHalf',
            'MhEnd',
            'initP',
            'value']
    col = len(kk)
    row = len(initParams) +3
    if row <= len(paramDict)/8:
        row = len(paramDict)/8
    outputTable = qti.app.newTable(tName, row, col)
    ##simParam parameter values in first column
    empty = 0
    for i in range(simParams.numRows()):
        if not emptyCell(simParams, param, i + 1):
            value = simParams.cell(param, i + 1)
            outputTable.setCellData(1, i +1 - empty, value)
        else:
            empty = empty + 1
    for i in range(6, 11):
        outputTable.setColTextFormat(i)
    ## derivative characteristics float values and success messages
    for i in range(len(paramDict)/8):
        for ki, k in enumerate(kk[1:9]):
            outputTable.setCellData(ki + 2, i + 1, paramDict[k+'_%i' % (i + 1)])
    ## initial parameters
    for yi, y in enumerate(initParams):
        outputTable.setCellData(10, yi + 1, y)
        outputTable.setCellData(11, yi + 1, initParams[y])
    ## writing date and time in the table


    outputTable.setCellData(10, len(initParams) + 2, 'date/time')
    outputTable.setCellData(10, len(initParams) + 3, time.strftime("%Y-%m-%d %H:%M:%S"))
    ## setting column names
    for yi, y in enumerate(kk):
        outputTable.setColName(yi + 1, y)
    return tName
## Function for generating original table names.
#  @param graphName Str : name for checking
#  @return  Str : original table name
def existGraphName(graphName):
    i = 1
    while i < 1000:
        to = -1
        if i > 10:
            to = -2
        if i > 100:
            to = -3
        graphName = graphName[:to] + str(i)
        try:
            g = qti.app.graph(graphName)
            l = g.activeLayer()
        except:
            break
        i=i+1
    return graphName



## Function for creating temperature-pressure phase diagram fit result table.
#  @param x Array : all x data points
#  @param funcname Function : mathematical function for mathematical formula ex. Gaussian
#  @param params Object : parameters values
#  @param experm_table Str : experimental data table name
#  @param selcolindexes Array : 2D array of indexes of selected columns
#  @param fitTableName="fit_result_1" String : output table name
#  @return Str : name of fit result table
def fitResultElip(x, funcname, params, other_params, experm_table, selcolindexes, fitTableName = "fit_result_1"):
    fitTableName = existTableName(fitTableName)
    cols = len(selcolindexes[1])
    xcol = selcolindexes[0][0]
    t1 = qti.app.newTable(fitTableName, len(x), 2 * cols + 1)
    t2 = experm_table
    #selCols = t2.selectedColumns()
    #sc = selcolindexes
    # sc = []
    # for k in selCols:
        # index = t2.colIndex(k)
        # sc.append(index)
    qti.app.updateLog('cols = ' +str(cols)+' \n ') 
    sc = selcolindexes[1]
    for i in range(len(x)):
        t1.setCellData(1, i+1, x[i])
        for k in range(cols):
            #if i == 0 or (i > 0 and x[i] >= x[i-1]):
              #  y = funcname(params, x[i], k, other_params)
            #elif i > 0 and x[i] < x[i-1]:
               # y = funcname(params, x[i], k, other_params)
            y = funcname(params, x[i], k, other_params)
            t1.setCellData((k+1) * 2 + 1, i+1, y)
    for i in range(len(x)):
        for yk, k in enumerate(sc):
            y = t2.cell(k + 1, i+1)
            t1.setCellData((yk + 1) * 2 , i+1, y)
    ##  setting column values
    t1.setColName(1, t2.colName(xcol + 1))
    for yi, i in enumerate(sc):
        value = t2.colName(i + 1)
        t1.setColName((yi + 1 ) * 2 + 1, 'fit%i' % (yi+1))
        t1.setColName(yi + 2, value)
    return fitTableName

## Function for deleting all graphs.
def deleteGraphs():
    f = qti.app.activeFolder()
    lst = f.windows()
    for w in lst:
        if isinstance(w, qti.Graph):
            name = w.objectName()
            t = qti.app.graph(name)
            try:
                t.confirmClose(False)
            except:
                continue
            t.confirmClose(False)
            t.close()
## Function for deleting all fit graphs.
def deleteFitGraphs():
    f = qti.app.activeFolder()
    lst = f.windows()
    for w in lst:
        if isinstance(w, qti.Graph):
            name = w.objectName()
            if name.startsWith('Gfit'):
                t = qti.app.graph(name)
                try:
                    t.confirmClose(False)
                except:
                    continue
                t.confirmClose(False)
                t.close()
## Function for deleting all tables except "initParams", "fitWizard" and "SimParam".
def deleteTables():
    f = qti.app.activeFolder()
    lst = f.windows()
    for w in lst:
        if isinstance(w, qti.Table):
            name = w.objectName()
            if name == "initParams" or name == "fitWizard" or name == "SimParam" or name =="permExample":
                continue
            t = qti.app.table(name)
            try:
                t.confirmClose(False)
            except:
                continue
            t.confirmClose(False)
            t.close()
## Function for deleting all tables fit output tables.
def deleteFitTables():
    f = qti.app.activeFolder()
    lst = f.windows()
    for w in lst:
        if isinstance(w, qti.Table):
            name = w.objectName()
            if name.startsWith('Tfit'):
                t = qti.app.table(name)
                try:
                    t.confirmClose(False)
                except:
                    continue
                t.confirmClose(False)
                t.close()

## Function for deleting all graphs and tables (except "initParams", "fitWizard" and "SimParam").
def deleteGraphsTables():
    deleteGraphs()
    deleteTables()

## Function for deleting all graphs and tables (except "initParams", "fitWizard" and "SimParam").
def deleteFitGraphsTables():
    deleteFitGraphs()
    deleteFitTables()


'''
Code below is used for working with mathematical functions and datasets.
'''
## Function for applying Brents numerical methon .
#  @param func function for Brent method
#  @param y0 Float : Brent method y value
#  @param dict Dictionary:  arguments keys and values
#  @param dict['xa'] Float : Brent method minimum x value
#  @param dict['xb'] Float : Brent method maximum x value
#  @param dict['errTol'] Float : value of y selection tolerance
#  @param dict['max_iter'] Float : number of Brent method loops
#  @param dict['xtol'] Float : value of x selection tolerance
#  @param i=0 Integer : parameter for global fit parameters indexing
#  @return Float : Brent method x value
def brent(func, y0, dict, ii=0, other_params = {}):
    try:
        dic = dict.valuesdict()
    except:
        dic = dict

    xa =                    other_params['xa']
    xb =                    other_params['xb']
    errTol =                other_params['errTol']
    max_iter =              other_params['max_iter']
    xtol =                  other_params['xtol']


    fa = func(dic, xa, ii, other_params) - y0
    fb = func(dic, xb, ii, other_params) - y0

# Function is not defined if in both interval endings is positive
    if(fa * fb >= 0):
#               print "function is not bracketed"
        if(fa < fb): #return(min)
            return xa
        else:
            return xb
# variables switch places
    if(abs(fa) < abs(fb)):
        xa, xb = xb, xa
        fa, fb = fb, fa
# xc assigned as xa and mflag = 1
    xc, fc = xa, fa
    mflag = True

    for i in range(int(max_iter)):
# try to calculat `xs` with reverse quadratic interpolation
# if we can not use Secant rool
        if((fa != fc) and (fb != fc)):
            xs = (xa * fb * fc / ((fa - fb) * (fa - fc)) +
            xb * fa * fc / ((fb - fa) * (fb - fc)) +
            xc * fa * fb / ((fc - fa) * (fc - fb)))
# Secant rool
        else:
            xs = xb - (fb * (xb - xa) / (fb - fa))

        temp2 = (3*xa+xb)/4

# check if xs in interval [(3*xa+xb)/4; b]
        if(((xs < ((temp2)/4)) and (xs > xb)) or
                (mflag == True  and (abs(xs-xb)) >= (abs(xb - xc)/2)) or
                (mflag == False and (abs(xs-xb)) >= (abs(xc -  d)/2)) or
                (mflag == True  and (abs(xb-xc)) < xtol) or
                (mflag == False and (abs(xc - d)) < xtol)):
        #bisection method
            xs = (xa + xb) / 2

            mflag = True
        else:
            mfalg = False

#               fs = Lt(Aa, Bb, xs)
        logic_float = True

        fs = func(dict, xs, ii,other_params) - y0

        if(abs(fs) < errTol):
#                       return (xs, fs, True)
            return xs
        if(abs(xb - xa) < xtol):
#                       return (xs, fs, False)
            return xs
        d = xc
        xc, fc = xb, fb

        if(fa * fs < 0):
            xb, fb = xs, fs
        else:
            xa, fa = xs, fs

        if(abs(fa) < abs(fb)):
            xa, xb = xb, xa
            fa, fb = fb, fa
    #qti.app.updateLog('Lt: ' + str("%.4g"%y0) + '\n')
    #qti.app.updateLog('Pm: ' + str("%.4g"%xs) + '\n')
    return xs

## Function for calculating R^2 in the QtiPlot table.
#  @param table Qti Table Object : input tabble, sequence of columns : x, exp1, fit1, exp2...
#  @return 1D array : list with R^2 values
def Rsquare(table):
    data = tableToArray(table)[1]
    n_data, n_x = data.shape
    i = 0
    answer = []
    while i < n_data:
        gen_sum = 0
        gen_av = 0
        for k in range(n_x):
            gen_sum = gen_sum + data[i][k]
        gen_av = gen_sum / (n_x)
        s_yy = 0
        s_rs = 0
        for k in range(n_x):
            s_yy = s_yy + (data[i][k] - gen_av)**2
            s_rs = s_rs + (data[i + 1][k] - gen_av)**2
        r_sq = s_rs / s_yy
        answer.append(r_sq)
        i = i + 2
    return answer



## Function for calculating logarithmic(10) derivative dy/d(log(x)).
#  @param x Array : x values
#  @param data 2D Array : y values
#  @return 2D Numpy Array : values of logarithmic(10) derivative
def logDerivative(x, data):
    ndata, nx = data.shape
    derivative = []
    for k in data:
        der = []
        for pos, y in enumerate(k):
            if pos < (len(k) - 1):
                d1 = (k[pos+1] - k[pos])/(log(x[pos+1], 10) - log(x[pos], 10))
            d2 = (k[pos] - k[pos-1])/(log(x[pos], 10) - log(x[pos-1], 10))
            if pos == 0:
                der.append(d1)
            elif pos == (len(k) - 1):
                der.append(d2)
            else:
                der.append(0.5*(d1 + d2))
        derivative.append(der)
    derivative = np.array(derivative)
    return derivative


## Function for finding nearest array value to the chosen one.
#  @param array 1 DArray : list of float values
#  @param value Float : chosen value
#  @return Float : nearest value in array
def findNearest(array, value):
    n = [abs(i-value) for i in array]
    idx = n.index(min(n))
    return array[idx]

## Function for making active buttons
#  @param rr String array : list of all parameters which can be changed
#  @return String|String array : selected parameter|parameters
def makeButtons(rr):
    ret = False
    msgBox = QtGui.QMessageBox()
    msgBox.setText('Select parameter')
    buttons = {}
    answer = []
    for i in rr:
        buttons[i] = msgBox.addButton(msgBox.tr(i), QtGui.QMessageBox.ActionRole)
    all = msgBox.addButton(msgBox.tr("ALL"), QtGui.QMessageBox.ActionRole)
    cancel = msgBox.setStandardButtons(QtGui.QMessageBox.Close)
    msgBox.exec_();
    for i in rr:
        if msgBox.clickedButton() == buttons[i]:
            ret = [i]
    if msgBox.clickedButton() == all:
        ret = rr

    return ret

## Function for making data closer to reality.
#  @param data Numpy 2D Array : array with experimental data
#  @return 3D array : array with edited input data in place 0 and value of longest array lenght
def smartData(data):
    answer = []
    ans = []
    max = 0
    for i in data:
        for k in range(1,len(i) - 1):
        #if i[k] > i[k - 1] and i[k] < i [k + 1]:
                #del i[k:]

            if i[k] < i[k - 1]  and i[k + 1] > i[k] and i[k + 1] - i[k] > 0.5 :
                ans.append(i[:k + 1])
                if k +1 > max:
                    max = k + 1
                break

            elif i[k] > i[k - 1]  and i[k + 1] < i[k] and i[k] - i[k + 1] > 0.5 :
                ans.append(i[:k + 1])
                if k +1 > max:
                    max = k + 1
                break

            elif k == len(i) -2:
                ans.append(i)
                max = k + 2

    ans = np.array(ans)
    answer.append(ans)
    answer.append(max)
    return answer


## Function for adding white noise on data.
#  @param array Numpy 2D Array : array with experimental data
#  @param amplitude Float : value of white noise amplitude
#  @return 2D array : experimental data with noise added
def addNoise(array, amplitude):
    ndata, nx = array.shape
    noise =[]
    for i in range(ndata):
        ar = []
        for k in range(nx):
            nois = amplitude * rnd.uniform(-1.0,1.0)
            ar.append(nois)
        noise.append(ar)
    noise = np.array(noise)
    array = array + noise
    return array

## Function for checking Qtilot's active folder objects types.
#  @param nameList String 1D Array : objects names for checking
#  @return Dict : key - object name, value - type of object, if False - there is no such object on current Folder
def objectChecker(nameList):
    f = qti.app.activeFolder()
    lst = f.windows()
    infoDict = {}
    objectDict = {}
    for ob in lst:
        name = str(ob.objectName())
        typ = str(type(ob))
        objectDict.update({name : typ})
    for name in nameList:
        if name in objectDict:
            infoDict.update({name : objectDict[name]})
        else:
            infoDict.update({name : False})
    return infoDict


## Function for creating message if indicated object is not found.
#  @param name String : object name
#  @param type String : object type. Possible options: 'Table', 'Matrix', 'Note', 'Graph'
#  @return QtGui.QMessageBox() : message box with noted missing object
#  @return Boolean : True if object is found, False if it is not
def objectMessanger(name, type):
    trueType = "<class 'qti." + type + "'>"
    if bp.objectChecker([name])[name] != trueType:
        msgBox = QtGui.QMessageBox()
        msgBox.setText("%s  named  '%s'  is not found" %(type, name))
        ret = msgBox.exec_()
        return False
    else:
        return True


## Function for calculating runs test P value.
#  @param table Qti Table Object : input tabble, sequence of columns : x, exp1, fit1, exp2...
#  @return 1 D Array : probabilities
def runs_test(table):
    data = bp.tableToArray(table)[1]
    n_data, n_x = data.shape
    i = 0
    answer = []
    while i < n_data:
        res_all = []
        res_pos = 0.0
        res_neg = 0.0
        seq_num = 1.0
        for k in range(n_x):
            res = data[i + 1][k] - data[i][k]
            if res > 0:
                res_pos = res_pos + 1.0
            else:
                res_neg = res_neg + 1.0
            if k > 0:
                if res * res_all[-1] < 0:
                    seq_num = seq_num + 1.0
            res_all.append(res)
        e_r = (2 * res_neg * res_pos / (res_pos + res_neg)) + 1
        e_v = 2 * res_neg * res_pos * \
            (2 * res_neg * res_pos - res_pos - res_neg) \
            / (((res_pos + res_neg)**2) * (res_pos + res_neg - 1))
        if e_v != 0:
		    z_sd = (seq_num - e_r) / sqrt(e_v)
        else:
            z_sd = 100000000000000
        p_prob = 2 * st.norm.sf(abs(z_sd))
        answer.append(p_prob)
        i = i + 2
    return answer

## Function for calculating runs test P value.
#  @param table Qti Table Object : input tabble, sequence of columns : x, exp1, fit1, exp2...
#  @return 1 D Array : probabilities
def chi_squares(table, n_varys):
    data = bp.tableToArray(table)[1]
    n_data, n_x = data.shape
    i = 0
    answer = [[], []]
    while i < n_data:
        chi_sqr = 0
        for k in range(n_x):
            res = data[i + 1][k] - data[i][k]
            chi_sqr = (chi_sqr + res**2)
        answer[0].append(chi_sqr)
        answer[1].append(chi_sqr / (n_x - n_varys))
        i = i + 2
    return answer

## Fuction for adding values to Qtiplot's resultsLog
#  @param log_inform=dict Float : amplitude of white noise
def results_logger(log_inform=dict):
    qti.app.resultsLog().append('\n')
    qti.app.resultsLog().append("Date: " +time.strftime("%Y-%m-%d %H:%M:%S") + '\n')
    if log_inform != {}:
        for key, value in log_inform.items():
            qti.app.resultsLog().append(key + ": " +str(value) + '\n')

            
## Function for one Gaus model.
#  @param dict Dictionary : Pm function arguments
#  @param dict['xMin'] Float : minimum value
#  @param dict['xMax'] Float : maximum value
#  @param dict['rate'] Int :  sample rate
#  @param scale='log10' String : step scale
#  @return numpy 1D array : list of steped values
def xGen(dict, scale = 'log10'):
    xMin = dict['Lt_min']
    xMax = dict['Lt_max']
    rate = dict['Lt_range']
    answer = []
    b1 = log(xMin, 10)
    b2 = log(xMax, 10)
    if scale == 'log10':
        for i in range(int(rate)):
            value = 10 ** (b1 + (i) * (b2 - b1) / (rate - 1.))
            answer.append(value)
    elif scale == 'linear':
        for i in range(int(rate)):
            value = xMin + i*(xMax-xMin)/(rate-1.)
            answer.append(value)
    answer = np.array(answer)
    return answer