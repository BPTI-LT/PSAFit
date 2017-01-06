"""
    v0.6.1  
    This is a package  genPmLt  which holds function used in protein dose 
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

## Function for calculating derivative characteristics.
#  @param x 1D Numpy Array : list of x(Ligand concentration) Float values
#  @param data 2D Numpy Array : list of y(melting pressure) Float values
#  @return dic Dictionary : derivative characteristics
#  @return dic['hMax_*'] Float : maximum derivative value
#  @return dic['xHmax_*'] Float : x at maximum derivative value
#  @return dic['sHalf_*'] Float : half width at half maximum
#  @return dic['hEnd_*'] Float : derivative value at the end of the curve
#  @return dic['MhMax_*'] String : success message of hMax
#  @return dic['MxHmax_*'] String : success message of xHmax
#  @return dic['MsHalf_*'] String : success message of sHalf
#  @return dic['MhEnd_*'] String : success message of hEnd
def derivPar(x, data):
    ## messages of parameters calculations success
    answer = bp.smartData(data)
    ## editing data curves - values with false brent method are eliminated
    sdata = answer[0]
    sx = x[:answer[1]]
    ## calculating logarithmic derivative of
    derivative  = bp.logDerivative(sx, sdata)
    answer = bp.smartData(data)
    dic = {}
    for yi, y in enumerate(derivative):
        dic.update({'MhMax_%i' % (yi + 1):'OK'})
        dic.update({'MxHmax_%i' % (yi + 1):'OK'})
        dic.update({'MsHalf_%i' % (yi + 1):'OK'})
        dic.update({'MhEnd_%i' % (yi +1 ):'OK'})
        dic.update({'hMax_%i' % (yi + 1): max(y)})
        indexMax = y.tolist().index(dic['hMax_%i' % (yi + 1)])
        dic.update({'xHmax_%i' % (yi + 1) : sx[indexMax]})
        array = y[:indexMax]
        hHalf = bp.findNearest(array, dic['hMax_%i' % (yi+1)]/2)
        indexHalf = y.tolist().index(hHalf)
        xHhalf = sx[indexHalf]
        dic.update({'sHalf_%i' % (yi + 1) : (dic['xHmax_%i' % (yi + 1)] - xHhalf)})
        dic.update({'hEnd_%i' % (yi + 1) : y[-1]})
        if dic['hMax_%i' % (yi+1)] == y[-1]:
            message = "The maximum of derivative value might be unreached"
            if y[-1] == y[-2]:
                message = "There is a possible linear growth in the end of the curve. "
            dic['MhMax_%i' % (yi + 1)] = message
            dic['MxHmax_%i' % (yi + 1)] = message
        if (indexMax - indexHalf) < 6:
            dif = indexMax - indexHalf
            dic['MsHalf_%i' % (yi + 1)] = "There is just %i samples between maximum and half of the maximum derivative value. "  \
                                    "It is recommended to look at the graphical output. " %dif
        if len(sx) != len(x):
            dic['MhEnd_%i' % (yi + 1)] = "Curve has been generated from " + "{:.2e}".format(sx[0]) + " to " + "{:.2e}".format(sx[-1]) \
                                     + " (instead of " + "{:.2e}".format(x[-1]) + ") of x values."
    return dic


## Function for showing logarithmic derivative parameters in QtiPlot table.
#  @param SimParam Qti Table Object : template table
#  @param initParams Qti Table Object : template table
#  @return Qti Table Object : table with derivative parameters
def derivatParams(model = 'N'):
    #                                Parameters initialization
    t = qti.app.table("initParams")
    SimParam = qti.app.table("SimParam")
    SimParams = bp.simParamDict(SimParam)
    output_table_name = "TderPar" + model + "_1"
    genlog = ''
    # Code. Do not change following text
    rr = bp.colNames(SimParam)
    if model == 'N':
        ## erase ligand binding to protein unfolded state parameters from options
        for i in rr:
            if i[1:3] == 'bu':
                rr.remove(i)
    d = bp.initParamsDict(t)
    x = bp.xGen(d)
    for k in range(100):
        r = bp.makeButtons(rr)
        d = bp.initParamsDict(t)
        if model == 'N':
            ## erase ligand binding to protein unfolded state parameters from initial parameters
            for i in ['DbuV', 'Kbu', 'DbuBeta']:
                d.pop(i, None)
        
        if r == False:
            break
        else:
            for param in r:
                
                ## generation of data curves
                if model == 'N':
                    data = bp.PmLtUpper(x, d, param, SimParams, d)
                else:
                    data = bp.PmLtFullSim(x, d, param, SimParams,  d)
                answer = bp.derivPar(x, data)
                output_table_name = bp.outTabDeriv(answer, param, d, SimParam, output_table_name)
                genlog += "    " + param + ": \n" + '      values: ' + str(SimParams[param]) + '\n'
                genlog += '        output table: ' + output_table_name + '\n' 
                
    log = "Date: " +time.strftime("%Y-%m-%d %H:%M:%S") + '\n'   
    log += "[[Generation of dose curve derivative parameters]] \n"
    log += "  Bindig model: " + model + '\n'
    log += "[[Varied parameters]] \n"
    log += genlog
    qti.app.resultsLog().append(log)
## Function for generating Pm(Lt) curves.
#  @param noiseAmplitude Float : amplitude of white noise
#  @param SimParam Qti Table Object : template table
#  @param initParams Qti Table Object : template table
#  @return Qti Table Object : results table(data curves)
#  @return Qti Graph Object : results graph(plot curves)
def generPmLt(noiseAmplitude, model = 'N'):
    #                                Parameters initialization
    t = qti.app.table("initParams")
    SimParam = qti.app.table("SimParam")
    SimParams = bp.simParamDict(SimParam)

    pdc ={
        'y_name' : "<i>P<sub>m</sub></i>, MPa",
        'x_name' : "<i>L<sub>t</sub></i>, M",
        'change_x' : True,
        'x_scale' : 1,
        'x_num_format' : 5,
        'Dots_on' : False}

    #                                                               Code. Do not change following text
    
  
    rr = bp.colNames(SimParam)
    if model == 'N':
        ## erase ligand binding to protein unfolded state parameters from options
        for i in rr:
            if i[1:3] == 'bu':
                rr.remove(i)
    d = bp.initParamsDict(t)
    x = bp.xGen(d)
    for k in range(100):
        r = bp.makeButtons(rr)
        d = bp.initParamsDict(t)
        if model == 'N':
            ## erase ligand binding to protein unfolded state parameters from initial parameters
            for i in ['DbuV', 'Kbu', 'DbuBeta']:
                d.pop(i, None)
        if r == False:
            break
        else:
            genlog = ''
            for param in r:
                
                genlog += "    " + param + ": \n" + '      values: ' + str(SimParams[param]) + '\n'
                if model == 'N':
                    data = bp.PmLtUpper(x, d, param, SimParams, d)
                else:
                    data = bp.PmLtFullSim(x, d, param, SimParams, d)
                answer = bp.smartData(data)
                sdata = answer[0]
                sdata1 = sdata
                # add noise
                sdata = bp.addNoise(sdata, noiseAmplitude)
                sx = x[:answer[1]]
                pdc.update({'x_min' : sx[0],
                        'x_max': sx[-1],
                        'graphName': 'G' + param + model + '_1'})
                output_table = bp.generatorResultTable(sx, sdata, d, param,  'T'+param+model+'_1')
                output_graph = bp.Geditor(dict=pdc, fromTable=True, table=output_table)
                genlog += '        output table: ' + output_table.objectName() + '\n'
                genlog += '        output graph: ' + output_graph.objectName() + '\n'
               # bp.results_logger(log_inform)
       
    log = "Date: " +time.strftime("%Y-%m-%d %H:%M:%S") + '\n'   
    log += "[[Generation of dose curves]] \n"
    log += "  Bindig model: " + model + '\n'
    log += "  Noise amplitude: " + str(noiseAmplitude) + " MPa \n"
    log += "[[Varied parameters]] \n"
    log += genlog
    qti.app.resultsLog().append(log)

## Function for calculating logarithmic derivative.
#  @param SimParam Qti Table Object : template table
#  @param initParams Qti Table Object : template table
#  @return Qti Table Object : results table(data curves)
#  @return Qti Graph Object : results graph(plot curves)
def logDerivat(model = 'N'):
    input_table = qti.app.currentTable()
    output_table_name = "Tderiv" + model + "_1"
    a = bp.logDerivTable(tName=output_table_name , table=input_table)
    pdc ={
        'y_name':"<i>P<sub>m</sub></i>, MPa",
        'x_name':"<i>L<sub>t</sub></i>, M",
        'y2_name':"d<i>P<sub>m</sub></i> / d[log<sub>10</sub><i>(L<sub>t</sub>)</i>]",
        'change_x': True,
        'x_scale' : 1,
        'x_num_format' : 5,
        'DotsLine' : 0,
        'Dots_on' : False,
        'x_max' : 100000,
        'twoYAxis' : True,
        'changeVerticalAxis' : 2,
        'graphName' : 'GderivN' + model + "_1"}
    log_inform = {
        'Action' : 'Dose curves differentiation',
        'Input table' : input_table.objectName(),
        'Output table' : output_table_name,
        'Output graph': 'Gderiv' + model + "_1"
    }
    bp.Geditor(dict=pdc, fromTable=True, table=a)
    bp.results_logger(log_inform)
