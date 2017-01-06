"""
    v0.6.1  
    This is a package  PmLt  which holds function used in protein dose curves 
    programs.

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

## Function for calculating Pm(Lt) curve .
#  @param Lt Float : value of ligand concentration
#  @param dict Dictionary : Pm function arguments
#  @param param="none" String : name of parameter from SimParam table
#  @param SimParam Dict : dictionary with varied parameter values
#  @return Float : value of protein/ligand solution melting pressure
#  @return Numpy 2D Array : values of protein/ligand solution melting pressure, if input is array
def PmLtFullSim(Lt, dict, param="none", SimParam={}, other_params={}):

    Pms = []
    PmAll = []

    if param == "none":
        Pm_range_calc = findPmRange(Pm_low, LtFull, dict, ii=0, other_params=other_params)
        xa = Pm_range_calc['low']
        xb = Pm_range_calc['high']
        dict.update({'xa':xa})
        dict.update({'xb':xb})
        if isinstance(Lt, int) or isinstance(Lt, float):
            Pm = bp.brent(LtFull, Lt, dict, ii=0, other_params= other_params)
            return Pm
        else:
            for i in range(len(Lt)):
                Pm = bp.brent(LtFull, Lt[i], dict, ii=0, other_params= other_params)
                Pms.append(Pm)
            PmAll.append(Pms)
            answer = np.array(PmAll)
            return answer
    else:
        for i in SimParam[param]:
            Pms = []
            if param != "none":
                old = dict[param]
                dict.update({param:i})
                Pm_range_calc = findPmRange(Pm_low, LtFull, dict, ii=0, other_params=other_params)
                xa = Pm_range_calc['low']
                xb = Pm_range_calc['high']
                dict.update({'xa':xa})
                dict.update({'xb':xb})

            if isinstance(Lt, int) or isinstance(Lt, float):
                Pm = bp.brent(LtFull, Lt, dict, ii=0, other_params= other_params)
                return Pm
            else:
                for k in range(len(Lt)):
                    Pm = bp.brent(LtFull, Lt[k], dict, ii=0, other_params= other_params)
                    Pms.append(Pm)
            PmAll.append(Pms)
            # update dict to old value
            dict.update({param:old})
        answer = np.array(PmAll)
        return answer



## Function for calculating ligand concentration proteint melting pressure curves Lt(Pm).
#  @param x Float : protein melting pressure value
#  @param dict Dictionary : modeling parameters keys and values
#  @param dict['Pr'] Float : reference pressure
#  @param dict['R'] Float : gas constant
#  @param dict['T'] Float : temperature
#  @param dict['Mt'] Float : protein concentration
#  @param dict['DuV'] Float : change of protein unfolding volume
#  @param dict['DbnV'] Float : change of ligand - native protein binding volume
#  @param dict['DbnV'] Float : change of ligand - unfolded protein binding volume
#  @param dict['DuBeta'] Float : change of protein unfolding isothermic compressibility
#  @param dict['DbnBeta'] Float : change of protein ligand - native protein binding isothermic compressibility
#  @param dict['DbuBeta'] Float : change of protein ligand - unfolded protein binding isothermic compressibility
#  @param dict['Kbn'] Float : ligand - native protein binding equilibrium constant; it is possible to use dict['DbnG'] insted of this
#  @param dict['Kbu'] Float : ligand - unfolded protein binding equilibrium constant; it is possible to use dict['DbnG'] insted of this
#  @param dict['DuG'] Float :  change of protein unfolding Gibbs free energy;
#  @param dict['DbnG'] Float :  change of ligand - native protein binding Gibbs free energy; this parameter is active if there is no dict['Kbn'] parameter value
#  @param dict['DbuG'] Float :  change of ligand - unfolded protein binding Gibbs free energy; this parameter is active if there is no dict['Kbu'] parameter value
#  @param i=0 Int :  iterations parameter
#  @return Float : ligand concentration value
def LtFull(dict, x, i=0, other_params={}):
    try:
        dic = dict.valuesdict()
    except:
        dic = dict
    try:
        Pr =                    other_params['Pr']
        R =                     other_params['R']
        T =                     other_params['T']
        Mt =                    dic['Mt']

        DuV =                   dic['DuV']
        DbnV =                  dic['DbnV']
        DbuV =                  dic['DbuV']

        DuBeta =                dic['DuBeta']
        DbnBeta =               dic['DbnBeta']
        DbuBeta =               dic['DbuBeta']
        DuG =                   dic['DuG']
        Kbn =                   dic['Kbn']
        DbnG =                  (-1*log(Kbn)*R*T)
        Kbu =                   dic['Kbu']
        DbuG =                  (-1*log(Kbu)*R*T)

    except KeyError:
        Pr =                    other_params['Pr']
        R =                     other_params['R']
        T =                     other_params['T']
        Mt =                    dic['Mt_%i' % (i+1)]

        DuV =                   dic['DuV_%i' % (i+1)]
        DbnV =                  dic['DbnV_%i' % (i+1)]
        DbuV =                  dic['DbuV_%i' % (i+1)]

        DuBeta =                dic['DuBeta_%i' % (i+1)]
        DbnBeta =               dic['DbnBeta_%i' % (i+1)]
        DbuBeta =               dic['DbuBeta_%i' % (i+1)]
        DuG =                   dic['DuG_%i' % (i+1)]
        Kbn =                   dic['Kbn_%i' % (i+1)]
        DbnG =                  (-1*log(Kbn)*R*T)
        Kbu =                   dic['Kbu_%i' % (i+1)]
        DbuG =                  (-1*log(Kbu)*R*T)

    ## Equilibrium constants
    # @param KbnFunc Float : ligand binding to protein native state
    # @param KbuFunc Float : ligand binding to protein unfolded state
    # @param KuFunc Float : protein unfolding
    Kbnpower = -((DbnG + DbnV*(x-Pr)+(DbnBeta/2)*(x-Pr)*(x-Pr))/(R*T))
    Kbupower = -((DbuG+DbuV*(x-Pr)+DbuBeta/2*(x-Pr)*(x-Pr))/(R*T))
    Kupower = -((DuG + DuV*(x-Pr)+(DuBeta/2)*(x-Pr)*(x-Pr))/(R*T))
    #qti.app.resultsLog().append("Kbnpower = %s" %Kbnpower)
    #qti.app.resultsLog().append("Kbupower = %s" %Kbupower)
    #qti.app.resultsLog().append("Kupower = %s" %Kupower)
    if(Kbnpower > 700 or Kbupower > 700 or Kupower > 700):
        return 0.1
    KbnFunc = exp(Kbnpower)
    KbuFunc = exp(Kbupower)
    KuFunc = exp(Kupower)
    yyy = KuFunc*KbuFunc-KbnFunc
    nnn = KuFunc*(KbuFunc-KbnFunc)
    #qti.app.resultsLog().append("KuFunc*KbuFunc-KbnFunc = %s" %yyy)
    #qti.app.resultsLog().append("KuFunc*(KbuFunc-KbnFunc) = %s" %nnn)
    if yyy == 0 or nnn == 0:
        return 0.1
    return ((1-KuFunc)*(Mt/2.*(KbnFunc+KbuFunc*KuFunc)/(KuFunc*(KbuFunc-KbnFunc))+1/(KuFunc*KbuFunc-KbnFunc)))

## Lt(Pm) function with protein unfolding and ligand-native protein interactions parameters.
#  @param dict Dictionary|Class : function arguments
#  @param x Float : Pm(melting pressure) value
#  @param x Int : iterations parameter
#  @return Float : value of Lt
def LtUpper(params, x, i=0, other_params={}):
    a={}
    gl=0
    try:
        a = params.valuesdict()
    except:
        a = params
    #qti.app.resultsLog().append("dicLtUpper =" + str(a)) 
    try:
        Pr =            other_params['Pr']
        R =             other_params['R']
        T =             other_params['T']
        DuG =           a['DuG']
        DuV =           a['DuV']
        DuBeta=         a['DuBeta']
        DbnV =          a['DbnV']
        DbnBeta =       a['DbnBeta']
        Mt      =               a['Mt']
        Kbn =           a['Kbn']
    except:
        Pr =            other_params['Pr']
        R =             other_params['R']
        T =             other_params['T']
        DuG =           a['DuG_%i' % (i+1)]
        DuV =           a['DuV_%i' % (i+1)]
        DuBeta=         a['DuBeta_%i' % (i+1)]
        DbnV =          a['DbnV_%i' % (i+1)]
        DbnBeta =       a['DbnBeta_%i' % (i+1)]
        Mt      =               a['Mt_%i' % (i+1)]
        Kbn =           a['Kbn_%i' % (i+1)]
    DbnG =          (-1*log(Kbn)*R*T)
    Kbpower = -((DbnG + DbnV*(x-Pr)+(DbnBeta/2)*(x-Pr)*(x-Pr))/(R*T))
    Kupower = -((DuG + DuV*(x-Pr)+(DuBeta/2)*(x-Pr)*(x-Pr))/(R*T))
        #qti.app.updateLog("Kbnpower = %s" %x + '\n')
    #qti.app.resultsLog().append("Kbupower = %s" %Kbupower)
    #qti.app.resultsLog().append("Kupower = %s" %Kupower)
    #To avoid math overflow
    PowerLimit = 300.0
    if(Kbpower > PowerLimit and Kupower > PowerLimit):
        Kpower = Kupower - Kbpower
        if(Kpower > PowerLimit):
            Kpower = PowerLimit
        model = exp(Kpower)
        #qti.app.resultsLog().append("model1= %s" %model)
    elif(Kupower > PowerLimit):
        Kupower = PowerLimit
        KbFunc = exp(Kbpower)
        KuFunc = exp(Kupower)
        if(KbFunc <> 0.):
            model = KuFunc/KbFunc
        else:
            model = KuFunc
        #qti.app.resultsLog().append("model2= %s" %model)
    elif(Kbpower > PowerLimit):
        KuFunc = exp(Kupower)
        if(KuFunc <> 0.):
            model = (KuFunc - 1.0) * (Mt/2.0/KuFunc)
        else:
            model = exp(PowerLimit)
            #qti.app.resultsLog().append("Pm errr for Lt = %s" %x)
        #qti.app.resultsLog().append("model3= %s" %model)
    else:
        KbFunc = exp(Kbpower)
        KuFunc = exp(Kupower)
        if(KbFunc <> 0. and KuFunc <> 0.):
            model = (((Mt/(2.0*KuFunc))+(1.0/KbFunc))*(KuFunc - 1.0))
        else:
            model = exp(PowerLimit)
            #qti.app.resultsLog().append("Pm errr for Lt = %s" %x)
        #qti.app.resultsLog().append("model4= %s" %model)
   
    return model
#                                                                                                               Sigmoidal function                                                      #



## Function for finding protein melting Pressure(equal native and unfolded protein concetration).
#  @param dict Dictionary : Pm function arguments
#  @param dict['DuG'] Float : change of protein unfolding Gibbs free energy
#  @param dict['DuV'] Float : change of protein unfolding volume change
#  @param dict['DuBeta'] Float : change of protein unfolding isotermic compressibility
#  @param dict['Pr'] Float : reference - atmospheric pressure
#  @param i=0 Int : iterations parameter
#  @return Float : value of melting pressure
def Pm_low(dict, i=0, other_params = {}):
    try:
        DuG =                   dict['DuG']
        DuV =                   dict['DuV']
        DuBeta=                 dict['DuBeta']
        Pr =                    other_params['Pr']
    except KeyError:
        DuG =                   dict['DuG_%i' % (i+1)]
        DuV =                   dict['DuV_%i' % (i+1)]
        DuBeta=                 dict['DuBeta_%i' % (i+1)]
        Pr =                    other_params['Pr']
    if DuBeta == 0:
        if DuV == 0:
            return 0
        else:
            return -1*(DuG/DuV)+Pr
    elif (DuV*DuV-2*DuBeta*DuG) < 0:
        return 0
    else:
        return -1*((sqrt(DuV*DuV-2*DuBeta*DuG)+DuV)/DuBeta)+Pr










## Function for Pm(Lt) curve finding (for fit).
#  @param Lt Float : value of ligand concentration
#  @param params lmfit.parameters|Dictionary : Pm function arguments
#  @param params['DuG'] Float : protein unfolding free Gibbs energy change
#  @param params['DuV'] Float : protein unfoldind volume change
#  @param params['DuBeta'] Float : protein unfolding isotermic compressibility change
#  @param params['Kbn'] Float : ligand and native protein binding equilibrium constant
#  @param params['DbnG'] Float : ligand and native protein binding free Gibbs energy change
#  @param params['DbnV'] Float : ligand and native protein binding volume change
#  @param params['DbnBeta'] Float : ligand and native protein binding isotermic compressibility change
#  @param params['Kbu'] Float : ligand and unfolded protein binding equilibrium constant
#  @param params['DbuG'] Float : ligand and unfolded protein binding free Gibbs energy change
#  @param params['DbuV'] Float : ligand and unfolded protein binding volume change
#  @param params['DbuBeta'] Float : ligand and unfolded protein binding isotermic compressibility change
#  @param params['xa'] Float : minimum pressure in case of downer function
#  @param params['xb'] Float : maximum pressure in case of upper function
#  @param params['Pr'] Float : atmospheric pressure
#  @param params['T'] Float : temperature
#  @param params['R'] Float : gas constant
#  @param i Integer : iterations parameter for more then one curve
def PmLtFull(Lt, params, i=0, other_params={}):
    Pm_range_calc = findPmRange(Pm_low, LtFull, params, i, other_params)
    xa = Pm_range_calc['low']
    xb = Pm_range_calc['high']
    other_params.update({'xa':xa})
    other_params.update({'xb':xb})
    Pm = bp.brent(LtFull, Lt, params, i, other_params)
    return Pm



## Function for calculating Pm(Lt) curve (Ligand-Native protein binding model).
#  @param Lt Float : value of ligand concentration
#  @param dict Dictionary : Pm function arguments
#  @param param="none" String : name of parameter from SimParam table
#  @param SimParam Dict : dictionary with varied parameters values
#  @return Float : value of protein/ligand solution melting pressure, if input is float
#  @return Numpy 2D Array : values of protein/ligand solution melting pressure, if input is array
def PmLtUpper(Lt, dict, param="none", SimParam={}, other_params = {}):

    Pms = []
    PmAll = []

    
    #if we do not vary parameters
    if param == "none":

        if isinstance(Lt, int) or isinstance(Lt, float):
            Pm = bp.brent(LtUpper, Lt, dict, ii=0, other_params = other_params)

            return Pm
        else:
            for i in range(len(Lt)):
                Pm = bp.brent(LtUpper, Lt[i], dict, ii=0, other_params = other_params)
                Pms.append(Pm)
            PmAll.append(Pms)
            answer = np.array(PmAll)
            return answer
    #if we vary one or all parameters
    else:
        for i in SimParam[param]:
            Pms = []
            if param != "none":
                old = dict[param]
                dict.update({param:i})

            if isinstance(Lt, int) or isinstance(Lt, float):
                Pm = bp.brent(LtUpper, Lt, dict, 0, other_params = other_params)
                return Pm
            else:
                for k in range(len(Lt)):
                    Pm = bp.brent(LtUpper, Lt[k], dict, 0, other_params = other_params)
                    Pms.append(Pm)
            PmAll.append(Pms)
            # update dict to old value
            dict.update({param:old})
        answer = np.array(PmAll)
        return answer



## function for finding P_m range.
#  @param func_Pm_low Function : function for P_m value without ligand
#  @param func_Lt Function : function L_t(P_m)
#  @param dic['DuG'] Float : protein unfolding free Gibbs energy change
#  @param dic['DuV'] Float : protein unfoldind volume change
#  @param dic['DuBeta'] Float : protein unfolding isotermic compressibility change
#  @param dic['Kbn'] Float : ligand and native protein binding equilibrium constant
#  @param dic['DbnG'] Float : ligand and native protein binding free Gibbs energy change
#  @param dic['DbnV'] Float : ligand and native protein binding volume change
#  @param dic['DbnBeta'] Float : ligand and native protein binding isotermic compressibility change
#  @param dic['Kbu'] Float : ligand and unfolded protein binding equilibrium constant
#  @param dic['DbuG'] Float : ligand and unfolded protein binding free Gibbs energy change
#  @param dic['DbuV'] Float : ligand and unfolded protein binding volume change
#  @param dic['DbuBeta'] Float : ligand and unfolded protein binding isotermic compressibility change
#  @param dic['xa'] Float : minimum pressure in case of downer function
#  @param dic['xb'] Float : maximum pressure in case of upper function
#  @param dic['Pr'] Float : atmospheric pressure
#  @param dic['T'] Float : temperature
#  @param dic['R'] Float : gas constant
#  @param ii=0 : iterations parameter
#  @return range Dict : dictionary with upper or downer function pressure range values
def findPmRange(func_Pm_low, func_Lt, dic, ii=0, other_params = {}):
    try:
        dic = dic.valuesdict()
    except:
        dic = dic
    try:
        R =                             other_params['R']
        T =                                     other_params['T']
        DuG =                   dic['DuG']
        DuV = dic['DuV']
        DuBeta = dic['DuBeta']
        Pr = other_params['Pr']
        xb = other_params['xb']
        DbnBeta = dic['DbnBeta']
        DbuBeta = dic['DbuBeta']
        DbnV = dic['DbnV']
        DbuV = dic['DbuV']
        Kbn =                   dic['Kbn']
        Kbu =                   dic['Kbu']


    except KeyError:
        R =                             other_params['R']
        T =                                     other_params['T']
        DuG =                   dic['DuG_%i' % (ii+1)]
        DuV = dic['DuV_%i' % (ii+1)]
        DuBeta = dic['DuBeta_%i' % (ii+1)]
        Pr = other_params['Pr']
        xb = other_params['xb']
        DbnBeta = dic['DbnBeta_%i' % (ii+1)]
        DbuBeta = dic['DbuBeta_%i' % (ii+1)]
        DbnV = dic['DbnV_%i' % (ii+1)]
        DbuV = dic['DbuV_%i' % (ii+1)]
        Kbn =                   dic['Kbn_%i' % (ii+1)]
        Kbu =                   dic['Kbu_%i' % (ii+1)]



    dPm_p = func_Pm_low(dic, ii, other_params) + (func_Pm_low(dic, ii, other_params)/10e5)
    dPm_n = func_Pm_low(dic, ii, other_params) - (func_Pm_low(dic, ii, other_params)/10e5)
    Ltp = func_Lt(dic, dPm_p, ii, other_params)
    Ltn = func_Lt(dic, dPm_n, ii, other_params)
    Pm = func_Pm_low(dic, ii, other_params)
    #extraordinary Pm values for NU model

    DbnG = (-1.*log(Kbn)*R*T)
    DbuG = (-1.*log(Kbu)*R*T)
    A1 = 0.5*(DbnBeta - DbuBeta)
    # qti.app.resultsLog().append("DbnBeta = %s" %DbnBeta)
    # qti.app.resultsLog().append("DbuBeta = %s" %DbuBeta)
    B1 = DbnV - DbuV
    C1 = DbnG - DbuG
    A2 = A1-0.5*(DuBeta)
    B2 = B1 - DuV
    C2 = C1 - DuG
    # qti.app.resultsLog().append("A1 = %s" %A1)
    # qti.app.resultsLog().append("B1 = %s" %B1)
    # qti.app.resultsLog().append("C1 = %s" %C1)
    # qti.app.resultsLog().append("A2 = %s" %A2)
    # qti.app.resultsLog().append("B2 = %s" %B2)
    # qti.app.resultsLog().append("C2 = %s" %C2)
    extraPm = []
    #qti.app.resultsLog().append("Pm0 = %s" %Pm)
    if (A1 <> 0):
        D = B1*B1 -4.*A1*C1
        if D > 0.:
            x1 = -1.*(sqrt(B1*B1 -4.*A1*C1) +B1)/2./A1
            x2 = (sqrt(B1*B1 -4.*A1*C1) - B1)/2./A1
            extraPm.append(x1)
            extraPm.append(x2)
        #qti.app.resultsLog().append("D1 = %s" %D)
        # qti.app.resultsLog().append("x2 = %s" %x2)
    elif(B1 <> 0):
        x1 = -C1/B1
        extraPm.append(x1)
        #qti.app.resultsLog().append("extrax1 = %s" %x1)

    if (A2 <> 0):
        D = B2*B2 -4.*A2*C2
        if D > 0.:
            x3 = -1*(sqrt(B2*B2 -4.*A2*C2) +B2)/2./A2
            x4 = (sqrt(B2*B2 -4.*A2*C2) - B2)/2./A2
            extraPm.append(x3)
            extraPm.append(x4)
        #qti.app.resultsLog().append("D2 = %s" %D)
        # qti.app.resultsLog().append("x4 = %s" %x4)

    elif(B2 <> 0):
        x1 = -C2/B2
        extraPm.append(x1)
        #qti.app.resultsLog().append("extrax21 = %s" %x1)
    temparray = []
    tempPmmin = 1000000.
    for temp in extraPm:
        if temp > Pm:
            if temp < tempPmmin:
                tempPmmin = temp


    #END extraordinary Pm values for NU model
    LTl = func_Lt(dic, Pm, ii, other_params)
    r1 = Ltp - LTl
    r2 = Ltn - LTl
    range = dict(low=0, high=0)
    if((r1 > 0) and (r2 < 0)):
        range['low'] = Pm
        #Atimam 0.00000000001, kad negautume begalybes sitame rezyje
        range['high'] = tempPmmin + Pr - 0.00000000001
        #qti.app.resultsLog().append("Pmhigh = %s N binder" %range['high'])
    elif((r1 < 0) and (r2 > 0)):
        range['low'] = Pr
        range['high'] = Pm
        #qti.app.resultsLog().append("binder = %s U binder" %r1)
    else:
        qti.app.resultsLog().append("binder error")

    return range
