"""
    v0.6.1 
    This is a package  psaFit  which manages other packages of protein dose 
    curves and pressure-temperature phase diagram fit and generation modules.
    
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
import sys

paths = sys.path
for y in paths:
    if y[-7:] == 'qtiplot':
        path = y
        break
    else:
        path = paths[-1]

## Add path to Modules folder from QtiPlot Python Configuration Files directory
directory = path +'/modules'


sys.path.append(directory)
message = {'basic':'checking',
                        'PmLt':'checking',
                        'genPmLt':'checking',
                        'fitPmLt':'checking',
                        'fitPmTm':'checking'                            }
from basic import*

try:
    from PmLt import*
    message.update({'PmLt':'imported'})
except ImportError:
    message.update({'PmLt':'no module'})
try:
    from genPmLt import*
    message.update({'genPmLt':'imported'})
except ImportError:
    message.update({'genPmLt':'no module'})
try:
    from fitPmLt import*
    message.update({'fitPmLt':'imported'})
except ImportError:
    message.update({'fitPmLt':'no module'})
try:
    from fitPmTm import*
    message.update({'fitPmTm':'imported'})
except ImportError:
    message.update({'fitPmTm':'no module'})
