#!/usr/bin/env python

"""
emissions.py 

Module for estimating emissions based on road network traffic conditions.
Estimations are based on: Bureau of Public Roads (BPR) curve, COPERT_v5, and TRANSYT7f.

Usage:
    Initialize_Speed_Interpolation_Functions(Method='old')

    (cost, unit) = ExternalCost(Speed, Speed_Units, Method = 'TRANSYT7f')
    (EmissionsFactor, unit) = COPERT5(Speed, Speed_Units, Polutant)
    (EmissionsFactor, unit) = TRANSYT7f(Speed, Speed_Units, Polutant)

    (speed, unit) = Speed_via_Density(Density,Density_Units,FreeFlowSpeed,FreeFlowSpeed_Units)
    (TravelTime, unit) = TravelTime_BPR_via_Flow(Flow, Flow_Units, FreeFlowSpeed, FreeFlowSpeed_Units, 
                                                 LinkLength, LinkLength_Units, Method='old')
    
v0.1    June 6, 2019    initial release

_author = 	"Sinan Salman (sinan.salman@zu.ac.ae)"
_version = 	"Version: 0.1"
_date = 	"Date: 2019/06/06"
_copyright= "Copyright (c)2019 Sinan Salman"
_license = 	"GPLv3"
"""

import numpy as _np
import scipy.optimize as _sp_o
import scipy.interpolate as _sp_i

_ft_per_km = 3280.84  # 1km = 3280.84ft
_sec_per_h = 3600    # 1h = 3600sec

_Density = _np.concatenate([ _np.arange(0,100,1),_np.arange(100,250,5)]) # density veh/(km.ln)
_MaxDensity = _Density.max()
_SpeedAtMaxDensity = 1 # km/h
_Speed_min, _Speed_max = 1, 120  # parameters for setting up valid Speed range for root search

Polutants = ['CO','NOx','VOC']

# COPERT v5 units: Speed:km/h --> Emissions:g/(km.veh)
# Ntziachristos, L., & Samaras, Z. (2017). EMEP/EEA air pollutant emission inventory guidebook 2016. Retrieved from https://www.eea.europa.eu/publications/emep-eea-guidebook-2016
_COPERT5_coef  = {  'CO' :{ 'Alpha': 0.000445265807540281, 
                            'Beta':-0.102075931519888,    
                            'Gamma':6.87692822338964,   
                            'Delta':10.3838494701837,   
                            'Epsilon':0.00162104796148859,  
                            'Zeta':-0.43756305596665,   
                            'Eta':30.3373327734632,  
                            'RF':0},
                    'NOx':{ 'Alpha':-0.000314540672240406, 
                            'Beta': 0.103056777252331,    
                            'Gamma':0.239056575951972,  
                            'Delta':-0.3392791490405,   
                            'Epsilon':0.0345358283618128,   
                            'Zeta': 1.98601274556741,   
                            'Eta': 1.26376326022528, 
                            'RF':0},
                    'VOC':{ 'Alpha': 3.81411947410715E-06, 
                            'Beta':-0.000707365446956988,
                            'Gamma':0.0452491977889404, 
                            'Delta': 0.173073521070119, 
                            'Epsilon':0.000069899000382836, 
                            'Zeta':-0.0475381691356608, 
                            'Eta': 6.21205348188937, 
                            'RF':0}}

# TRANSYT-7F units: Speed:ft/sec --> Emissions:g/ft.veh  (A:g/(ft.veh) B:sec/ft C:sec/ft)
# Penic, M. A., & Upchurch, J. (1992). TRANSYT-7F: enhancement for fuel consumption, pollution emissions, and user costs. Transportation Research Record, (1360).
_TRANSYT7f_coef = { 'CO' :{ 'A':3.3963, 'B':0.014561, 'C':1000},
                    'NOx':{ 'A':1.5718, 'B':0.040732, 'C':10000},
                    'VOC':{ 'A':2.7843, 'B':0.015062, 'C':10000}}

# medean air pollution external costs $/g (in 1992)
# Matthews, H. S. (n.d.). The External Costs of Air Pollution and the Environmental Impact of the Consumer in the U.S. Economy. 203.
Costs = {   'CO' : 0.00052,
            'NOx': 0.00106,
            'VOC': 0.00140}
TimeAdjustedCosts = {k:v*1.7898 for k,v in Costs.items()} # adjusted deom 1992 to 2018 using CPI change of 178.98%; https://www.officialdata.org/us/inflation/1992?amount=100
# look info using GDP-deflator as Matthews did

# BASED ON ABU DHABI ROAD NETWORK ASSUMPTIONS
# FFS is assumed to equal MaxSpeed
# 90km/h < Speed < 120km/h assumed to be freeway
# 60km/h < Speed <  90km/h assumed to be multilane highway
_BPR_coef = {120.0:{'c': 2400, 'a':0.39, 'b':6.3, 'type':'freeway'}, # C is based on Exhibts 21-2 & 23-2
             100.0:{'c': 2300, 'a':0.25, 'b':9.0, 'type':'freeway'}, # minimal difference with 104km/h
              80.0:{'c': 2000, 'a':0.07, 'b':6.0, 'type':'Multilane HW'},
              60.0:{'c': 1800, 'a':0.06, 'b':6.0, 'type':'Multilane HW'}} # a and b estimated using regression line


def COPERT5(Speed, Speed_Units, Polutant):
    """Estimate emissions using COPERT v5
    (EmissionsFactor, unit) = TRANSYT7f(Speed, Speed_Units, Polutant)
        Polutants include: CO, NOx, and VOC
        Speed must be in km/h
    results will be in g/(km.veh)"""
    c = _COPERT5_coef.get(Polutant)
    if c == None:
        raise ValueError(f'Unrecognized polutant:{Polutant}')
    if Speed_Units != 'km/h':
        raise ValueError(f'Speed must be provided in km/h units')
    return \
        ( c['Alpha']*Speed**2 + c['Beta']*Speed + c['Gamma'] + c['Delta']/Speed ) / \
        ( c['Epsilon']*Speed**2 + c['Zeta']*Speed + c['Eta']) * (1-c['RF']) , \
        'g/(km.veh)'


def TRANSYT7f(Speed, Speed_Units, Polutant):
    """Estimate emissions using TRANSYT-7F
    (EmissionsFactor, unit) = TRANSYT7f(Speed, Speed_Units, Polutant)
        Polutants include: CO, NOx, and VOC
        Speed must be in km/h
    results will be in g/(km.veh)"""
    c = _TRANSYT7f_coef.get(Polutant)
    if c == None:
        raise ValueError(f'Unrecognized polutant:{Polutant}')
    if Speed_Units != 'km/h':
        raise ValueError(f'Speed must be provided in km/h units')
    Speed_ft_per_sec = Speed*_ft_per_km/_sec_per_h
    return (c['A'] * _np.exp(c['B']*Speed_ft_per_sec) / (c['C']*Speed_ft_per_sec)) * _ft_per_km , \
           'g/(km.veh)'


def ExternalCost(Speed, Speed_Units, Method = 'TRANSYT7f'):
    """Estimate emissions external cost using (Matthews 1999)
    (cost, unit) = ExternalCost(Speed, Speed_Units, Method = 'TRANSYT7f')
        Speed must be in km/h
    results will be in $/(km.veh)"""
    if Method.upper() == 'TRANSYT7F':
        EF = TRANSYT7f
    elif Method.upper() == 'COPERT5':
        EF = COPERT5
    else:
        raise ValueError(f'Unrecognized emissions estimation method: {Method}')
    if Speed_Units != 'km/h':
        raise ValueError(f'Speed must be provided in km/h units')
    return sum([EF(Speed, Speed_Units, p)[0] * TimeAdjustedCosts[p] for p in Polutants]), '$/(km.veh)'


def TravelTime_BPR_via_Flow(Flow, Flow_Units, 
                            FreeFlowSpeed, FreeFlowSpeed_Units, 
                            LinkLength, LinkLength_Units, 
                            Method='old'):
    """Estimate link travel time using Bureau of Public Roads (BPR) curve
    (TravelTime, unit) = TravelTime_BPR_via_Flow(Flow, Flow_Units, 
                                                 FreeFlowSpeed, FreeFlowSpeed_Units, 
                                                 LinkLength, LinkLength_Units, 
                                                 Method='old')
        Flow in veh/(h.lane)
        FreeFlowSpeed is average speed in free flow conditions (km/h)
        LinkLength in km
        Method can be: 'old' for the original BPR 1960's parameters (a=0.15, b=4), or 
                       'new' for Abu Dhabi assumptions
    results will be in hours"""
    data = _BPR_coef.get(FreeFlowSpeed)
    if data == None:
        raise ValueError(f'Unsupported free flow speed:{FreeFlowSpeed}')
    if Method == 'old':
        a, b, c = 0.15, 4, data['c']
    elif Method == 'new':
        a, b, c = data['a'], data['b'], data['c']
    else:
        raise ValueError(f'Unsupported calculation method:{Method}')
    return LinkLength/FreeFlowSpeed * ( 1 + a * ( Flow / c )**b ) , 'h'


def _Speed_BPR_via_Density(Speed, FreeFlowSpeed, Density, Method='old'):
    """This function should not be used directly; it is used internally in finding the root via scipy.optimize.brentq().
    Estimate link speed using Bureau of Public Roads (BPR) curve, however basecd on Density instead of Flow
    TravelTime_Polynomial_Evaluation = _Speed_BPR_via_Density(Speed, FreeFlowSpeed, Density, Method='old')
    Speed is the average speed at 'Density' conditions (km/h)
    FreeFlowSpeed is average speed in free flow conditions (km/h)
    Density in veh/(km.lane)
    Method can be: 'old' for the original BPR 1960's parameters (a=0.15, b=4), or 
                   'new' for Abu Dhabi assumptions"""
    data = _BPR_coef.get(FreeFlowSpeed)
    if data == None:
        raise ValueError(f'Unsupported free flow speed:{FreeFlowSpeed}')
    if Method == 'old':
        a, b, c = 0.15, 4, data['c']
    elif Method == 'new':
        a, b, c = data['a'], data['b'], data['c']
    else:
        raise ValueError(f'Unsupported calculation method:{Method}')
    Dc = c*(1+a)/FreeFlowSpeed
    return 1/Speed**(b+1) - (1/FreeFlowSpeed)*(1/Speed**b) - (a * (1+a)**b * (Density/Dc)**b * 1/FreeFlowSpeed**(b+1))


_Speed_functions = {}  # place holder for TravelTime(via Density) interpolation functions initialized later

def Initialize_Speed_Interpolation_Functions(Method='old'):
    """Setup Density_to_Speed hash tables for interpolation and create interpolation classes used in Speed_via_Density()
    Method can be: 'old' for the original BPR 1960's parameters (a=0.15, b=4), or 
                   'new' for Abu Dhabi assumptions"""
    for ffs in _BPR_coef.keys():
        Speed = [_sp_o.brentq(_Speed_BPR_via_Density,a=_Speed_min,b=_Speed_max,args=(ffs,d,Method)) for d in _Density]
        _Speed_functions[ffs] = _sp_i.interp1d(_Density,Speed,kind='quadratic')


def Speed_via_Density(Density,Density_Units,FreeFlowSpeed,FreeFlowSpeed_Units):
    """Estimate link speed using Bureau of Public Roads (BPR) curve, however basecd on Density instead of Flow
    (speed, unit) = Speed_via_Density(Density,Density_Units,FreeFlowSpeed,FreeFlowSpeed_Units)
        Density in veh/(km.lane)
        FreeFlowSpeed is average speed in free flow conditions (km/h)
    Results will be in (km/h)"""
    if Density_Units != 'veh/(km.lane)':
        raise ValueError(f'Speed must be provided in veh/(km.lane) units')
    if FreeFlowSpeed_Units != 'km/h':
        raise ValueError(f'Speed must be provided in km/h units')
    if isinstance(Density, list):
        Density=_np.array(Density)
    if isinstance(FreeFlowSpeed, list):
        FreeFlowSpeed=_np.array(FreeFlowSpeed)
    if not isinstance(Density, _np.ndarray):
        Density=_np.array([Density])
    if not isinstance(FreeFlowSpeed, _np.ndarray):
        FreeFlowSpeed=_np.array([FreeFlowSpeed])
    # assert Density.shape==FreeFlowSpeed.shape,f'Density and FreeFlowSpeed have different lengths: {Density.shape} != {FreeFlowSpeed.shape}'
    return _np.array([_Speed_functions[ffs](D) if D<_MaxDensity else _SpeedAtMaxDensity for ffs, D in zip(FreeFlowSpeed,Density)]), 'km/h'


if __name__ == '__main__':
    print('emissions.py is designed to operate as module')
else:
    Initialize_Speed_Interpolation_Functions()
