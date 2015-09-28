#!/usr/bin/python
import numpy as np
import scipy as sci
from scipy.interpolate import interp1d

def FlatPlateCollectorBackLoss(kins, L1,L2, backinsthick,Tplate,Tamb):
    return kins*L1*L2*(Tplate-Tamb)/backinsthick

def FlatPlateCollectorSideLoss(kins, L1,L2,L3, sideinsthick,Tplate,Tamb):
    return kins*(2*(L1+L2)*L3)*(Tplate-Tamb)/(2*sideinsthick)

#def FlatPlateSingleCoverglassTopLoss
