#!/usr/bin/python
import numpy as np
import scipy as sci
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rc
from scipy.interpolate import interp1d
from scipy.optimize import brentq
#from pylab import figure, plot, xlabel, grid, hold, legend, title, savefig
from pylab import *
import FluidProperties as fp
import CorrelationFunctions as cf
import FlatPlateCollectorFunctions as FPCfns

print "Water:"
print "Thermal conductivity is ", fp.waterprops(300.0)[0]
print "Dynamic viscosity is ", fp.waterprops(300.0)[3]
print "Density is ", fp.waterprops(300.0)[1]
print "Kinematic viscosity is ", fp.waterprops(300.0)[4]
print "Thermal diffusivity is ", fp.waterprops(300.0)[5]

print "Air: "
print "Thermal conductivity is ", fp.airprops(300.0)[0]
print "Dynamic viscosity is ", fp.airprops(300.0)[3]
print "Density is ", fp.airprops(300.0)[1]
print "Kinematic viscosity is ", fp.airprops(300.0)[4]
print "Thermal diffusivity is ", fp.airprops(300.0)[5]
print "\n"


AbsPlateL1 = 1.9
AbsPlateL2 = 1.0
AbsCovGlassL = 0.04
AbsPlateEmiss = 0.92
CovGlassPlateEmiss = 0.88
CollectorTilt = 20.0 # in degrees
CollectorTilt = CollectorTilt*np.pi/180.0
AbsPlateMeanTemp = 70.0 # degree centigrade
AmbAirTemperature = 24.0 # degree centigrade
WindSpeed = 2.5 # in m/s
BackInsulThick = 0.08
SideInsulThick = 0.04
kInsul = 0.05 # W/m/K

AbsPlateL3 = 2*AbsCovGlassL + BackInsulThick

#Initialise count, heat flows and temperatures
count=0
Qin1=20.0
Qout1=10.0
Qin2= 10.0
Qout2=20.0
Qinp=200.0*AbsPlateL1*AbsPlateL2 
Qoutp=30.0
tol = 0.01
tc1 = 40.0
tc2=30.0

#Assume the energy absorbed by the absorber plate is convective heat transfer from a dummy temperature
tdummy=200

#Calculate cover glass and absorber plate temperatures
while ((np.abs((Qin1-Qout1)/(0.5*(Qin1+Qout1)))>tol) or (np.abs((Qin2-Qout2)/(0.5*(Qin2+Qout2)))>tol) or (np.abs((Qinp-Qoutp)/(0.5*(Qinp+Qoutp)))>tol)):
    
    count+=1
    
    #Calculate h-conv due to wind
    hwind = cf.HTCTopCover(WindSpeed)
    
    #Calculate sky temperature
    tsky = AmbAirTemperature - 6
    
    #Calculate mean temperatures for calculating h-conv between 1. plate and tc1, and 2. tc1 and tc2
    Mean1 = (AbsPlateMeanTemp+tc1)/2
    Mean2 = (tc2+tc1)/2
    
    #Calculate air properties for both interfaces (between 1. plate and tc1, and 2. tc1 and tc2)
    airproperties1 = fp.airprops(Mean1+273)
    airproperties2 = fp.airprops(Mean2+273)
    
    #Calculate AirCTE for both interfaces (between 1. plate and tc1, and 2. tc1 and tc2)
    AirCTE1 = 1.0/(Mean1+273)
    AirCTE2 = 1.0/(Mean2+273)
    
    #Calculate Rayleigh Numbers for both interfaces (between 1. plate and tc1, and 2. tc1 and tc2)
    Ra1 = cf.RayleighNum(9.8, AirCTE1, np.abs(AbsPlateMeanTemp-Mean1),\
                        AbsCovGlassL,airproperties1[4],airproperties1[5])
    Ra2 = cf.RayleighNum(9.8, AirCTE2, np.abs(tc1-Mean2),\
                        AbsCovGlassL,airproperties2[4],airproperties2[5])
    
    #Calculate effective Rayleigh Numbers for both interfaces (between 1. plate and tc1, and 2. tc1 and tc2)
    EffRa1 = cf.EffectiveRayleighNum(Ra1, CollectorTilt)
    EffRa2 = cf.EffectiveRayleighNum(Ra2, CollectorTilt)
    
    #Calculate Nusselt Number for both interfaces (between 1. plate and tc1, and 2. tc1 and tc2)
    NussL1 = cf.NuForInclinedParallelSurfs(EffRa1)
    NussL2 = cf.NuForInclinedParallelSurfs(EffRa2)
    
    #Calculate h-conv for both interfaces (between 1. plate and tc1, and 2. tc1 and tc2)
    h1 = NussL1*airproperties1[0]/AbsCovGlassL
    h2 = NussL2*airproperties2[0]/AbsCovGlassL
    
    #Calculate individual heat flows
    conv_plate_cg1=cf.convtrans(AbsPlateMeanTemp,tc1,h1)*AbsPlateL1*AbsPlateL2 
    rad_plate_cg1=cf.radtranspp(AbsPlateMeanTemp+273,tc1+273,AbsPlateEmiss,CovGlassPlateEmiss)*AbsPlateL1*AbsPlateL2 
    
    conv_cg1_cg2=cf.convtrans(tc1,tc2,h2)*AbsPlateL1*AbsPlateL2 
    rad_cg1_cg2=cf.radtranspp(tc1+273,tc2+273,CovGlassPlateEmiss,CovGlassPlateEmiss)*AbsPlateL1*AbsPlateL2 
    
    conv_cg2_amb=cf.convtrans(tc2,AmbAirTemperature,hwind)*AbsPlateL1*AbsPlateL2 
    rad_cg2_sky=cf.radtransObjAmb(tc2+273,tsky+273,CovGlassPlateEmiss)*AbsPlateL1*AbsPlateL2 
    
    #Calculate power input and output for cover glass 1    
    Qin1 = conv_plate_cg1 + rad_plate_cg1
    Qout1 = conv_cg1_cg2 + rad_cg1_cg2
    
    #Calculate power input and output for cover glass 2    
    Qin2 = conv_cg1_cg2 + rad_cg1_cg2
    Qout2 = conv_cg2_amb + rad_cg2_sky
    
    #Calculate Back and Side and Top loss from absorber plate
    BackLoss = FPCfns.FlatPlateCollectorBackLoss(kInsul,AbsPlateL1,AbsPlateL2,BackInsulThick,AbsPlateMeanTemp,AmbAirTemperature)
    SideLoss = FPCfns.FlatPlateCollectorSideLoss(kInsul,AbsPlateL1,AbsPlateL2,AbsPlateL3,SideInsulThick,AbsPlateMeanTemp,AmbAirTemperature)
    TopLoss = conv_plate_cg1 + rad_plate_cg1
    
    #Calculate power output from plate
    Qoutp = TopLoss + BackLoss + SideLoss
    
    #Calculate effective h-conv, assuming each power input/output to be purely convective
    effhQin1 = Qin1/(AbsPlateMeanTemp-tc1)
    effhQout1 = Qout1/(tc1-tc2)
    effhQin2 = Qin2/(tc1-tc2)
    effhQout2 = Qout2/(tc2-AmbAirTemperature)
    effhQinp= Qinp/(tdummy-AbsPlateMeanTemp)
    effhQoutp=Qoutp/(AbsPlateMeanTemp-tc1)
    
    #Update tc1 and tc2
    tc1 = (effhQin1*AbsPlateMeanTemp+effhQout1*tc2)/(effhQin1+effhQout1)
    tc2 = (effhQin2*tc1+effhQout2*AmbAirTemperature)/(effhQin2+effhQout2)
    
    #Update plate temperature
    AbsPlateMeanTemp = (effhQinp*tdummy+effhQoutp*tc1)/(effhQinp+effhQoutp)
    
    #Print values
    print "Iteration Number: %d" % count
    print "Assumed absorption plate temperature: %f C" % AbsPlateMeanTemp
    print "Assumed lower cover glass temperature: %f C" % tc1
    print "Assumed upper cover glass temperature: %f C" % tc2
    print "\n"
    
#Print final values
print "INDIVIDUAL HEAT FLOWS:"
print "Convection from plate to lower cover glass: %f W" % conv_plate_cg1
print "Radiation from plate to lower cover glass: %f W" % rad_plate_cg1
print "Convection from lower to upper cover glass: %f W" % conv_cg1_cg2
print "Radiation from lower to upper cover glass: %f W" % rad_cg1_cg2
print "Convection from upper cover glass to ambient: %f W" % conv_cg2_amb
print "Radiation from upper cover glass to sky: %f w" % rad_cg2_sky
print "\nCOMPONENT HEAT FLOWS:"
print "Heat flow into plate: %f W" % Qinp
print "Top loss from plate: %f W"% TopLoss
print "Back loss: %f W" % BackLoss
print "Side loss: %f W" % SideLoss
print "Therefore, total heat flow out of plate: %f W" % Qoutp
print "Heat flow into lower cover glass: %f W" % Qin1
print "Heat flow out of lower cover glass: %f W" % Qout1
print "Heat flow into upper cover glass: %f W" % Qin2
print "Heat flow out of upper cover glass: %f W" % Qout2
print "\nFINAL TEMPERATURES:"                                                         
print "Number of iterations: %d" % count
print "Absorption plate temperature = %f C" % AbsPlateMeanTemp
print "Lower cover glass temperature: %f C" % tc1
print "Upper cover glass temperature: %f C" % tc2

