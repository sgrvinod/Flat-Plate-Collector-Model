#!/usr/bin/python
import numpy as np
import scipy as sci
from scipy.interpolate import interp1d
from scipy.optimize import brentq

def waterprops(T):
    import numpy as np
    import scipy as sci
    #Returns k, rho, cp, mu, nu, alpha in SI units

    waterdata = np.genfromtxt('WaterAtPressure100kPa.txt', delimiter = '\t')

    thermalcondfn = interp1d(waterdata[:,0],waterdata[:,12])
    dynviscosityfn = interp1d(waterdata[:,0],waterdata[:,11])
    densityfn = interp1d(waterdata[:,0],waterdata[:,2])
    SpHeatCpfn = interp1d(waterdata[:,0],waterdata[:,8])

    return thermalcondfn(T), densityfn(T), 1000*SpHeatCpfn(T), \
           dynviscosityfn(T), dynviscosityfn(T)/densityfn(T), \
           thermalcondfn(T)/(densityfn(T)*1000*SpHeatCpfn(T))

def airprops(T):
    import numpy as np
    import scipy as sci
    #Returns k, rho, cp, mu, nu, alpha in SI units

    airdata = np.genfromtxt('AirPropertiesAtPressure100kPa.txt', delimiter = '\t')

    thermalcondfn = interp1d(airdata[:,0],airdata[:,3]*0.01)
    dynviscosityfn = interp1d(airdata[:,0],airdata[:,4]*1e-5)
    densityfn = interp1d(airdata[:,0],airdata[:,1])
    SpHeatCpfn = interp1d(airdata[:,0],airdata[:,2])

    return thermalcondfn(T), densityfn(T), 1000*SpHeatCpfn(T), \
           dynviscosityfn(T), dynviscosityfn(T)/densityfn(T), \
           thermalcondfn(T)/(densityfn(T)*1000*SpHeatCpfn(T))
           
