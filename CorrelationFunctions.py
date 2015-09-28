#!/usr/bin/python
import numpy as np
import scipy as sci
from scipy.interpolate import interp1d

sb = 5.67e-8 #Stefan Boltzmann constant

def RayleighNum(grav, CTE, deltaT, lengthscale, nu, alpha):
    return grav*CTE*deltaT*lengthscale**3/(nu*alpha)

def EffectiveRayleighNum(RaNum, inclin):
    #inclin is inclination to horizontal in radians
    return RaNum*np.cos(inclin)

def NuForInclinedParallelSurfs(EffRaNumber):
    import numpy as np
    import scipy as sci

    if EffRaNumber <= 1708:
        NuL = 1.0
    elif 1708 < EffRaNumber <= 5900:
        NuL = 1 + 1.446*(1-1708.0/EffRaNumber)
    elif 5900 < EffRaNumber <= 9.23e4:
        NuL = 0.229*(EffRaNumber**0.252)
    elif 9.23e4 < EffRaNumber <= 1.0e6:
        NuL = 0.157*(EffRaNumber**0.285)

    return NuL

def HTCTopCover(windvelocity):
    #windvelocity is speed in m/s
    #returns HTC in SI units
    return 8.55 + 2.56*windvelocity

def radtranspp(T1, T2, eps1, eps2):
    #radiative transfer between two parallel plates given their temperatures and
    #emissivities
    return sb*(T1**4-T2**4)/(1.0/eps1+1.0/eps2-1)

def radtransObjAmb(T1, T2, eps1):
    #radiative transfer between two parallel plates given their temperatures and
    #emissivities
    return sb*eps1*(T1**4-T2**4)

def convtrans(T1, T2, htc):
    #convective transfer per unit area between two objects given their 
    #temperatures and heat transfer coefficient
    return htc*(T1-T2)
    
