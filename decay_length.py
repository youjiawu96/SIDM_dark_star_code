import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

hbar = 6.582*10**(-25) #GeV*s
me = 0.511*10**(-3) #GeV
mmu = 0.1056 #GeV
alpha = 1.0/137.0 #fine structure constant
mpi = 0.135 #GeV
mK = 0.495 #GeV
c = 3*10**8 #m/s
AU = 1.496* 10**11 #m
mu = 1.8*10**(-3) #GeV
md = 4.3*10**(-3) #GeV
ms = 92*10**(-3) #GeV

def Gamma_lepton(mA, ml, gy):
    # all masses in GeV
    if mA > 2*ml:
        return (2.94*10**(-6)*ml/me*gy)**2/8.0/np.pi * mA * (1 - 4*ml**2/mA**2)**0.5 / hbar
    else:
        return 0

def F(tau):
    if tau < 1.0:
        return 2/tau*(np.arcsin(tau**0.5))**2.0 + 0.0j
    else:
        return 2/tau*(-0.25*(np.log((1+(1-1/tau)**0.5)/(1-(1-1/tau)**0.5))-3.1415926j)**2)

def Gamma_photon(mA, gy):
    x = 3*4.0/9.0*2.94*10**(-6)*gy/mu*F(mA**2/4/mpi**2) + 3*1.0/9.0*2.94*10**(-6)*gy/md*F(mA**2/4/mpi**2) + 3*4.0/9.0*2.94*10**(-6)*gy/ms*F(mA**2/4/mK**2)
    return alpha**2.0*mA**3/(256*np.pi**3) * (x.real**2 + x.imag**2) / hbar

def decay_length(mA, mD, gy):
    # decay length in AU
    return mD/mA * c / (Gamma_lepton(mA, me, gy)+Gamma_lepton(mA, mmu, gy) + Gamma_photon(mA, gy))/AU


def eqn(mA, mD, gy):
    return decay_length(mA, mD, gy) - 1


D_mass = np.arange(-2.4, 4, 0.05)
D_mass = 10**D_mass

exclude1 = np.zeros(len(D_mass))
exclude2 = np.zeros(len(D_mass))

for i in range(len(D_mass)):
    lim1 = fsolve(eqn, 0.001, args=(D_mass[i], 8*10**(-5)))
    lim2 = fsolve(eqn, 0.001, args=(D_mass[i], 10**(-5)))
    exclude1[i] = lim1[0]
    exclude2[i] = lim2[0]

data = np.vstack((D_mass, exclude1, exclude2))
data = np.transpose(data)
print(data)

np.savetxt("scalar_exclude_lim.txt", data, delimiter=' ')