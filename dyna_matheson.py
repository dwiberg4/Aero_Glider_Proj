import numpy as np
from scipy.linalg import eig
import json


input = 'BaselineGlider.json'

def printEigen(eigval,eigvec,mode,lat,long,dimensional):
    """Prints an Eigenvalue Summary to the terminal with data depending on the information given to the Eigenvalue
    The program must be given an eigenvalue to be able to do any printing. The eigenvector is only important to print
    the amplitude and phase table. The mode is determined by the mode string variable. The Lat and Long variables determined
    which variables to use when printing the table of amplitude and phase tables. The final variable dimensional is whether to
    multiply 2Vo/c or 2Vo/b to get final results.

    Args:
        eigval (np complex): Complex number that represents the eigen value of the system.
        eigvec (1D np array dtype = complex): Vector of complex numberst that represent the eigenvector that corresponds to the value in eigval.
        mode (str): A string that gives the name of mode this eigval corresponds to. 
        lat (Bool): If true, the lateral amplitude and phase table will be printed. 
        long (Bool): if True, the longitudinal amplitude and phase table will be printed. Requires an eigvector.
        dimensional (float): A float to represent the value to dimensionalize with. 
    """
    print('The {0:<20s} is given by {1:9.7f} ± {2:9.7f} i.'.format(mode,eigval.real,np.abs(eigval.imag)))
    sig = -eigval.real*dimensional
    damped = (0>eigval.real)
    damp_nat_freq = np.abs(eigval.imag)*dimensional
    damp_rate = -eigval.real

    if damp_nat_freq > 0.000001:
        # If there is an imaginary part calculate and print all of these values.
        eig1 = np.complex(sig,damp_nat_freq)
        eig2 = np.complex(sig,-damp_nat_freq)
        undamp_freq = np.sqrt(eig1*eig2)*dimensional
        period = 2*np.pi/damp_nat_freq
        damp_ratio = -(eig1 + eig2)/(2*np.sqrt(eig1*eig2))
        if lat:
            print(''.ljust(56,'='))
            print('{0:<15s}{1:>15s}{2:>10s}'.format('Component','Amplitude','Phase'))
            print('{0:<20s}{1:>9.7f}{2:>10.2f}°'.format('Δμ',np.abs(eigvec[0]),np.arctan2(eigvec[0].imag,eigvec[0].real)))
            print('{0:<20s}{1:>9.7f}{2:>10.2f}°'.format('Δα',np.abs(eigvec[1]),np.arctan2(eigvec[1].imag,eigvec[1].real)))
            print('{0:<20s}{1:>9.7f}{2:>10.2f}°'.format('Δqbar',np.abs(eigvec[2]),np.arctan2(eigvec[2].imag,eigvec[2].real)))
            print('{0:<20s}{1:>9.7f}{2:>10.2f}°'.format('Δξx',np.abs(eigvec[3]),np.arctan2(eigvec[3].imag,eigvec[3].real)))
            print('{0:<20s}{1:>9.7f}{2:>10.2f}°'.format('Δξz',np.abs(eigvec[4]),np.arctan2(eigvec[4].imag,eigvec[4].real)))
            print('{0:<20s}{1:>9.7f}{2:>10.2f}°'.format('Δθ',np.abs(eigvec[5]),np.arctan2(eigvec[5].imag,eigvec[5].real)))
            print(''.ljust(56,'='))
        if long:
            print(''.ljust(56,'='))        
            print('{0:<15s}{1:>15s}{2:>10s}'.format('Component','Amplitude','Phase'))
            print('{0:<20s}{1:>9.7f}{2:>10.2f}°'.format('Δβ',np.abs(eigvec[0]),np.arctan2(eigvec[0].imag,eigvec[0].real)))
            print('{0:<20s}{1:>9.7f}{2:>10.2f}°'.format('Δpbar',np.abs(eigvec[1]),np.arctan2(eigvec[1].imag,eigvec[1].real)))
            print('{0:<20s}{1:>9.7f}{2:>10.2f}°'.format('Δrbar',np.abs(eigvec[2]),np.arctan2(eigvec[2].imag,eigvec[2].real)))
            print('{0:<20s}{1:>9.7f}{2:>10.2f}°'.format('Δξy',np.abs(eigvec[3]),np.arctan2(eigvec[3].imag,eigvec[3].real)))
            print('{0:<20s}{1:>9.7f}{2:>10.2f}°'.format('Δφ',np.abs(eigvec[4]),np.arctan2(eigvec[4].imag,eigvec[4].real)))
            print('{0:<20s}{1:>9.7f}{2:>10.2f}°'.format('Δψ',np.abs(eigvec[5]),np.arctan2(eigvec[5].imag,eigvec[5].real)))
            print(''.ljust(56,'='))
        print('The {1:<20s} Period is {0:9.7f}'.format(period,mode))   
        print('The {1:<20s} Damped Frequency is {0:9.7f}'.format(damp_nat_freq,mode))                   
    else:
        # If there is not an imaginary part then do not print phase in the table.
        if lat:
            print(''.ljust(56,'='))
            print('{0:<15s}{1:>15s}'.format('Component','Amplitude'))
            print('{0:<20s}{1:>9.7f}'.format('Δμ',np.abs(eigvec[0])))
            print('{0:<20s}{1:>9.7f}'.format('Δα',np.abs(eigvec[1])))
            print('{0:<20s}{1:>9.7f}'.format('Δqbar',np.abs(eigvec[2])))
            print('{0:<20s}{1:>9.7f}'.format('Δξx',np.abs(eigvec[3])))
            print('{0:<20s}{1:>9.7f}'.format('Δξz',np.abs(eigvec[4])))
            print('{0:<20s}{1:>9.7f}'.format('Δθ',np.abs(eigvec[5])))
            print(''.ljust(56,'='))
        if long:
            print(''.ljust(56,'='))        
            print('{0:<15s}{1:>15s}'.format('Component','Amplitude'))
            print('{0:<20s}{1:>9.7f}'.format('Δβ',np.abs(eigvec[0])))
            print('{0:<20s}{1:>9.7f}'.format('Δpbar',np.abs(eigvec[1])))
            print('{0:<20s}{1:>9.7f}'.format('Δrbar',np.abs(eigvec[2])))
            print('{0:<20s}{1:>9.7f}'.format('Δξy',np.abs(eigvec[3])))
            print('{0:<20s}{1:>9.7f}'.format('Δφ',np.abs(eigvec[4])))
            print('{0:<20s}{1:>9.7f}'.format('Δψ',np.abs(eigvec[5])))
            print(''.ljust(56,'='))
    if damped:
        # If the system is damped, print a 99% Damping Time.
        damp_time = -np.log(0.01)/(sig)          # 99% damping time
        print('The {0:<20s} 99% Damping Time is {1:9.7f}'.format(mode,damp_time))
    else:
        # If the system is not damped, print a Doubling Time instead.
        damp_time = -np.log(2)/sig             # Doubling Time
        print('The {0:<20s} Doubling Time is {1:9.7f}'.format(mode,damp_time))
    print('The {1:<20s} Damping Rate is {0:9.7f}'.format(sig,mode))
    print()
    print()

    




with open(input,'r') as f:
    aircraft = json.load(f)
f.close()

## Read in all values from the JSON file.
Sw = aircraft["aircraft"]["wing_area[ft^2]"]
b = aircraft["aircraft"]["wing_span[ft]"]
W = aircraft["operating"]["weight[lbf]"]
gamma = (aircraft["operating"]["climb[deg]"]) * (np.pi/ 180)
rho = aircraft["operating"]["density[slugs/ft^3]"]
CLo = aircraft["reference"]["CL"]
Ixx = aircraft["reference"]["Ixx[slugs*ft^2]"]
Iyy = aircraft["reference"]["Iyy[slugs*ft^2]"]
Izz = aircraft["reference"]["Izz[slugs*ft^2]"]
Ixy = aircraft["reference"]["Ixy[slugs*ft^2]"]
Ixz = aircraft["reference"]["Ixz[slugs*ft^2]"]
Iyz = aircraft["reference"]["Iyz[slugs*ft^2]"]
hx = aircraft["reference"]["hx"]
hy = aircraft["reference"]["hy"]
hz = aircraft["reference"]["hz"]
CD0 = aircraft["reference"]["CD0"]
CD1 = aircraft["reference"]["CD1"]
CD2 = aircraft["reference"]["CD2"]
CL_a = aircraft["reference"]["CL,a"]
Cl_a = aircraft["reference"]["Cl,a"]
Cm_a = aircraft["reference"]["Cm,a"]
Cn_a = aircraft["reference"]["Cn,a"]
CL_ahat = aircraft["reference"]["CL,ahat"]
CD_ahat = aircraft["reference"]["CD,ahat"]
Cm_ahat = aircraft["reference"]["Cm,ahat"]
CL_uhat = aircraft["reference"]["CL,uhat"]
CD_uhat = aircraft["reference"]["CD,uhat"]
Cm_uhat = aircraft["reference"]["Cm,uhat"]
CY_b = aircraft["reference"]["CY,b"]
Cl_b = aircraft["reference"]["Cl,b"]
Cm_b = aircraft["reference"]["Cm,b"]
Cn_b = aircraft["reference"]["Cn,b"]
CY_p = aircraft["reference"]["CY,p"]
Cl_p = aircraft["reference"]["Cl,p"]
Cn_p = aircraft["reference"]["Cn,p"]
CL_q = aircraft["reference"]["CL,q"]
CD_q = aircraft["reference"]["CD,q"]
Cm_q = aircraft["reference"]["Cm,q"]
CY_r = aircraft["reference"]["CY,r"]
Cl_r = aircraft["reference"]["Cl,r"]
Cn_r = aircraft["reference"]["Cn,r"]

thrust_v = 0.0
alpha_0 = 0.0
theta = 0*(np.pi/180)
g = 32.17

c = Sw / b
Vo = np.sqrt((2*W*np.cos(theta)/(rho*Sw*CLo)))      # Equation 10.45 Ch. 7 Overview
CD_o = CD0 + (CD1*CLo) + ((CD2 * (CLo**2)))         # Equation 10.46 Ch. 7 Overview
CD_a = CD1*CL_a + 2*CD2*CLo*CL_a                    # Equation 10.57 Ch. 7 Overview
CT_V = thrust_v/ (0.5*rho*Vo*Sw)                    # Equation 10.74 Ch. 7 Overview
Cm_o = 0.0
z_To = 0.00

R_gx = g*c / (2*Vo**2)                               # Equation 10.71 Ch. 7 Overview
R_gy = g*b / (2*Vo**2)                               # Equation 10.71 Ch. 7 Overview
R_rhox = 4*W/g / (rho*Sw*c)                          # Equation 10.72 Ch. 7 Overview
R_rhoy = 4*W/g / (rho*Sw*b)                          # Equation 10.72 Ch. 7 Overview
R_xx = 8*Ixx / (rho*Sw*b**3)                         # Equation 10.73 Ch. 7 Overview
R_yy = 8*Iyy / (rho*Sw*c**3)                         # Equation 10.73 Ch. 7 Overview
R_zz = 8*Izz / (rho*Sw*b**3)                         # Equation 10.73 Ch. 7 Overview
R_xz = 8*Ixz / (rho*Sw*b**3)                         # Equation 10.73 Ch. 7 Overview

# Longitudinal Equations from 10.75 Ch. 7 Overview
A = np.zeros((6,6))
B = np.zeros((6,6))


A[0,0] = (-2*CD_o) + CT_V*np.cos(alpha_0)
A[1,0] = (-2*CLo) - (CT_V*np.sin(alpha_0))
A[2,0] = 2*Cm_o + CT_V*(z_To/c)
A[3,0] = np.cos(theta)
A[4,0] = -np.sin(theta)

A[0,1] = (CLo - CD_a)
A[1,1] = (-CL_a - CD_o)
A[2,1] = Cm_a
A[3,1] = np.sin(theta)
A[4,1] = np.cos(theta)

A[0,2] = - CD_q
A[1,2] = -CL_q + R_rhox
A[2,2] = Cm_q
A[5,2] = 1

A[0,5] = -R_rhox*R_gx*np.cos(theta)
A[1,5] = -R_rhox*R_gx*np.sin(theta)
A[3,5] = -np.sin(theta)
A[4,5] = -np.cos(theta)

B[0,0] = R_rhox + CD_uhat
B[1,0] = CL_uhat
B[2,0] = -Cm_uhat

B[0,1] = CD_ahat
B[1,1] = R_rhox + CL_ahat
B[2,1] = -Cm_ahat

B[2,2] = R_yy
B[3,3] = 1
B[4,4] = 1
B[5,5] = 1

C = np.matmul(np.linalg.inv(B),A)
long_eigvals, long_eigvecs = eig(C)


# Calculate the magnitude of each Eigenvalue
mag = np.abs(long_eigvals)

# Get index that sorts the magnitudes from largest to smallest.
idx = mag.argsort()[::-1]
#sort both values and vectors by that index.
long_eigvals = long_eigvals[idx]
long_eigvecs = long_eigvecs[:,idx]

## These are commented out and replaced by the printEigen function, but I left them for checking purposes.
# print('Short Period Mode is given by {0:9.7f} ± {1:9.7f} i.'.format(long_eigvals[0].real,np.abs(long_eigvals[0].imag)))
# print(''.ljust(56,'='))
# sp_sig = -(long_eigvals[0].real)*2*Vo/c
# sp_99damp = np.log(0.01)/-sp_sig
# sp_dampnatfreq = np.abs(long_eigvals[0].imag)*2*Vo/c
# sp_dampratio = -(long_eigvals[0] +long_eigvals[1])/(2*np.sqrt(long_eigvals[0]*long_eigvals[1]))
# sp_undampnatfreq = 2*Vo/c*np.sqrt(long_eigvals[0]*long_eigvals[1])
# sp_period = 2*np.pi/sp_dampnatfreq
# print('{0:<15s}{1:>15s}{2:>10s}'.format('Component','Amplitude','Phase'))
# print('{0:<20s}{1:>9.7f}{2:>10.2f}°'.format('Δμ',np.abs(long_eigvecs[0,0]),np.arctan2(long_eigvecs[0,0].imag,long_eigvecs[0,0].real)))
# print('{0:<20s}{1:>9.7f}{2:>10.2f}°'.format('Δα',np.abs(long_eigvecs[1,0]),np.arctan2(long_eigvecs[1,0].imag,long_eigvecs[1,0].real)))
# print('{0:<20s}{1:>9.7f}{2:>10.2f}°'.format('Δqbar',np.abs(long_eigvecs[2,0]),np.arctan2(long_eigvecs[2,0].imag,long_eigvecs[2,0].real)))
# print('{0:<20s}{1:>9.7f}{2:>10.2f}°'.format('Δξx',np.abs(long_eigvecs[3,0]),np.arctan2(long_eigvecs[3,0].imag,long_eigvecs[3,0].real)))
# print('{0:<20s}{1:>9.7f}{2:>10.2f}°'.format('Δξz',np.abs(long_eigvecs[4,0]),np.arctan2(long_eigvecs[4,0].imag,long_eigvecs[4,0].real)))
# print('{0:<20s}{1:>9.7f}{2:>10.2f}°'.format('Δθ',np.abs(long_eigvecs[5,0]),np.arctan2(long_eigvecs[5,0].imag,long_eigvecs[5,0].real)))
# print(''.ljust(56,'='))
# print('The Short Period Damping Rate is {0:9.7f}'.format(sp_sig))
# print('The Short Period 99% Damping Time is {0:9.7f}'.format(sp_99damp))
# print('The Short Period Damped Frequency is {0:9.7f}'.format(sp_dampnatfreq))
# print('The Short Period Period is {0:9.7f}'.format(sp_period))
# print()
# print()
   
# print('Long Period Mode is given by {0:9.7f} ± {1:9.7f} i.'.format(long_eigvals[2].real,np.abs(long_eigvals[2].imag)))
# print(''.ljust(56,'='))
# sp_sig = -(long_eigvals[2].real)*2*Vo/c
# sp_99damp = np.log(0.01)/-sp_sig
# sp_dampnatfreq = np.abs(long_eigvals[2].imag)*2*Vo/c
# sp_dampratio = -(long_eigvals[2] +long_eigvals[3])/(2*np.sqrt(long_eigvals[2]*long_eigvals[3]))
# sp_undampnatfreq = 2*Vo/c*np.sqrt(long_eigvals[2]*long_eigvals[3])
# sp_period = 2*np.pi/sp_dampnatfreq
# print('{0:<15s}{1:>15s}{2:>10s}'.format('Component','Amplitude','Phase'))
# print('{0:<20s}{1:>9.7f}{2:>10.2f}°'.format('Δμ',np.abs(long_eigvecs[0,2]),np.arctan2(long_eigvecs[0,2].imag,long_eigvecs[0,2].real)))
# print('{0:<20s}{1:>9.7f}{2:>10.2f}°'.format('Δα',np.abs(long_eigvecs[1,2]),np.arctan2(long_eigvecs[1,2].imag,long_eigvecs[1,2].real)))
# print('{0:<20s}{1:>9.7f}{2:>10.2f}°'.format('Δqbar',np.abs(long_eigvecs[2,2]),np.arctan2(long_eigvecs[2,2].imag,long_eigvecs[2,2].real)))
# print('{0:<20s}{1:>9.7f}{2:>10.2f}°'.format('Δξx',np.abs(long_eigvecs[3,2]),np.arctan2(long_eigvecs[3,2].imag,long_eigvecs[3,2].real)))
# print('{0:<20s}{1:>9.7f}{2:>10.2f}°'.format('Δξz',np.abs(long_eigvecs[4,2]),np.arctan2(long_eigvecs[4,2].imag,long_eigvecs[4,2].real)))
# print('{0:<20s}{1:>9.7f}{2:>10.2f}°'.format('Δθ',np.abs(long_eigvecs[5,2]),np.arctan2(long_eigvecs[5,2].imag,long_eigvecs[5,2].real)))
# print(''.ljust(56,'='))
# print('The Long Period Mode Damping Rate is {0:9.7f}'.format(sp_sig))
# print('The Long Period Mode 99% Damping Time is {0:9.7f}'.format(sp_99damp))
# print('The Long Period Mode Damped Frequency is {0:9.7f}'.format(sp_dampnatfreq))
# print('The Long Period Mode Period is {0:9.7f}'.format(sp_period))
# print()
# print()

# All of that commented code is contained in these two function calls.
printEigen(long_eigvals[0],long_eigvecs[:,0],'Short Period Mode',False,True,2*Vo/c)
printEigen(long_eigvals[2],long_eigvecs[:,2],'Long Period Mode',False,True,2*Vo/c)



# Lateral Equations from 10.76 Ch. 7 Overview
D = np.zeros((6,6))
E = np.zeros((6,6))

D[0,0] = CY_b
D[1,0] = Cl_b
D[2,0] = Cn_b
D[3,0] = 1

D[0,1] = CY_p
D[1,1] = Cl_p
D[2,1] = Cn_p
D[4,1] = 1

D[0,2] = CY_r - R_rhoy
D[1,2] = Cl_r
D[2,2] = Cn_r
D[4,2] = np.tan(theta)
D[5,2] = 1/np.cos(theta)

D[0,4] = R_rhoy * R_gy * np.cos(theta)

D[3,5] = np.cos(theta)

E[0,0] = R_rhoy

E[1,1] = R_xx
E[2,1] = -R_xz

E[1,2] = -R_xz
E[2,2] = R_zz

E[3,3] = 1
E[4,4] = 1
E[5,5] = 1

F = np.matmul(np.linalg.inv(E),D)
lat_eigvals, lat_eigvecs = eig(F)

# Calculate the magnitude of each Eigenvalue if there is an imaginary part give it a 1000 penalty
mag = np.where(lat_eigvals.imag == 0,np.abs(lat_eigvals),np.abs(lat_eigvals)-1000)


# Sort by lat_eigvals, lat_eigvecs by largest magnitude if there is no imaginary part.
idx = mag.argsort()[::-1]
lat_eigvals = lat_eigvals[idx]
lat_eigvecs = lat_eigvecs[:,idx]

# print('The Roll Mode is given by {0:9.7f} ± {1:9.7f} i.'.format(lat_eigvals[0].real,lat_eigvals[0].imag))
# print(''.ljust(56,'='))
# sp_sig = -(lat_eigvals[0].real)*2*Vo/c
# sp_99damp = np.log(0.01)/-sp_sig
# sp_dampnatfreq = np.abs(lat_eigvals[0].imag)*2*Vo/c
# sp_dampratio = -(lat_eigvals[0] +lat_eigvals[1])/(2*np.sqrt(lat_eigvals[0]*lat_eigvals[1]))
# sp_undampnatfreq = 2*Vo/c*np.sqrt(lat_eigvals[0]*lat_eigvals[1])
# sp_period = 2*np.pi/sp_dampnatfreq
# print('{0:<15s}{1:>15s}{2:>10s}'.format('Component','Amplitude','Phase'))
# print('{0:<20s}{1:>9.7f}{2:>10.2f}°'.format('Δβ',np.abs(lat_eigvecs[0,0]),np.arctan2(lat_eigvecs[0,0].imag,lat_eigvecs[0,0].real)))
# print('{0:<20s}{1:>9.7f}{2:>10.2f}°'.format('Δpbar',np.abs(lat_eigvecs[1,0]),np.arctan2(lat_eigvecs[1,0].imag,lat_eigvecs[1,0].real)))
# print('{0:<20s}{1:>9.7f}{2:>10.2f}°'.format('Δrbar',np.abs(lat_eigvecs[2,0]),np.arctan2(lat_eigvecs[2,0].imag,lat_eigvecs[2,0].real)))
# print('{0:<20s}{1:>9.7f}{2:>10.2f}°'.format('Δξy',np.abs(lat_eigvecs[3,0]),np.arctan2(lat_eigvecs[3,0].imag,lat_eigvecs[3,0].real)))
# print('{0:<20s}{1:>9.7f}{2:>10.2f}°'.format('Δφ',np.abs(lat_eigvecs[4,0]),np.arctan2(lat_eigvecs[4,0].imag,lat_eigvecs[4,0].real)))
# print('{0:<20s}{1:>9.7f}{2:>10.2f}°'.format('Δψ',np.abs(lat_eigvecs[5,0]),np.arctan2(lat_eigvecs[5,0].imag,lat_eigvecs[5,0].real)))
# print(''.ljust(56,'='))
# print('The Roll Mode Damping Rate is {0:9.7f}'.format(sp_sig))
# print('The Roll Mode 99% Damping Time is {0:9.7f}'.format(sp_99damp))
# #print('The Roll Mode Frequency is {0:9.7f}'.format(sp_undampnatfreq))
# #print('The Roll Mode Period is {0:9.7f}'.format(sp_period))
# print()
# print()
   
# print('Spiral Mode is given by {0:9.7f} ± {1:9.7f} i.'.format(lat_eigvals[1].real,lat_eigvals[1].imag))
# print(''.ljust(56,'='))
# sp_sig = -(lat_eigvals[1].real)*2*Vo/c
# sp_99damp = np.log(0.01)/-sp_sig
# sp_dampnatfreq = np.abs(lat_eigvals[2].imag)*2*Vo/c
# sp_dampratio = -(lat_eigvals[2] +lat_eigvals[3])/(2*np.sqrt(lat_eigvals[2]*lat_eigvals[3]))
# sp_undampnatfreq = 2*Vo/c*np.sqrt(lat_eigvals[2]*lat_eigvals[3])
# sm_double = np.log(2)/sp_sig
# sp_period = 2*np.pi/sp_dampnatfreq
# print('{0:<15s}{1:>15s}{2:>10s}'.format('Component','Amplitude','Phase'))
# print('{0:<20s}{1:>9.7f}{2:>10.2f}°'.format('Δβ',np.abs(lat_eigvecs[0,1]),np.arctan2(lat_eigvecs[0,1].imag,lat_eigvecs[0,1].real)))
# print('{0:<20s}{1:>9.7f}{2:>10.2f}°'.format('Δpbar',np.abs(lat_eigvecs[1,1]),np.arctan2(lat_eigvecs[1,1].imag,lat_eigvecs[1,1].real)))
# print('{0:<20s}{1:>9.7f}{2:>10.2f}°'.format('Δrbar',np.abs(lat_eigvecs[2,1]),np.arctan2(lat_eigvecs[2,1].imag,lat_eigvecs[2,1].real)))
# print('{0:<20s}{1:>9.7f}{2:>10.2f}°'.format('Δξy',np.abs(lat_eigvecs[3,1]),np.arctan2(lat_eigvecs[3,1].imag,lat_eigvecs[3,1].real)))
# print('{0:<20s}{1:>9.7f}{2:>10.2f}°'.format('Δφ',np.abs(lat_eigvecs[4,1]),np.arctan2(lat_eigvecs[4,1].imag,lat_eigvecs[4,1].real)))
# print('{0:<20s}{1:>9.7f}{2:>10.2f}°'.format('Δψ',np.abs(lat_eigvecs[5,1]),np.arctan2(lat_eigvecs[5,1].imag,lat_eigvecs[5,1].real)))
# print(''.ljust(56,'='))
# print('The Spiral Mode Damping Rate is {0:9.7f}'.format(sp_sig))
# print('The Spiral Mode 99% Damping Time is {0:9.7f}'.format(sp_99damp))
# print('The Spiral Mode doubling time is {0:9.7f}'.format(sp_dampnatfreq))
# #print('The Spiral Mode Period is {0:9.7f}'.format(sp_period))
# print()
# print()


# print('Dutch Roll Mode is given by {0:9.7f} ± {1:9.7f} i.'.format(lat_eigvals[4].real,lat_eigvals[4].imag))
# print(''.ljust(56,'='))
# sp_sig = -(lat_eigvals[4].real)*2*Vo/c
# sp_99damp = np.log(0.01)/-sp_sig
# sp_dampnatfreq = np.abs(lat_eigvals[4].imag)*2*Vo/c
# sp_dampratio = -(lat_eigvals[4] +lat_eigvals[5])/(2*np.sqrt(lat_eigvals[4]*lat_eigvals[5]))
# sp_undampnatfreq = 2*Vo/c*np.sqrt(lat_eigvals[4]*lat_eigvals[5])
# sp_period = 2*np.pi/sp_dampnatfreq
# print('{0:<15s}{1:>15s}{2:>10s}'.format('Component','Amplitude','Phase'))
# print('{0:<20s}{1:>9.7f}{2:>10.2f}°'.format('Δβ',np.abs(lat_eigvecs[0,4]),np.arctan2(lat_eigvecs[0,4].imag,lat_eigvecs[0,4].real)))
# print('{0:<20s}{1:>9.7f}{2:>10.2f}°'.format('Δpbar',np.abs(lat_eigvecs[1,4]),np.arctan2(lat_eigvecs[1,4].imag,lat_eigvecs[1,4].real)))
# print('{0:<20s}{1:>9.7f}{2:>10.2f}°'.format('Δrbar',np.abs(lat_eigvecs[2,4]),np.arctan2(lat_eigvecs[2,4].imag,lat_eigvecs[2,4].real)))
# print('{0:<20s}{1:>9.7f}{2:>10.2f}°'.format('Δξy',np.abs(lat_eigvecs[3,4]),np.arctan2(lat_eigvecs[3,4].imag,lat_eigvecs[3,4].real)))
# print('{0:<20s}{1:>9.7f}{2:>10.2f}°'.format('Δφ',np.abs(lat_eigvecs[4,4]),np.arctan2(lat_eigvecs[4,4].imag,lat_eigvecs[4,4].real)))
# print('{0:<20s}{1:>9.7f}{2:>10.2f}°'.format('Δψ',np.abs(lat_eigvecs[5,4]),np.arctan2(lat_eigvecs[5,4].imag,lat_eigvecs[5,4].real)))
# print(''.ljust(56,'='))
# print('The Dutch Roll Mode Damping Rate is {0:9.7f}'.format(sp_sig))
# print('The Dutch Roll Mode 99% Damping Time is {0:9.7f}'.format(sp_99damp))
# print('The Dutch Roll Mode Damped Frequency is {0:9.7f}'.format(sp_dampnatfreq))
# print('The Dutch Roll Mode Period is {0:9.7f}'.format(sp_period))
# print()
# print()

# All of the above code is replaced by these values.
printEigen(lat_eigvals[0],lat_eigvecs[:,0],'Roll Mode',True,False,2*Vo/b)
printEigen(lat_eigvals[1],lat_eigvecs[:,1],'Spiral Mode',True,False,2*Vo/b)
printEigen(lat_eigvals[4],lat_eigvecs[:,4],'Dutch Roll Mode',True,False,2*Vo/b)





print('APPROXIMATIONS:')
print('Short Period Approximation:')
Asp = R_yy*(R_rhox + CL_ahat)
Bsp = R_yy*(CL_a + CD_o) - Cm_q*(R_rhox + CL_ahat) - Cm_ahat*(R_rhox - CL_q)
Csp = -Cm_q*(CL_a + CD_o) - Cm_a*(R_rhox - CL_q)
sigma_sp = Vo/c*Bsp/Asp
omega_sp = Vo/c*np.abs(np.lib.scimath.sqrt(Bsp**2-4*Asp*Csp)/Asp)
eigenval_sp = c/(2*Vo)*np.complex(-sigma_sp,omega_sp)
printEigen(eigenval_sp,0,'Short Period Mode',False,False,2*Vo/c)


print('Long Period Approximation:')
sigma_d = g/Vo*CD_o/CLo
sigma_q = g/Vo*np.abs((CLo - CD_a)*Cm_q/(R_rhox*Cm_a + (CD_o + CL_a)*Cm_q))
Rs = R_rhox*Cm_a/(R_rhox*Cm_a + (CD_o + CL_a)*Cm_q)
sigma_psi = -g/Vo*R_gx*Rs*(R_rhox*Cm_q - R_yy*(CD_o + CL_a))/(R_rhox*Cm_a + (CD_o + CL_a)*Cm_q)
sigma_p = sigma_d + sigma_q + sigma_psi
omega_p = np.sqrt(2*(g/Vo)**2*Rs - (sigma_d + sigma_q)**2)
eigenval_p = c/(2*Vo)*np.complex(-sigma_p,omega_p)
printEigen(eigenval_p,0,'Long Period Mode',False,False,2*Vo/c)


print('Roll Mode Approximation: ')
sigma_r = (rho*Sw*b**2*Vo*Cl_p)/(4*Ixx)*b/(2*Vo)
printEigen(np.complex(sigma_r,0),0,'Roll Mode',False,False,2*Vo/b)

print('Spiral Mode Approximation: ')
sigma_s = -(g/Vo)*(Cl_b*Cn_r - Cl_r*Cn_b)/(Cl_b*Cn_p - Cl_p*Cn_b)*b/(2*Vo)
printEigen(np.complex(sigma_s,0),0,'Spiral Mode',False,False,2*Vo/b)


print('Dutch Roll Mode Approximation: ')
R_Ds = (Cl_b*(R_gy*R_rhoy*R_zz - (R_rhoy - CY_r)*Cn_p) - CY_b*Cl_r*Cn_p)/(R_rhoy*R_zz*Cl_p)
sigma_dr = Vo/b*(CY_b/R_rhoy + Cn_r/R_zz - Cl_r*Cn_p/(Cl_p*R_zz) + R_gy*(Cl_r*Cn_b - Cl_b*Cn_r)/(Cl_p*(Cn_b + CY_b*Cn_r/R_rhoy) - R_xx*R_Ds/Cl_p))*b/(2*Vo)
omega_dr = np.sqrt((1-CY_r/R_rhoy)*(Cn_b/R_zz) + CY_b*Cn_r/(R_rhoy*R_zz) + R_Ds - (1/4)*(CY_b/R_rhoy + Cn_r/R_zz)**2)
printEigen(np.complex(sigma_dr,omega_dr),0,'Dutch Roll Mode', False,False,2*Vo/b)