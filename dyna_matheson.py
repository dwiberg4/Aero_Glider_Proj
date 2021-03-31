import numpy as np
from scipy.linalg import eig
import json
input = 'BaselineGlider.json'

with open(input,'r') as f:
    aircraft = json.load(f)
f.close()

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
R_xx = 8*Ixx / (rho*Sw*b**3)                          # Equation 10.73 Ch. 7 Overview
R_yy = 8*Iyy / (rho*Sw*c**3)                          # Equation 10.73 Ch. 7 Overview
R_zz = 8*Izz / (rho*Sw*b**3)                          # Equation 10.73 Ch. 7 Overview
R_xz = 8*Ixz / (rho*Sw*b**3)                          # Equation 10.73 Ch. 7 Overview

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
eigvals, eigvecs = eig(C)

# Calculate the magnitude of each Eigenvalue
mag = np.abs(eigvals)

# Sort by eigvals, eigvecs by largest magnitude
idx = mag.argsort()[::-1]
eigvals = eigvals[idx]
eigvecs = eigvecs[:,idx]
print('Short Period Mode is given by {0:9.7f} ± {1:9.7f} i.'.format(eigvals[0].real,eigvals[0].imag))
print(''.ljust(56,'='))
sp_sig = -(eigvals[0].real)*2*Vo/c
sp_99damp = np.log(0.01)/-sp_sig
sp_dampnatfreq = np.abs(eigvals[0].imag)*2*Vo/c
sp_dampratio = -(eigvals[0] +eigvals[1])/(2*np.sqrt(eigvals[0]*eigvals[1]))
sp_undampnatfreq = 2*Vo/c*np.sqrt(eigvals[0]*eigvals[1])
sp_period = 2*np.pi/sp_dampnatfreq
print('{0:<15s}{1:>15s}{2:>10s}'.format('Component','Amplitude','Phase'))
print('{0:<20s}{1:>9.7f}{2:>10.2f}°'.format('Δμ',np.abs(eigvecs[0,0]),np.arctan2(eigvecs[0,0].imag,eigvecs[0,0].real)))
print('{0:<20s}{1:>9.7f}{2:>10.2f}°'.format('Δα',np.abs(eigvecs[1,0]),np.arctan2(eigvecs[1,0].imag,eigvecs[1,0].real)))
print('{0:<20s}{1:>9.7f}{2:>10.2f}°'.format('Δqbar',np.abs(eigvecs[2,0]),np.arctan2(eigvecs[2,0].imag,eigvecs[2,0].real)))
print('{0:<20s}{1:>9.7f}{2:>10.2f}°'.format('Δξx',np.abs(eigvecs[3,0]),np.arctan2(eigvecs[3,0].imag,eigvecs[3,0].real)))
print('{0:<20s}{1:>9.7f}{2:>10.2f}°'.format('Δξz',np.abs(eigvecs[4,0]),np.arctan2(eigvecs[4,0].imag,eigvecs[4,0].real)))
print('{0:<20s}{1:>9.7f}{2:>10.2f}°'.format('Δθ',np.abs(eigvecs[5,0]),np.arctan2(eigvecs[5,0].imag,eigvecs[5,0].real)))
print(''.ljust(56,'='))
print('The Short Period Damping Rate is {0:9.7f}'.format(sp_sig))
print('The Short Period 99% Damping Time is {0:9.7f}'.format(sp_99damp))
print('The Short Period Damped Frequency is {0:9.7f}'.format(sp_dampnatfreq))
print('The Short Period Period is {0:9.7f}'.format(sp_period))
print()
print()
   
print('Long Period Mode is given by {0:9.7f} ± {1:9.7f} i.'.format(eigvals[2].real,eigvals[2].imag))
print(''.ljust(56,'='))
sp_sig = -(eigvals[2].real)*2*Vo/c
sp_99damp = np.log(0.01)/-sp_sig
sp_dampnatfreq = np.abs(eigvals[2].imag)*2*Vo/c
sp_dampratio = -(eigvals[2] +eigvals[3])/(2*np.sqrt(eigvals[2]*eigvals[3]))
sp_undampnatfreq = 2*Vo/c*np.sqrt(eigvals[2]*eigvals[3])
sp_period = 2*np.pi/sp_dampnatfreq
print('{0:<15s}{1:>15s}{2:>10s}'.format('Component','Amplitude','Phase'))
print('{0:<20s}{1:>9.7f}{2:>10.2f}°'.format('Δμ',np.abs(eigvecs[0,2]),np.arctan2(eigvecs[0,2].imag,eigvecs[0,2].real)))
print('{0:<20s}{1:>9.7f}{2:>10.2f}°'.format('Δα',np.abs(eigvecs[1,2]),np.arctan2(eigvecs[1,2].imag,eigvecs[1,2].real)))
print('{0:<20s}{1:>9.7f}{2:>10.2f}°'.format('Δqbar',np.abs(eigvecs[2,2]),np.arctan2(eigvecs[2,2].imag,eigvecs[2,2].real)))
print('{0:<20s}{1:>9.7f}{2:>10.2f}°'.format('Δξx',np.abs(eigvecs[3,2]),np.arctan2(eigvecs[3,2].imag,eigvecs[3,2].real)))
print('{0:<20s}{1:>9.7f}{2:>10.2f}°'.format('Δξz',np.abs(eigvecs[4,2]),np.arctan2(eigvecs[4,2].imag,eigvecs[4,2].real)))
print('{0:<20s}{1:>9.7f}{2:>10.2f}°'.format('Δθ',np.abs(eigvecs[5,2]),np.arctan2(eigvecs[5,2].imag,eigvecs[5,2].real)))
print(''.ljust(56,'='))
print('The Short Period Damping Rate is {0:9.7f}'.format(sp_sig))
print('The Short Period 99% Damping Time is {0:9.7f}'.format(sp_99damp))
print('The Short Period Damped Frequency is {0:9.7f}'.format(sp_dampnatfreq))
print('The Short Period Period is {0:9.7f}'.format(sp_period))




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




