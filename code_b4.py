# Baseline code for developping the Dynamic Analysis code
# for the Aero 2.0 Glider Project

# B4 - B7
# Both Longitudinal and Lateral Matrices and eigenvalue calculations

import numpy as np 
import json
import eigensovler as eig



#### ------------------------------------------------------ ####
#### --------------------Unpack the JSON------------------- ####
input_file = 'BaselineGlider.json'

with open(input_file, "r") as f:
    aircraft = json.load(f)
Sw = aircraft["aircraft"]["wing_area[ft^2]"]
b = aircraft["aircraft"]["wing_span[ft]"]
W = aircraft["operating"]["weight[lbf]"]
theta = (aircraft["operating"]["climb[deg]"]) * (np.pi/ 180)
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

g = 32.17     #ft/s2

# Non-existent Thrust values for this case
T_V = 0
alpha_To = 0
z_T = 0
x_T = 0

Cm_o = 0



#### ------------------------------------------------------ ####
#### -----------------Preliminary Vals Calcs--------------- ####

c = Sw / b
Vo = np.sqrt((2*W*np.cos(theta)) / (rho*Sw*CLo))
CD_o = CD0 + (CD1 * CLo) + (CD2 * (CLo**2))
CD_a = (CD1 * CL_a) + (2 * CD2 * CLo * CL_a)
CT_V = T_V / (0.5 * rho * Vo * Sw)
z_To = (z_T * np.cos(alpha_To)) + (x_T * np.sin(alpha_To))

Rgx = (g * c) / (2 * (Vo**2))
Rgy = (g * b) / (2 * (Vo**2))
Rrhox = ((4 * W)/g) / (rho * Sw * c)
Rrhoy = ((4 * W)/g) / (rho * Sw * b)
Rxx = (8 * Ixx) / (rho * Sw * (b**3))
Ryy = (8 * Iyy) / (rho * Sw * (c**3))
Rzz = (8 * Izz) / (rho * Sw * (b**3))
Rxz = (8 * Ixz) / (rho * Sw * (b**3))


#### ------------------------------------------------------ ####
#### --------------------Matrix Filling-------------------- ####

# A and B for Longitudinal
A = np.zeros((6,6))
B = np.zeros((6,6))

A[0,0] = (-2 * CD_o) + (CT_V * np.cos(alpha_To))
A[1,0] = (-2 * CLo) - (CT_V * np.sin(alpha_To))
A[2,0] = (2 * Cm_o) + (CT_V * (z_To/c))
A[3,0] = np.cos(theta)
A[4,0] = -np.sin(theta)

A[0,1] = (CLo - CD_a)
A[1,1] = (-CL_a - CD_o)
A[2,1] = Cm_a
A[3,1] = np.sin(theta)
A[4,1] = np.cos(theta)

A[0,2] = -CD_q
A[1,2] = (-CL_q + Rrhox)
A[2,2] = Cm_q
A[5,2] = 1

A[0,5] = (-Rrhox * Rgx * np.cos(theta))
A[1,5] = (-Rrhox * Rgx * np.sin(theta))
A[3,5] = -np.sin(theta)
A[4,5] = -np.cos(theta)

B[0,0] = (Rrhox + CD_uhat)
B[1,0] = CL_uhat
B[2,0] = -Cm_uhat

B[0,1] = CD_ahat
B[1,1] = (Rrhox + CL_ahat)
B[2,1] = -Cm_ahat

B[2,2] = Ryy
B[3,3] = 1
B[4,4] = 1
B[5,5] = 1

# D and E for Lateral
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

D[0,2] = (CY_r - Rrhoy)
D[1,2] = Cl_r
D[2,2] = Cn_r
D[4,2] = np.tan(theta)
D[5,2] = 1/ (np.cos(theta))

D[0,4] = (Rrhoy * Rgy * np.cos(theta))

D[3,5] = np.cos(theta)

E[0,0] = Rrhoy

E[1,1] = Rxx
E[2,1] = -Rxz

E[1,2] = -Rxz
E[2,2] = Rzz

E[3,3] = 1
E[4,4] = 1
E[5,5] = 1



print(A)
print(B)
print(D)
print(E)


#### ------------------------------------------------------ ####
#### -------------------Eigenprob Calcs-------------------- ####

eigvals, eigvecs, vals = eig.eig_solve(A,B, char = True)
