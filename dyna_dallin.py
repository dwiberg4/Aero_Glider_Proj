# Baseline code for developping the Dynamic Analysis code
# for the Aero 2.0 Glider Project

# B4 - B7
# Both Longitudinal and Lateral Matrices and eigenvalue calculations

import numpy as np 
import json
import eigensovler as myeig
import eigen_approx as approx
import mode_approx as mapprox


#### ------------------------------------------------------ ####
#### --------------------Unpack the JSON------------------- ####
# input_file = 'BaselineGlider.json'
# input_file = 'edited_0000.json'
input_file = '2001.json'
# input_file = 'book_example.json'
if True:
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

#### ------------------------------------------------------ ####
#### -----------------Preliminary Vals Calcs--------------- ####
if True:
    g = 32.17     #ft/s2

    # Non-existent Thrust values for this case
    T_V = 0
    alpha_To = 0
    z_T = 0
    x_T = 0

    Cm_o = 0
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
if True:

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

    # G and H for Longitudinal APPROXIMATION
    G = np.zeros((2,2))
    H = np.zeros((2,2))

    G[0,0] = (-CL_a - CD_o)
    G[1,0] = Cm_a
    G[0,1] = (-CL_q + Rrhox)
    G[1,1] = Cm_q

    H[0,0] = (Rrhox + CL_ahat)
    H[1,0] = -Cm_ahat
    H[0,1] = 0
    H[1,1] = Ryy

    # J and K for Lateral APPROXIMATION
    J = np.zeros((2,2))
    K = np.zeros((2,2))


    # print(A)
    # print(B)
    # print(D)
    # print(E)


#### ------------------------------------------------------ ####
#### -------------------Eigenprob Calcs-------------------- ####
# Longitudinal Analysis
Long_eigs, Long_vecs, Long_vals, Long_dim_eigs = myeig.eig_solve(Vo,b,c,A,B, char = True, file = True, title = "Longitudinal")

# #_________0______1_______2_____3_______4____5_____6_____7_______8______9_____10__
# vals = [w_n_d, sigmas, taus, halves, nines,dubs, w_n, zetas, periods, amps, phas]


# for i in range(eigvals.shape[0]):
#     print("------------------------------------------------------------")
#     print("\n\nEigenvalue number: %d \t %0.12f %0.12fj" % ((i+1),eigvals[i].real,eigvals[i].imag))
#     print("\n\tThe component Amplitudes of the associated Eigenvector:")
#     for j in range(eigvecs.shape[1]):
#         print("\t", vals[9][i][j])
#     print("\n\tThe component Phases of the associated Eigenvector:")
#     for k in range(eigvecs.shape[1]):
#         print("\t", vals[10][i][j])
#     print("\n\tThe Damping Rate: \n\t", vals[1][i])

#     if eigvals[i].real > 0.0:
#         print("\n++The Mode for this Eigenvalue is Divergent\n")
#         print("\n\tThe Doubling Time is: \n\t",vals[5][i])
#         if eigvals[i].imag != 0.0:
#             print("\n++The Mode for this Eigenvalue is Oscillatory")
#             print("\n\tThe Damped Natural Frequency: \n\t", vals[0][i])
#             print("\n\tThe Period: \n\t", vals[8][i])
#         else:
#             print("\n++The mode is Non-Oscillatory")

#     elif eigvals[i].real < 0.0:
#         print("\n++The Mode for this Eigenvalue is Convergent\n")
#         print("\n\tThe 99% Damping Time is: \n\t",vals[4][i])
#         if eigvals[i].imag != 0.0:
#             print("\n++The Mode for this Eigenvalue is Oscillatory")
#             print("\n\tThe Damped Natural Frequency: \n\t", vals[0][i])
#             print("\n\tThe Period: \n\t", vals[8][i])
#         else:
#             print("\n++The mode is Non-Oscillatory")

#     else:
#         print("\n++The mode is Rigid Body Displacement Mode")


# Lateral Analysis       
Lat_eigs, Lat_vecs, Lat_vals, Lat_dim_eigs = myeig.eig_solve(Vo,b,c,D,E, char = True, file = True, title = "Lateral")


#### ------------------------------------------------------ ####
#### ------------------Mode APPROXIMATIONS----------------- ####

# One Eig Approx for Short Period... using the eig solver.
eigvals, eigvecs, vals = approx.eig_solve(Vo,b,c,G,H, char = True, file = True, title = "Longitudinal")


# SHORT PERIOD Approx
sp_eig1,sp_eig2,sigma_sp,w_d_sp = mapprox.spapprox(Vo,b,c,CL_a,CD_o,Cm_a,CL_q,Rrhox,Cm_q,CL_ahat,Cm_ahat,Ryy,True)
# PHUGOID Approx
ph_eig1,ph_eig2,sigma_ph,w_d_ph = mapprox.phapprox(g,Vo,b,c,CD_o,CLo,CD_a,Cm_q,Rrhox,Cm_a,CL_a,Rgx,Ryy,True)

# ROLL Approx
rl_eig,sigma_rl = mapprox.rlapprox(Vo,b,Sw,Ixx,rho,Cl_p,Rxx,True)
# SPIRAL Approx
sr_eig,sigma_sr = mapprox.srapprox(g,Vo,b,Cl_b,Cn_r,Cl_r,Cn_b,Cn_p,Cl_p,True)
# DUTCH ROLL Approx
dr_eig1,dr_eig2,sigma_dr,w_d_dr = mapprox.drapprox(Vo,b,CY_b,Cn_r,Cl_r,Cn_p,Cn_b,Cl_b,Cl_p,CY_r,Rgy,Rrhoy,Rzz,Rxx,True)



#### ------------------------------------------------------ ####
#### ------------------Handling Qualities------------------ ####
#### --------------------For Category B-------------------- ####
if True:
    # Manually input the indices corresponding to each MODE
    i_sp = (2,0)
    i_ph = (4,1)
    i_dr = (3,0)
    i_rl = 2        # 2 for BaselineGlider  5 for edited
    i_sr = 5        # 5 for BaselineGlider  2 for edited
    i_lp = (7,7)    #Not applicable for this case

    # SHORT PERIOD
    CW = W / (0.5*rho*(Vo**2)*Sw)
    acc_sens = CL_a / CW
    #print(acc_sens)
    w_n_SP = float(Long_vals[6][i_sp[1]])
    print(w_n_SP)

    CAP = (w_n_SP**2) / acc_sens
    print("\n\nThe CAP value for the Short Period is: \n\t\t\t",CAP)

    zeta_SP = float(Long_vals[7][i_sp[1]])
    #print("The zeta_SP value is: \t\t\t",zeta_SP)

    if ((CAP <= 0.038) or (zeta_SP <= 0.15)):
        hq_SP = 4
    if ((CAP > 0.038) and (zeta_SP > 0.15)):
        hq_SP = 3
    if ((0.038 < CAP < 10) and (0.2 < zeta_SP < 2.0)):
        hq_SP = 2
    if (0.085 < CAP < 3.6) and (0.3 < zeta_SP < 2.0):
        hq_SP = 1

    # PHUGOID
    zeta_PH = float(Long_vals[7][i_ph[1]])
    dubs_PH = float(Long_vals[5][i_ph[0]])

    if (zeta_PH > 0.04):
        hq_PH = 1
    elif (0 < zeta_PH <= 0.04):
        hq_PH = 2
    else:
        if (dubs_PH > 55):
            hq_PH = 3
        else:
            hq_PH = 4

    # DUTCH ROLL
    w_n_DR = float(Lat_vals[6][i_dr[1]])
    zeta_DR = float(Lat_vals[7][i_dr[1]])

    prod = w_n_DR * zeta_DR
    if (zeta_DR > 0.08 and prod > 0.15 and w_n_DR > 0.4):
        hq_DR = 1
    elif ((0.02 < zeta_DR <= 0.08) and (0.05 < zeta_DR <= 0.15) and (w_n_DR > 0.4) ):
        hq_DR = 2
    elif ((0.0 < zeta_DR <= 0.02) and (0.05 < zeta_DR <= 0.15) and (w_n_DR > 0.4) ):
        hq_DR = 3
    else:
        hq_DR = 4

    # ROLL
    tau_RL = float(Lat_vals[2][i_rl])

    if tau_RL < 1.4:
        hq_RL = 1
    elif 1.4 <= tau_RL < 3:
        hq_RL = 2
    elif 3 <= tau_RL < 10:
        hq_RL = 3
    else:
        hq_RL = 4

    # SPIRAL
    dubs_SR = float(Lat_vals[5][i_sr])

    if ((dubs_SR > 20) or (dubs_SR < 0)):
        hq_SR = 1
    elif 12 < dubs_SR <= 20:
        hq_SR = 2
    elif 4 < dubs_SR <= 12:
        hq_SR = 3
    else:
        hq_SR = 4

    # LATERAL PHUGOID
    if i_lp[0] < 6:
        w_n_LP = float(Lat_vals[6][i_lp[1]])
        zeta_LP = float(Lat_vals[7][i_lp[1]])

        sigma = w_n_LP * zeta_LP
        if (sigma > 0.5):
            hq_LP = 1
        elif (0.3 < zeta_LP <= 0.5):
            hq_LP = 2
        elif (0.15 < zeta_LP <= 0.3):
            hq_LP = 3
        else:
            hq_LP = 4
    else:
        hq_LP = "NA"


    print("\nHANDLING QUALITY SCORES:")
    print("Short Period:\t\t",hq_SP)
    print("Phugoid:\t\t",hq_PH)
    print("Roll:\t\t\t",hq_RL)
    print("Spiral:\t\t\t",hq_SR)
    print("Dutch Roll:\t\t",hq_DR)
    print("Lateral Phugoid:\t",hq_LP)

