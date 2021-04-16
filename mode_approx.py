# Approximation Methods for Dynamics Modes

import numpy as np
import cmath as cm

# SHORT PERIOD MODE Approximation
def spapprox(Vo,b,c,CL_a,CD_o,Cm_a,CL_q,Rrhox,Cm_q,CL_ahat,Cm_ahat,Ryy,printer):
    A = Ryy * (Rrhox + CL_ahat)
    B = (Ryy* (CL_a + CD_o)) - (Cm_q* (Rrhox + CL_ahat)) - (Cm_ahat* (Rrhox - CL_q))
    C = (-Cm_q* (CL_a + CD_o)) - (Cm_a* (Rrhox - CL_q))

    # tmpB = np.complex(-B,0)
    # print(tmpB)
    eig1 = ( (-B) + cm.sqrt((B**2)-(4*A*C)) ) / (2*A)
    eig2 = ( (-B) - cm.sqrt((B**2)-(4*A*C)) ) / (2*A)
    # print(A)
    # print(B)
    # print(C)
    # thise = cm.sqrt((B**2)-(4*A*C))
    # print("complex part: ",thise)
    # print(tmpB+thise)
    # print(tmpB-thise)
    # print(2*A)
    # print((tmpB+thise)/(2*A))
    # print((tmpB-thise)/(2*A))
    # print(np.complex(1.23423,23.1234)


    sigma = (-(2*Vo)/c) * eig1.real
    w_d = ((2*Vo)/c) * abs(eig1.imag)

    nn = (-np.log(0.01)) / sigma
    period = ((2*np.pi) / w_d) 

    s2 = (Vo/c) * (B/A)
    w2 = (Vo/c) * abs( (cm.sqrt((B**2)-(4*A*C)))/ A)
    eig12 = (c/(2*Vo)) * np.complex(-s2, w2)
    eig22 = (c/(2*Vo)) * np.complex(-s2, -w2)
    nn2 = (-np.log(0.01)) /s2
    period2 = ((2*np.pi) / w2)

    if printer:
        print("\n\n\n\nSHORT PERIOD MODE APPROXIMATION-------------")
        print("--------------------------------------------")
        print("A: \t\t\t",A)
        print("B: \t\t\t",B)
        print("C: \t\t\t",C)
        print("\neig1: \t\t\t",eig12)
        print("eig2: \t\t\t",eig22)
        print("\nDamping Rate: \t\t",s2)
        print("Damped Nat Frequency: \t",w2)
        print("99% Damping Time: \t",nn2)
        print("Period: \t\t",period2)

        print("\n----------Second Method:")
        print("\tCOMPLEX MATH ISSUES ******** SOMETIMES DOESN'T WORK")
        print("\teig12: \t\t\t",eig1)
        print("\teig22: \t\t\t",eig2)
        print("\n\tDamping Rate: \t\t",sigma)
        print("\tDamped Nat Frequency: \t",w_d)
        print("\t99% Damping Time: \t",nn)
        print("\tPeriod: \t\t",period)
    
    return eig12,eig22,s2,w2

# PHUGOID MODE Approximation
def phapprox(g,Vo,b,c,CD_o,CLo,CD_a,Cm_q,Rrhox,Cm_a,CL_a,Rgx,Ryy,printer):
    Rps = (Rrhox * Cm_a) / ( (Rrhox * Cm_a) + ((CD_o + CL_a) * Cm_q) )
    sigma_D = (g/Vo) * (CD_o/CLo)
    sigma_q = (g/Vo) * ( ((CLo - CD_a) * Cm_q) / ( (Rrhox * Cm_a) + ((CD_o + CL_a) * Cm_q) ) )
    sigma_phi = (-g/Vo) * Rgx * Rps * \
        ( ((Rrhox * Cm_q) - (Ryy * (CD_o + CL_a)) ) / ( (Rrhox * Cm_a) + ((CD_o + CL_a) * Cm_q) ))
    
    sigma_p = sigma_D + sigma_q + sigma_phi
    w_d_p = np.sqrt((2 * ((g/Vo)**2) * Rps) - ((sigma_D + sigma_q) **2) )

    eig1 = (c/(2*Vo)) * np.complex(-sigma_p, w_d_p)
    eig2 = (c/(2*Vo)) * np.complex(-sigma_p, -w_d_p)
    # print("\n\n\nsigma_p is :",sigma_p)
    # print("w_d_p is: ",w_d_p)
    # print("dim time is: ",(c/(2*Vo)))
    # print("complex part is:", (-sigma_p + w_d_p))
    # print("eig1 is:", eig1)
    # print("eig2 is:", eig2)

    nn = (-np.log(0.01)) / sigma_p
    period = ((2*np.pi) / w_d_p) 

    if printer:
        print("\n\n\n\nPHUGOID MODE APPROXIMATION------------------")
        print("--------------------------------------------")
        print("Rps: \t\t\t",Rps)
        print("sigma_D: \t\t",sigma_D)
        print("sigma_q: \t\t",sigma_q)
        print("sigma_phi: \t\t",sigma_phi)
        print("\neig1: \t\t\t",eig1)
        print("eig2: \t\t\t",eig2)
        print("\nDamping Rate: \t\t",sigma_p)
        print("Damped Nat Frequency \t",w_d_p)
        print("99% Damping Time: \t",nn)
        print("Period: \t\t",period)
    
    return eig1,eig2,sigma_p,w_d_p
        
# ROLL MODE Approximation
def rlapprox(Vo,b,Sw,Ixx,rho,Cl_p,Rxx,printer):
    eig = np.complex((Cl_p / Rxx),0)

    sigma = ((-rho * Sw * (b**2) * Vo)/ (4 * Ixx)) * Cl_p

    nn = (-np.log(0.01)) / sigma

    if printer:
        print("\n\n\n\nROLL MODE APPROXIMATION---------------------")
        print("--------------------------------------------")
        print("\neig: \t\t\t",eig)
        print("\nDamping Rate: \t\t",sigma)
        print("99% Damping Time: \t",nn)
    
    return eig,sigma

# SPIRAL MODE Approximation
def srapprox(g,Vo,b,Cl_b,Cn_r,Cl_r,Cn_b,Cn_p,Cl_p,printer):
    coeffs = ( ((Cl_b*Cn_r) - (Cl_r*Cn_b)) / ((Cl_b*Cn_p) - (Cl_p*Cn_b)) )

    eig = np.complex( (((-g * b) / (2* (Vo**2))) * coeffs), 0)

    sigma = (g/Vo) * coeffs

    nn = (-np.log(0.01)) / sigma
    dub = (-np.log(2)) / sigma

    if printer:
        print("\n\n\n\nSPIRAL MODE APPROXIMATION-------------------")
        print("--------------------------------------------")
        print("coeffs: \t\t",coeffs)
        print("\neig: \t\t\t",eig)
        print("\nDamping Rate: \t\t",sigma)
        print("99% Damping Time: \t",nn)
        print("Doubling Time: \t\t",dub)
    
    return eig,sigma

# DUTCH ROLL MODE Approximation
def drapprox(Vo,b,CY_b,Cn_r,Cl_r,Cn_p,Cn_b,Cl_b,Cl_p,CY_r,Rgy,Rrhoy,Rzz,Rxx,printer):
    top = (Cl_b * ((Rgy*Rrhoy*Rzz) - ((Rrhoy - CY_r)* Cn_p))) - (CY_b*Cl_r*Cn_p)
    RDs = top / (Rrhoy * Rzz * Cl_p)

    one = (1- (CY_r/Rrhoy)) * (Cn_b / Rzz)
    two = ((CY_b*Cn_r) / (Rrhoy*Rzz))
    tre = ((CY_b/Rrhoy) + (Cn_r/Rzz))
    tog = one + two + RDs - (0.25 * (tre**2))
    w_d = ((2*Vo) / b) * np.sqrt(tog)

    one = (CY_b / Rrhoy) + (Cn_r / Rzz)
    two = ((Cl_r*Cn_p) / (Cl_p*Rzz))
    ctop = Rgy * ((Cl_r*Cn_b) - (Cl_b*Cn_r))
    cbot = Cl_p * (Cn_b + (CY_b* (Cn_r/Rrhoy)) )
    four = Rxx * (RDs/Cl_p)
    tog = one - two + (ctop/cbot) - four
    sigma = (-Vo/b) * tog


    eig1 = (b/(2*Vo)) * np.complex(-sigma, w_d)
    eig2 = (b/(2*Vo)) * np.complex(-sigma, -w_d)

    nn = (-np.log(0.01)) / sigma
    period = ((2*np.pi) / w_d) 

    if printer:
        print("\n\n\n\nDUTCH ROLL MODE APPROXIMATION-------------")
        print("--------------------------------------------")
        print("RDs: \t\t\t",RDs)
        print("\neig1: \t\t\t",eig1)
        print("eig2: \t\t\t",eig2)
        print("\nDamping Rate: \t\t",sigma)
        print("Damped Nat Frequency: \t",w_d)
        print("99% Damping Time: \t",nn)
        print("Period: \t\t",period)
    
    return eig1,eig2,sigma,w_d


