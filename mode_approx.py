# Approximation Methods for Dynamics Modes

import numpy as np
import cmath as cm

# SHORT PERIOD MODE Approximation
def spapprox(Vo,b,c,CL_a,CD_o,Cm_a,CL_q,Rrhox,Cm_q,CL_ahat,Cm_ahat,Ryy,printer):
    A = Ryy * (Rrhox + CL_ahat)
    B = (Ryy* (CL_a + CD_o)) - (Cm_q* (Rrhox + CL_ahat)) - (Cm_ahat* (Rrhox - CL_q))
    C = (-Cm_q* (CL_a + CD_o)) - (Cm_a* (Rrhox - CL_q))

    eig1 = ( (-B) + cm.sqrt((B**2)-(4*A*C)) ) / (2*A)
    eig2 = ( (-B) - cm.sqrt((B**2)-(4*A*C)) ) / (2*A)

    sigma = (-(2*Vo)/c) * eig1.real
    #sigma2 = (-(2*Vo)/c) * eig2.real
    w_d = ((2*Vo)/c) * abs(eig1.imag)
    #w_d2 = ((2*Vo)/c) * abs(eig2.imag)

    # s1 = (Vo/c) * (B/A)
    # s2 = (Vo/c) * (B/A)
    # w1 = (Vo/c) * abs( (cm.sqrt((B**2)-(4*A*C)))/ A)
    # w2 = (Vo/c) * abs( (cm.sqrt((B**2)-(4*A*C)))/ A) 

    nn = (-np.log(0.01)) / sigma
    period = ((2*np.pi) / w_d) 

    if printer:
        print("\n\n\n\nSHORT PERIOD MODE APPROXIMATION-----------\n")
        print("A: \t\t\t",A)
        print("B: \t\t\t",B)
        print("C: \t\t\t",C)
        print("\neig1: \t\t\t",eig1)
        print("eig2: \t\t\t",eig2)
        print("\nDamping Rate: \t\t",sigma)
        #print("Damping rate: \t\t",sigma2)
        print("Damped Nat Frequency: \t",w_d)
        #print("Damped Nat Frequency: \t",w_d2)
        print("99% Damping Time: \t",nn)
        print("Period: \t\t",period)

        # print(s1)
        # print(s2)
        # print(w1)
        # print(w2)
    
    return eig1,eig2,sigma,w_d

# PHUGOID MODE Approximation
def phapprox(g,Vo,b,c,CD_o,CLo,CD_a,Cm_q,Rrhox,Cm_a,CL_a,Rgx,Ryy,printer):
    Rps = (Rrhox * Cm_a) / ( (Rrhox * Cm_a) + ((CD_o + CL_a) * Cm_q) )
    sigma_D = (g/Vo) * (CD_o/CLo)
    sigma_q = (g/Vo) * ( ((CLo - CD_a) * Cm_q) / ( (Rrhox * Cm_a) + ((CD_o + CL_a) * Cm_q) ) )
    sigma_phi = (-g/Vo) * Rgx * Rps * \
        ( ((Rrhox * Cm_q) - (Ryy * (CD_o + CL_a)) ) / ( (Rrhox * Cm_a) + ((CD_o + CL_a) * Cm_q) ))
    
    sigma_p = sigma_D + sigma_q + sigma_phi
    w_d_p = np.sqrt((2 * ((g/Vo)**2) * Rps) - ((sigma_D + sigma_q) **2) )

    eig1 = (c/(2*Vo)) * (-sigma_p + w_d_p)
    eig2 = (c/(2*Vo)) * (-sigma_p - w_d_p)

    nn = (-np.log(0.01)) / sigma_p
    period = ((2*np.pi) / w_d_p) 

    if printer:
        print("\n\n\n\nPHUGOID PERIOD APPROXIMATION---------------\n")
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

