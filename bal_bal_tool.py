# Ballast Moment Balance Tool

# A function to take:
# 1. The non-ballast weight of the aircraft
# 2. The Total Cm produced by the wing(s)
# 3. The Total CL produced by the wing(s)
# 4. The Location of the CG
# 5. The x-location of the VSTAB

# And output one, given the other of:
# A. Given the available moment arm, output the ballast mass
# B. Given the ballast mass, output the necessary moment arm


import numpy as np
import glider_funcs as gfs
import matplotlib.pyplot as plt

def balbal(params,results,plot, **kwargs):
    Cm = results["Cm"]
    CL = results["CL"]
    rho = params["rho"]
    V = params["vel"]
    Sw = params["Sw"]

    # Unpack vars
    l_arm = -(params["CG"][0])                # Moment are for the lift wrt the CG


    # the origin is the 1/4 chord of the main wing
    # The CG can be fore or aft (usually fore) of this origin point
    # Moment Balance: (wrt CG)
    # SUM_M = 0 = +(CL*-l_arm)  -Cm  +(m_b*b_arm)

    #masses = np.linspace(0.25,1.75,13)
    #F_b = masses / (0.5*rho*(V**2)*Sw)
    arms = np.linspace(params["max_arm"],params["x_vstab"],40)

    C_F_b_2 = np.zeros((arms.size))
    masses2 = np.zeros((arms.size))
    for i in range(arms.size):
        C_F_b_2[i] = mmt_bal(Cm,CL,l_arm,rho,V,Sw, b_arm = arms[i], coeff = True)
        masses2 [i] = C_F_b_2[i] * (0.5*rho*(V**2)*Sw)
    
    if plot == True:
        title = "Ballast Weight as function of trimmed location"
        gfs.gen_plotter(arms,masses2,xlabel="The x-location on the dowel",ylabel="Ballast Amount [Mass]",title=title,save="ballast_balance.png")

    if "b_arm" in kwargs:
        b_arm = kwargs["b_arm"]
        traget_coeff = mmt_bal(Cm,CL,l_arm,rho,V,Sw, b_arm = b_arm, coeff = True)
        target_mass = mmt_bal(Cm,CL,l_arm,rho,V,Sw, b_arm = b_arm, mass = True)
        print("\nGiven this desired moment arm, the necessary ballast mass is: ", target_mass)
        return target_mass
    if "C_F_b" in kwargs:
        C_F_b = kwargs["C_F_b"]
        target_arm = mmt_bal(Cm,CL,l_arm,rho,V,Sw, C_F_b = C_F_b)
        return target_arm
    if "b_mass" in kwargs:
        b_mass = kwargs["b_mass"]
        target_arm = mmt_bal(Cm,CL,l_arm,rho,V,Sw, b_mass = b_mass)
        return target_arm


def mmt_bal(Cm,CL,l_arm,rho,V,Sw, **kwargs):
    if "b_arm" in kwargs:
        b_arm = kwargs["b_arm"]
        C_F_b = (Cm + (CL*l_arm)) / b_arm
        if "coeff" in kwargs and kwargs["coeff"] == True:
            return C_F_b
        if "mass" in kwargs and kwargs["mass"] == True:
            return C_F_b * (0.5*rho*(V**2)*Sw)
    if "C_F_b" in kwargs:
        C_F_b = kwargs["C_F_b"]
        b_arm = (Cm + (CL*l_arm)) / C_F_b
        return b_arm
    if "b_mass" in kwargs:
        b_mass = kwargs["b_mass"]
        b_arm = (Cm + (CL*l_arm)) / (b_mass/(0.5*rho*(V**2)*Sw))
        return b_arm




# def main():
#     W_a = 1.44
#     Cm = 0.0025
#     CL = 0.5
#     CG = [-0.0835,0,0]
#     x_vstab = -2.643
#     rho = 2.048e-3
#     V = 21.03
#     Sw = 5.75
#     balbal(W_a,Cm,CL,CG,x_vstab,rho,V,Sw)
#     min_mass = mmt_bal(Cm,CL,-CG[0],rho,V,Sw, b_arm = (6+x_vstab), mass = True)
#     print("The minimum possible mass is: ", min_mass)
    



# main()