import numpy as np
from scipy.linalg import eig
import json
import matplotlib.pyplot as plt



input = 'edited_0000.json'

def printEigen(eigval,eigvec,mode,lat,long,dimensional,filename = None):
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
    if filename != None:
        filename.write(mode.upper())
        filename.write('\n')
        filename.write('The dimensionless eigenvalue is given by {0:9.7f} ± {1:9.7f} i.\n'.format(eigval.real,np.abs(eigval.imag)))
        
    if eigvec is None:
        0
    else:
        if filename != None:
            filename.write('The dimensionless eigenvector is given by: \n {0:9.7f} \n {1:9.7f} \n {2:9.7f} \n {3:9.7f} \n {4:9.7f} \n {5:9.7f} \n'.format(eigvec[0],\
            eigvec[1],eigvec[2],eigvec[3],eigvec[4],eigvec[5]))
    sig = -eigval.real
    omega = eigval.imag
    damped = (0>eigval.real)
    damp_nat_freq = np.abs(eigval.imag)*dimensional
    damp_rate = -eigval.real*dimensional
    damp_ratio = sig*dimensional
    undamp_freq = 0.0
    if damp_nat_freq > 0.0001:
        # If there is an imaginary part calculate and print all of these values.
        eig1 = np.complex(-eigval.real,eigval.imag)*dimensional
        eig2 = np.complex(-eigval.real,-eigval.imag)*dimensional
        undamp_freq = (np.sqrt(eig1*eig2)).real
        period = 2*np.pi/(damp_nat_freq)
        damp_ratio = ((eig1 + eig2)/(2*np.sqrt(eig1*eig2))).real
        if long:
            if filename != None:
                filename.write('\n'.ljust(56,'='))
                filename.write('\n')
                filename.write('{0:<15s}{1:>15s}{2:>10s}\n'.format('Component','Amplitude','Phase'))
                filename.write('{0:<20s}{1:>9.7f}{2:>10.2f}°\n'.format('Δμ',np.abs(eigvec[0]),-np.degrees(np.arctan2(eigvec[0].imag,eigvec[0].real))))
                filename.write('{0:<20s}{1:>9.7f}{2:>10.2f}°\n'.format('Δα',np.abs(eigvec[1]),-np.degrees(np.arctan2(eigvec[1].imag,eigvec[1].real))))
                filename.write('{0:<20s}{1:>9.7f}{2:>10.2f}°\n'.format('Δqbar',np.abs(eigvec[2]),-np.degrees(np.arctan2(eigvec[2].imag,eigvec[2].real))))
                filename.write('{0:<20s}{1:>9.7f}{2:>10.2f}°\n'.format('Δξx',np.abs(eigvec[3]),-np.degrees(np.arctan2(eigvec[3].imag,eigvec[3].real))))
                filename.write('{0:<20s}{1:>9.7f}{2:>10.2f}°\n'.format('Δξz',np.abs(eigvec[4]),-np.degrees(np.arctan2(eigvec[4].imag,eigvec[4].real))))
                filename.write('{0:<20s}{1:>9.7f}{2:>10.2f}°\n'.format('Δθ',np.abs(eigvec[5]),-np.degrees(np.arctan2(eigvec[5].imag,eigvec[5].real))))
                filename.write('\n'.ljust(56,'='))
                filename.write('\n')
        if lat:
            if filename != None:
                filename.write('\n'.ljust(56,'='))
                filename.write('\n')
                filename.write('{0:<15s}{1:>15s}{2:>10s}\n'.format('Component','Amplitude','Phase'))
                filename.write('{0:<20s}{1:>9.7f}{2:>10.2f}°\n'.format('Δβ',np.abs(eigvec[0]),-np.degrees(np.arctan2(eigvec[0].imag,eigvec[0].real))))
                filename.write('{0:<20s}{1:>9.7f}{2:>10.2f}°\n'.format('Δpbar',np.abs(eigvec[1]),-np.degrees(np.arctan2(eigvec[1].imag,eigvec[1].real))))
                filename.write('{0:<20s}{1:>9.7f}{2:>10.2f}°\n'.format('Δrbar',np.abs(eigvec[2]),-np.degrees(np.arctan2(eigvec[2].imag,eigvec[2].real))))
                filename.write('{0:<20s}{1:>9.7f}{2:>10.2f}°\n'.format('Δξy',np.abs(eigvec[3]),-np.degrees(np.arctan2(eigvec[3].imag,eigvec[3].real))))
                filename.write('{0:<20s}{1:>9.7f}{2:>10.2f}°\n'.format('Δφ',np.abs(eigvec[4]),-np.degrees(np.arctan2(eigvec[4].imag,eigvec[4].real))))
                filename.write('{0:<20s}{1:>9.7f}{2:>10.2f}°\n'.format('Δψ',np.abs(eigvec[5]),-np.degrees(np.arctan2(eigvec[5].imag,eigvec[5].real))))
                filename.write('\n'.ljust(56,'='))
                filename.write('\n')
        if filename != None:
            filename.write('The {1:<20s} Period is {0:9.7f} s.\n'.format(period,mode))   
            filename.write('The {1:<20s} Damped Natural Frequency is {0:9.7f} rad/s\n'.format(damp_nat_freq,mode))    
            filename.write('The {1:<20s} Damping Ratio is {0:9.7f}\n'.format(damp_ratio,mode)) 
            filename.write('The {1:<20s} Undamped Natural Frequency is {0:9.7f} rad/s\n'.format(undamp_freq,mode))             
    else:
        # If there is not an imaginary part then do not filename.write phase in the table.
        if long:
            if filename != None:
                filename.write('\n'.ljust(56,'='))
                filename.write('\n')
                filename.write('{0:<15s}{1:>15s}\n'.format('Component','Amplitude'))
                filename.write('{0:<20s}{1:>9.7f}\n'.format('Δμ',np.abs(eigvec[0])))
                filename.write('{0:<20s}{1:>9.7f}\n'.format('Δα',np.abs(eigvec[1])))
                filename.write('{0:<20s}{1:>9.7f}\n'.format('Δqbar',np.abs(eigvec[2])))
                filename.write('{0:<20s}{1:>9.7f}\n'.format('Δξx',np.abs(eigvec[3])))
                filename.write('{0:<20s}{1:>9.7f}\n'.format('Δξz',np.abs(eigvec[4])))
                filename.write('{0:<20s}{1:>9.7f}\n'.format('Δθ',np.abs(eigvec[5])))
                filename.write('\n'.ljust(56,'='))
                filename.write('\n')
        if lat:
            if filename != None:
                filename.write('\n'.ljust(56,'=')) 
                filename.write('\n')       
                filename.write('{0:<15s}{1:>15s}\n'.format('Component','Amplitude'))
                filename.write('{0:<20s}{1:>9.7f}\n'.format('Δβ',np.abs(eigvec[0])))
                filename.write('{0:<20s}{1:>9.7f}\n'.format('Δpbar',np.abs(eigvec[1])))
                filename.write('{0:<20s}{1:>9.7f}\n'.format('Δrbar',np.abs(eigvec[2])))
                filename.write('{0:<20s}{1:>9.7f}\n'.format('Δξy',np.abs(eigvec[3])))
                filename.write('{0:<20s}{1:>9.7f}\n'.format('Δφ',np.abs(eigvec[4])))
                filename.write('{0:<20s}{1:>9.7f}\n'.format('Δψ',np.abs(eigvec[5])))
                filename.write('\n'.ljust(56,'='))
                filename.write('\n')
    if damped:
        # If the system is damped, filename.write a 99% Damping Time.
        damp_time = -np.log(0.01)/(sig*dimensional)          # 99% damping time
        if filename != None:
            filename.write('The {0:<20s} 99% Damping Time is {1:9.7f} s.\n'.format(mode,damp_time))
    else:
        # If the system is not damped, filename.write a Doubling Time instead.
        damp_time = -np.log(2)/(sig*dimensional)             # Doubling Time
        if filename != None:
            filename.write('The {0:<20s} Doubling Time is {1:9.7f} s.\n'.format(mode,damp_time))
    if filename != None:
        filename.write('The {1:<20s} Damping Rate is {0:9.7f} s^(-1)'.format(damp_rate,mode))
        filename.write('\n')
        filename.write('\n')
    return damp_ratio,damp_time,undamp_freq

    




with open(input,'r') as f:
    aircraft = json.load(f)
f.close()


with open('longitudinal.txt','w',encoding = 'utf-8') as longf, open('lateral.txt','w',encoding = 'utf-8') as latf:

## Read in all values from the JSON file.
    Sw = aircraft["aircraft"]["wing_area[ft^2]"]
    b = aircraft["aircraft"]["wing_span[ft]"]
    W = aircraft["operating"]["weight[lbf]"]
    # gamma = (aircraft["operating"]["climb[deg]"]) * (np.pi/ 180)
    rho = aircraft["operating"]["density[slugs/ft^3]"]
    #rho = 0.0023769
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
    #CD_o = 0.05
    CD_a = CD1*CL_a + 2*CD2*CLo*CL_a                    # Equation 10.57 Ch. 7 Overview
    #CD_a = 0.35
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
    longf.write('The longitudinal A Matrix is given as: \n')
    for line in A:
        longf.write(' {0:>20.12f} {1:>20.12f} {2:>20.12f} {3:>20.12f} {4:>20.12f} {5:>20.12f} \n'.format(*line))
    longf.write('\n')
    longf.write('The longitudinal B Matrix is given as: \n')
    for line in B:
        longf.write(' {0:>20.12f} {1:>20.12f} {2:>20.12f} {3:>20.12f} {4:>20.12f} {5:>20.12f} \n'.format(*line))
    longf.write('\n')

    # Calculate the magnitude of each Eigenvalue
    mag = np.abs(long_eigvals)

    # Get index that sorts the magnitudes from largest to smallest.
    idx = mag.argsort()[::-1]
    #sort both values and vectors by that index.
    long_eigvals = long_eigvals[idx]
    long_eigvecs = long_eigvecs[:,idx]

    sp_damp_ratio,unneed,sp_undamp_freq = printEigen(long_eigvals[0],long_eigvecs[:,0],'Short Period Mode',False,True,2*Vo/c,longf)
    ph_damp_ratio,ph_damp_time,unneed1 = printEigen(long_eigvals[2],long_eigvecs[:,2],'Long Period Mode',False,True,2*Vo/c,longf)

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
    latf.write('The lateral A Matrix is given as: \n')
    for line in D:
        latf.write(' {0:>20.12f} {1:>20.12f} {2:>20.12f} {3:>20.12f} {4:>20.12f} {5:>20.12f} \n'.format(*line))
    latf.write('\n')
    latf.write('The lateral B Matrix is given as: \n')
    for line in E:
        latf.write(' {0:>20.12f} {1:>20.12f} {2:>20.12f} {3:>20.12f} {4:>20.12f} {5:>20.12f} \n'.format(*line))
    latf.write('\n')

    # Calculate the magnitude of each Eigenvalue if there is an imaginary part give it a 1000 penalty
    mag = np.where(lat_eigvals.imag == 0,np.abs(lat_eigvals),np.abs(lat_eigvals)-1000)

    # Sort by lat_eigvals, lat_eigvecs by largest magnitude if there is no imaginary part.
    idx = mag.argsort()[::-1]
    lat_eigvals = lat_eigvals[idx]
    lat_eigvecs = lat_eigvecs[:,idx]


    roll_sig,unneed,unneed1 = printEigen(lat_eigvals[0],lat_eigvecs[:,0],'Roll Mode',True,False,2*Vo/b,latf)
    unneed,spiral_time,unneed1 = printEigen(lat_eigvals[1],lat_eigvecs[:,1],'Spiral Mode',True,False,2*Vo/b,latf)
    dutch_roll_damp_ratio,unneed,dutch_roll_ud_freq = printEigen(lat_eigvals[4],lat_eigvecs[:,4],'Dutch Roll Mode',True,False,2*Vo/b,latf)

    longf.write('APPROXIMATIONS:\n\n')

    Asp = R_yy*(R_rhox + CL_ahat)
    Bsp = R_yy*(CL_a + CD_o) - Cm_q*(R_rhox + CL_ahat) - Cm_ahat*(R_rhox - CL_q)
    Csp = -Cm_q*(CL_a + CD_o) - Cm_a*(R_rhox - CL_q)
    sigma_sp = Vo/c*Bsp/Asp
    omega_sp = Vo/c*np.abs(np.lib.scimath.sqrt(Bsp**2-4*Asp*Csp)/Asp)
    eigenval_sp = c/(2*Vo)*np.complex(-sigma_sp,omega_sp)
    printEigen(eigenval_sp,None,'Short Period Mode Approximation',False,False,2*Vo/c,longf)



    sigma_d = g/Vo*CD_o/CLo
    sigma_q = g/Vo*np.abs((CLo - CD_a)*Cm_q/(R_rhox*Cm_a + (CD_o + CL_a)*Cm_q))
    Rs = R_rhox*Cm_a/(R_rhox*Cm_a + (CD_o + CL_a)*Cm_q)
    sigma_psi = -g/Vo*R_gx*Rs*(R_rhox*Cm_q - R_yy*(CD_o + CL_a))/(R_rhox*Cm_a + (CD_o + CL_a)*Cm_q)
    sigma_p = sigma_d + sigma_q + sigma_psi
    omega_p = np.sqrt(2*(g/Vo)**2*Rs - (sigma_d + sigma_q)**2)
    eigenval_p = c/(2*Vo)*np.complex(-sigma_p,omega_p)
    printEigen(eigenval_p,None,'Phugoid Mode Approximation',False,False,2*Vo/c,longf)



    sigma_r = (rho*Sw*b**2*Vo*Cl_p)/(4*Ixx)*b/(2*Vo)
    printEigen(np.complex(sigma_r,0),None,'Roll Mode Approximation',False,False,2*Vo/b,latf)


    sigma_s = -(g/Vo)*(Cl_b*Cn_r - Cl_r*Cn_b)/(Cl_b*Cn_p - Cl_p*Cn_b)*b/(2*Vo)
    printEigen(np.complex(sigma_s,0),None,'Spiral Mode Approximation',False,False,2*Vo/b,latf)


    R_Ds = (Cl_b*(R_gy*R_rhoy*R_zz - (R_rhoy - CY_r)*Cn_p) - CY_b*Cl_r*Cn_p)/(R_rhoy*R_zz*Cl_p)
    sigma_dr = Vo/b*(CY_b/R_rhoy + Cn_r/R_zz - Cl_r*Cn_p/(Cl_p*R_zz) + R_gy*(Cl_r*Cn_b - Cl_b*Cn_r)/(Cl_p*(Cn_b + CY_b*Cn_r/R_rhoy)) - R_xx*R_Ds/Cl_p)*b/(2*Vo)
    omega_dr = np.sqrt((1-CY_r/R_rhoy)*(Cn_b/R_zz) + CY_b*Cn_r/(R_rhoy*R_zz) + R_Ds - (1/4)*(CY_b/R_rhoy + Cn_r/R_zz)**2)
    printEigen(np.complex(sigma_dr,omega_dr),None,'Dutch Roll Mode Approximation', False,False,2*Vo/b,latf)


    sp_eigval1 = long_eigvals[0]*2*Vo/c
    sp_eigval2 = long_eigvals[1]*2*Vo/c
    ph_eigval1 = long_eigvals[2]*2*Vo/c
    ph_eigval2 = long_eigvals[3]*2*Vo/c
    roll_eigval1 = lat_eigvals[0]*2*Vo/b
    spiral_eigval1 = lat_eigvals[1]*2*Vo/b
    dutch_eigval1 = lat_eigvals[4]*2*Vo/b
    dutch_eigval2 = lat_eigvals[5]*2*Vo/b


    plt.scatter([-13.2755786023,-13.2755786023],[-4.9716170082,4.9716170082],label='BG Short Period',c='lightblue')
    plt.scatter([-0.07791825611944855,-0.07791825611944855],[-1.2155409399418993,1.2155409399418993],label = 'BG Phugoid',c = 'orange')
    plt.scatter([-73.15940944802125,-73.15940944802125],[0,0],label = 'BG Roll Mode',c = 'lightgreen')
    plt.scatter([0.49542512183772963,0.49542512183772963],[0,0],label = 'BG Spiral Mode',c = 'lightcoral')
    plt.scatter([-3.2384848352138427,-3.2384848352138427],[-4.264866280273775,4.264866280273775],label = 'BG Dutch Roll',c = 'violet')
    plt.scatter([sp_eigval1.real,sp_eigval2.real],[sp_eigval1.imag,sp_eigval2.imag],label='Short Period Mode',c='darkblue')
    plt.scatter([ph_eigval1.real,ph_eigval2.real],[ph_eigval1.imag,ph_eigval2.imag],label='Phugoid Mode',c = 'darkorange')
    plt.scatter([roll_eigval1.real],[roll_eigval1.imag],label = 'Roll Mode',c = 'darkgreen')
    plt.scatter([spiral_eigval1.real],[spiral_eigval1.imag],label = 'Spiral Mode',c = 'darkred')
    plt.scatter([dutch_eigval1.real,dutch_eigval2.real],[dutch_eigval1.imag,dutch_eigval2.imag],label = 'Dutch Roll Mode',c = 'darkviolet')
    plt.legend()
    plt.plot([0,0],[6,-6],'--',)
    plt.title('Dimensional Eigenvalues')
    plt.xlabel('Real')
    plt.ylabel('Imaginary')
    plt.savefig('eigenvalues.png')
    plt.show()

    plt.scatter([-13.2755786023,-13.2755786023],[-4.9716170082,4.9716170082],label='BG Short Period',c='lightblue')
    plt.scatter([-0.07791825611944855,-0.07791825611944855],[-1.2155409399418993,1.2155409399418993],label = 'BG Phugoid',c = 'orange')
    plt.scatter([0.49542512183772963,0.49542512183772963],[0,0],label = 'BG Spiral Mode',c = 'lightcoral')
    plt.scatter([-3.2384848352138427,-3.2384848352138427],[-4.264866280273775,4.264866280273775],label = 'BG Dutch Roll',c = 'violet')
    plt.scatter([sp_eigval1.real,sp_eigval2.real],[sp_eigval1.imag,sp_eigval2.imag],label='Short Period Mode',c='darkblue')
    plt.scatter([ph_eigval1.real,ph_eigval2.real],[ph_eigval1.imag,ph_eigval2.imag],label='Phugoid Mode',c = 'darkorange')
    plt.scatter([spiral_eigval1.real],[spiral_eigval1.imag],label = 'Spiral Mode',c = 'darkred')
    plt.scatter([dutch_eigval1.real,dutch_eigval2.real],[dutch_eigval1.imag,dutch_eigval2.imag],label = 'Dutch Roll Mode',c = 'darkviolet')
    plt.legend()
    plt.plot([0,0],[6,-6],'--',)
    plt.title('Zoom-in of Dimensional Eigenvalues')
    plt.xlabel('Real')
    plt.ylabel('Imaginary')
    plt.savefig('zoomeigenvalues.png')
    plt.show()




    CW = W/(0.5*rho*Vo**2*Sw)
    CAP = sp_undamp_freq.real**2/(CL_a/CW)
    print('The CAP of this airplane is {0:10.7f} s^(-2)'.format(CAP))
    if (sp_damp_ratio >0.3) & (sp_damp_ratio<2.00) & (CAP > 0.038) & (CAP < 3.6):
        print('The aircraft is Level 1 for Short Period.')
    elif (sp_damp_ratio > 0.2) & (sp_damp_ratio<2.00) &  (CAP > 0.085) & (CAP < 10):
        print('The aircraft is Level 2 for Short Period.')
    elif (sp_damp_ratio > 0.15) & (CAP>10):
        print('The aircraft is Level 3 for Short Period.')
    else:
        print('The aircraft is Level 4 for Short Period.')


    if (ph_damp_ratio>0.04):
        print('The aircraft is Level 1 for Phugoid Mode.')
    elif (ph_damp_ratio>0):
        print('The aircraft is Level 2 for Phugoid Mode.')
    elif (ph_damp_time > 55):
        print('The aircraft is Level 3 for Phugoid Mode.')
    else:
        print('The aircraft is Level 4 for Phugoid Mode.')

    roll_time = 1/roll_sig

    if (roll_time < 1.4):
        print('The aircraft is Level 1 for Roll Mode.')
    elif (roll_time < 3):
        print('The aircraft is Level 2 for Roll Mode.')
    elif (roll_time < 10):
        print('The aircraft is Level 3 for Roll Mode.')
    else: 
        print('The aircraft is Level 4 for Roll Mode.')


    if (spiral_time > 20):
        print('The aircraft is Level 1 for Spiral Mode.')
    elif (spiral_time > 12):
        print('The aircraft is Level 2 for Spiral Mode.')
    elif (spiral_time > 4):
        print('The aircraft is Level 3 for Spiral Mode.')
    else:
        print('The aircraft is Level 4 for Spiral Mode.')

    product = dutch_roll_damp_ratio*dutch_roll_ud_freq
    if (dutch_roll_damp_ratio > 0.08) & (dutch_roll_ud_freq > 0.4) & (product > 0.15):
        print('The aircraft is Level 1 for Dutch Roll Mode.')
    elif (dutch_roll_damp_ratio > 0.02) & (dutch_roll_ud_freq > 0.4) & (product > 0.05):
        print('The aircraft is Level 2 for Dutch Roll Mode.')
    elif (dutch_roll_damp_ratio > 0) & (dutch_roll_ud_freq > 0.4):
        print('The aircraft is Level 3 for Dutch Roll Mode.')
    else:
        print('The aircraft is Level 4 for Dutch Roll Mode.')

