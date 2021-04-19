# Program to solve for the eigenvalues of a system of equations

# Can recieve either 1 Matrix to solve the Special Eigenproblem
# OR can recieve 2 Matrices to solve generally

import numpy as np 
import math
import cmath as cm
from scipy.linalg import eig

# Eigensolver that takes for arguments either 1 or 2 matrices
# If given a single matrix, the function will assume and run 
#   Special Eigenproblem. If given 2 matrices, the function
#   will assume and run the General Eigenproblem
# kwargs: char
#   Will return a dictionary of the behavior characteristics
def eig_solve(Vo,b,c,*args, **kwargs):
    if (len(args)==2):
        #print("\n\tWe will now begin the Generalize Eigensolver routine\n")
        C = np.matmul(np.linalg.inv(args[1]),args[0])
        eigvals, eigvecs = eig(C)
    else:
        #print("\n\tWe will now begin the Special Eignsolver routine\n")
        C = args[0]
        eigvals, eigvecs = eig(C)

    if "title" in kwargs and kwargs["title"] == "Longitudinal":
        t_dimless = (c/(2*Vo))
        dim_eigs = eigvals / t_dimless
    elif "title" in kwargs and kwargs["title"] == "Lateral":
        t_dimless = (b/(2*Vo))
        dim_eigs = eigvals / t_dimless
    else:
        print("\nNo type specified. What Non-Dimensional Time do you want?")
        dim_eigs = np.zeros((6,6))
    
    #print(eigvecs)
    #for i in range(eigvals.size):   
        #print("The %dth index eigenvalue is: %f" %(i,eigvals[i]))
        #print("\tThe corresponding eigenvector is:", eigvecs[:,i])
    
    if "char" in kwargs and kwargs["char"]:
        if "title" in kwargs and kwargs["title"] == "Longitudinal":
            # for Baseline: which = [(2,3),(4,5)]
            which = [(2,3),(4,5)]
        elif "title" in kwargs and kwargs["title"] == "Lateral":
            which = [(3,4)]
        vals = beh_char(eigvals,eigvecs,t_dimless,which)
    else:
        vals = 0

    if "file" in kwargs and "title" in kwargs and kwargs["file"]:
        res_out(args,C,eigvals,eigvecs,vals,dim_eigs,kwargs["title"])

    # #_________0______1_______2_____3_______4____5_____6_____7_______8______9_____10__
    # vals = [w_n_d, sigmas, taus, halves, nines,dubs, w_n, zetas, periods, amps, phas]
    return (eigvals,eigvecs,vals,dim_eigs)
    
    
def beh_char(eigvals,eigvecs,t,which):
    # Recieves vector of eigvals; returns all the behavior characteristics
    w_n_d = damp_nat_freq(eigvals,t)
    #print("\n\t\tThe Damped Natural Frequencies are: \n",w_n_d)
    sigmas = damp_rate(eigvals,t)
    #print("\n\t\tThe Damping Rates are: \n\n",sigmas)
    taus = tau(sigmas)
    #print("\n\t\tThe Time Constants are: \n\n",taus)
    halves = half(sigmas)
    #print("\n\t\tThe Times to Half are: \n\n", halves)
    nines = ninenine(sigmas)
    #print("\n\t\tThe Times to 99'%' are: \n\n", nines)
    dubs = double(sigmas)
    #print("\n\t\tThe Double Times are: \n\n", dubs)
    w_n = nat_freq(eigvals,t,which)
    #print("\n\t\tThe Natural frequencies for the eigenvalue pairs are: \n\n",w_n)
    zetas = damp_ratio(eigvals,which)
    #print("\n\t\tThe Damping Ratios for the eigenvalue pairs are: \n\n",zetas)

    periods = period_func(w_n_d)
    amps = amp_func(eigvecs)
    phas = phase_func(eigvecs)
    #_________0______1_______2_____3_______4____5_____6_____7_______8______9_____10__
    vals = [w_n_d, sigmas, taus, halves, nines,dubs, w_n, zetas, periods, amps, phas]
    return vals


# Results File writer Function
def res_out(args,C,eigvals,eigvecs,vals,dim_eigs,title):
    name = title + "_Eig_sol.txt"
    with open(name, 'w') as resout:
        if (len(args)==2):
            resout.write('A Matrix Values: \n\n')
            for i in range(args[0].shape[0]):
                for j in range(args[0].shape[0]):
                    str_out = '{:>20.12f}'.format(args[0][i,j])
                    resout.write(str_out)
                resout.write('\n')
            resout.write('\n')
            resout.write('B Matrix Values: \n\n')
            for i in range(args[1].shape[0]):
                for j in range(args[1].shape[0]):
                    str_out = '{:>20.12f}'.format(args[1][i,j])
                    resout.write(str_out)
                resout.write('\n')
            resout.write('\n')
            resout.write('C Matrix Values: \n\n')
            for i in range(C.shape[0]):
                for j in range(C.shape[0]):
                    str_out = '{:>20.12f}'.format(C[i,j])
                    resout.write(str_out)
                resout.write('\n')
            resout.write('\n')
        else:
            resout.write('C Matrix Values: \n\n')
            for i in range(C.shape[0]):
                for j in range(C.shape[0]):
                    str_out = '{:>20.12f}'.format(C[i,j])
                    resout.write(str_out)
                resout.write('\n')
            resout.write('\n')

        resout.write('Eigenvalues: \n\n')
        for i in range(eigvals.size):
            str_out = '{:>20.12f}'.format(eigvals[i])
            resout.write(str_out)
            resout.write('\n')
        resout.write('\n')

        resout.write('Eigenvectors: \n\n')
        for i in range(eigvecs.shape[0]):
            for j in range(eigvecs.shape[1]):
                str_out = '{:>20.6f}'.format(eigvecs[i,j])
                resout.write(str_out)
            resout.write('\n')
        resout.write('\n')

        resout.write('DIMENSIONAL Eigenvalues: \n\n')
        for i in range(dim_eigs.size):
            str_out = '{:>20.12f}'.format(dim_eigs[i])
            resout.write(str_out)
            resout.write('\n')
        resout.write('\n')

        resout.write('Amplitudes: \n\n')
        for i in range(vals[9].shape[0]):
            for j in range(vals[9].shape[1]):
                str_out = '{:>20.12f}'.format(vals[9][i,j])
                resout.write(str_out)
            resout.write('\n')
        resout.write('\n')

        resout.write('Phases: \n\n')
        for i in range(vals[10].shape[0]):
            for j in range(vals[10].shape[1]):
                str_out = '{:>20.12f}'.format(vals[10][i,j])
                resout.write(str_out)
            resout.write('\n')
        resout.write('\n')

        resout.write('Damping Rates: \n\n')
        for i in range(vals[1].size):
            str_out = '{:>20.12f}'.format(vals[1][i])
            resout.write(str_out)
            resout.write('\n')
        resout.write('\n')

        resout.write('99% Damping Times (For Convergent Modes): \n\n')
        for i in range(vals[4].size):
            str_out = '{:>20.12f}'.format(vals[4][i])
            resout.write(str_out)
            resout.write('\n')
        resout.write('\n')

        resout.write('Doubling Times (For Divergent Modes): \n\n')
        for i in range(vals[5].size):
            str_out = '{:>20.12f}'.format(vals[5][i])
            resout.write(str_out)
            resout.write('\n')
        resout.write('\n')

        resout.write('Damped Natural Frequencies (For Oscillatory Systems): \n\n')
        for i in range(vals[0].size):
            str_out = '{:>20.12f}'.format(vals[0][i])
            resout.write(str_out)
            resout.write('\n')
        resout.write('\n')

        resout.write('Period (For Oscillatory Systems): \n\n')
        for i in range(vals[8].size):
            str_out = '{:>20.12f}'.format(vals[8][i])
            resout.write(str_out)
            resout.write('\n')
        resout.write('\n')

        resout.write('Damping Ratio (For Complex Eigenvalue Pairs): \n\n')
        for i in range(vals[7].size):
            str_out = '{:>20.12f}'.format(float(vals[7][i]))
            resout.write(str_out)
            resout.write('\n')
        resout.write('\n')

        resout.write('Undamped Natural Frequency (For Complex Eigenvalue Pairs): \n\n')
        for i in range(vals[6].size):
            str_out = '{:>20.12f}'.format(float(vals[6][i]))
            resout.write(str_out)
            resout.write('\n')
        resout.write('\n')
    

def damp_nat_freq(eigvals,t):
    # Recieves vector of eigvals; returns corresponding w_d vector
    return (abs(eigvals.imag)) / t

def damp_rate(eigvals,t):
    # Recieves vector of eigvals; returns corresponding sigma vector
    return (-(eigvals.real)) / t

def tau(sigmas):
    # Recieves vector of sigmas; returns corresponding tau vector
    return (1/sigmas)

def half(sigmas):
    # Recieves vector of sigmas; returns corresponding Time to Half vector
    return (-np.log(0.5)) / sigmas

def ninenine(sigmas):
    # Recieves vector of sigmas; returns corresponding 99% damping time vector
    return (-np.log(0.01)) / sigmas

def double(sigmas):
    # Recieves vector of sigmas; returns corresponding Double Time vector
    return (-np.log(2)) / sigmas 

def nat_freq(eigvals,t,which):
    # Recieves vector of eigvals; returns 1/2 size vector of natural
    #   frequencies of eigenvalue pairs
    temps = np.zeros((len(which),1),dtype = 'complex_')
    w_n = np.zeros((len(which),1))
    for i in range(len(which)):
        temps[i] = (cm.sqrt(eigvals[which[i][0]]*eigvals[which[i][1]])) / t
        w_n[i] = float(temps[i].real)
    return w_n

def damp_ratio(eigvals,which):
    # Recieves vector of eigvals; returns 1/2 size vector of damping
    # ratios for eigenvalue pairs
    temps = np.zeros((len(which),1),dtype = 'complex_')
    zeta = np.zeros((len(which),1))
    for i in range(len(which)):
        top = -(eigvals[which[i][0]] + eigvals[which[i][1]])
        temps[i] = top / (2* (cm.sqrt(eigvals[which[i][0]]*eigvals[which[i][1]])) )
        zeta[i] = float(temps[i].real)
    return zeta

def period_func(w_n_d):
    # Recieves vector of damped natural frequencies; returns 1/2 size
    # vector of damping ratios for eigenvalue pairs
    return ((2*np.pi) / w_n_d ) 

def amp_func(eigvecs):
    # Receives eigvecs; returns square matrix of amplitude for each 
    # component of each eigenvector
    amp = np.zeros((eigvecs.shape))
    for i in range(eigvecs.shape[0]):
        for j in range(eigvecs.shape[1]):
            amp[i,j] = np.sqrt( (eigvecs[i,j].real**2) + (eigvecs[i,j].imag**2) )
    return amp

def phase_func(eigvecs):
    # Receives eigvecs; returns square matrix of phase for each 
    # component of each eigenvector
    pha = np.zeros((eigvecs.shape))
    for i in range(eigvecs.shape[0]):
        for j in range(eigvecs.shape[1]):
            pha[i,j] = math.atan2( (eigvecs[i,j].real),(eigvecs[i,j].imag) )
    return pha





# # Manual workspace____
# # For hardcoding the values of the A and B matrices

# m1 = 20
# m2 = 20
# c1 = 30
# c2 = 15
# k1 = 2
# k2 = 100

# A = np.matrix([[-c1,0,(-k1-k2),k2], \
#     [0,-c2,k2,-k2],\
#     [1,0,0,0],\
#     [0,1,0,0]])
# # print(A)
# B = np.matrix([[m1,0,0,0],\
#     [0,m2,0,0],\
#     [0,0,1,0],\
#     [0,0,0,1]])
# # print(B)

# C = [[-1.5,0.0,-5.1,5.0],\
#  [0.0,-0.75,5.0,-5.0],\
#  [1.0,0.0,0.0,0.0],\
#  [0.0,1.0,0.0,0.0,]]



# eig_solve(A,B, char = True)
# #eig_solve(C)