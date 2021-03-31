# Program to solve for the eigenvalues of a system of equations

# Can recieve either 1 Matrix to solve the Special Eigenproblem
# OR can recieve 2 Matrices to solve generally

import numpy as np 
from scipy.linalg import eig

# Eigensolver that takes for arguments either 1 or 2 matrices
# If given a single matrix, the function will assume and run 
#   Special Eigenproblem. If given 2 matrices, the function
#   will assume and run the General Eigenproblem
# kwargs: char
#   Will return a dictionary of the behavior characteristics
def eig_solve(*args, **kwargs):
    if (len(args)==2):
        print("\n\tWe will now begin the Generalize Eigensolver routine\n")
        C = np.matmul(np.linalg.inv(args[1]),args[0])
        eigvals, eigvecs = eig(C)
    else:
        print("\n\tWe will now begin the Special Eignsolver routine\n")
        C = args[0]
        eigvals, eigvecs = eig(C)

    for i in range(eigvals.size):
        print("The %dth index eigenvalue is: %f" %(i,eigvals[i]))
        print("\tThe corresponding eigenvector is:", eigvecs[:,i])
    
    if "char" in kwargs and kwargs["char"]:
        vals = beh_char(eigvals)
    else:
        vals = 0
    return (eigvals,eigvecs,vals)
    
    
def beh_char(eigvals):
    # Recieves vector of eigvals; returns all the behavior characteristics
    w_n_d = damp_nat_freq(eigvals)
    print("\n\t\tThe Damped Natural Frequencies are: ",w_n_d)
    sigmas = damp_rate(eigvals)
    print("\n\t\tThe Damping Rates are: ",sigmas)
    taus = tau(sigmas)
    print("\n\t\tThe Time Constants are: ",taus)
    halves = half(sigmas)
    print("\n\t\tThe Times to Half are: ", halves)
    nines = ninenine(sigmas)
    print("\n\t\tThe Times to 99'%' are: ", nines)
    dubs = double(sigmas)
    print("\n\t\tThe Double Times are: ", dubs)
    w_n = nat_freq(eigvals)
    print("\n\t\tThe Natural frequencies for the eigenvalue pairs are: ",w_n)
    zetas = damp_ratio(eigvals)
    print("\n\t\tThe Damping Ratios for the eigenvalue pairs are: ",zetas)

    vals = [w_n_d,sigmas,taus,halves,nines,dubs,w_n,zetas]
    return vals

def damp_nat_freq(eigvals):
    # Recieves vector of eigvals; returns corresponding w_d vector
    return abs(eigvals.imag)

def damp_rate(eigvals):
    # Recieves vector of eigvals; returns corresponding sigma vector
    return -(eigvals.real)

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

def nat_freq(eigvals):
    # Recieves vector of eigvals; returns 1/2 size vector of natural
    #   frequencies of eigenvalue pairs
    w_n = np.zeros(((int(eigvals.size/2)),1),dtype = 'complex_')
    for i in range(int(eigvals.size/2)):
        j = i*2
        w_n[i] = np.sqrt(eigvals[j]*eigvals[j+1])
    return w_n

def damp_ratio(eigvals):
    # Recieves vector of eigvals; returns 1/2 size vector of damping
    # ratios for eigenvalue pairs
    zeta = np.zeros(((int(eigvals.size/2)),1),dtype = 'complex_')
    for i in range(int(eigvals.size/2)):
        j = i*2
        top = -(eigvals[j] + eigvals[j+1])
        zeta[i] = top / (2* (np.sqrt(eigvals[j]*eigvals[j+1])) )
    return zeta





# Manual workspace____
# For hardcoding the values of the A and B matrices

m1 = 20
m2 = 20
c1 = 30
c2 = 15
k1 = 2
k2 = 100

A = np.matrix([[-c1,0,(-k1-k2),k2], \
    [0,-c2,k2,-k2],\
    [1,0,0,0],\
    [0,1,0,0]])
# print(A)
B = np.matrix([[m1,0,0,0],\
    [0,m2,0,0],\
    [0,0,1,0],\
    [0,0,0,1]])
# print(B)

C = [[-1.5,0.0,-5.1,5.0],\
 [0.0,-0.75,5.0,-5.0],\
 [1.0,0.0,0.0,0.0],\
 [0.0,1.0,0.0,0.0,]]



eig_solve(A,B, char = True)
#eig_solve(C)