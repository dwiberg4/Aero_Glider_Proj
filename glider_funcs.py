# Aero; Flight Mechanics Glider Proj

import numpy as np 
import matplotlib.pyplot as plt

# KE equation
def kin_energy(g,KE, **kwargs):
    if kwargs is not None and 'W' in kwargs:
        #print (kwargs['W'])
        W = kwargs['W']
        vel = np.sqrt((g/W)*2*KE)
        #print("The velocity is: ",vel)
        return vel
    elif kwargs is not None and 'V' in kwargs:
        V = kwargs['V']
        weight = (2*KE*g)/(V**2)
        #print("the weight is: ",weight)
        return weight
    else:
        print("No key words present in arguments")


# Foam strength/ tolerance equation
def foam_toler(J,W,W_f,b,Ra):
    top = (16/3)*J
    bot = 1 + ((W_f/(W-W_f))*b)
    func_Ra = top/bot
    if (Ra <= func_Ra):
        sound = True
    else:
        sound = False


# C_L equation
def C_L_func(rho,Sw, **kwargs):
    if 'W' in kwargs and 'V' in kwargs:
        W = kwargs['W']
        V = kwargs['V']
        C_L = W / (0.5*rho*(V**2)*Sw)
        return C_L
    elif 'C_L' in kwargs and 'V' in kwargs:
        C_L = kwargs['C_L']
        V = kwargs['V']
        W = (0.5*rho*(V**2)*Sw) * C_L
        return W
    elif 'C_L' in kwargs and 'W' in kwargs:
        C_L = kwargs['C_L']
        W = kwargs['W']
        V = np.sqrt( (2*W) / (rho*Sw*C_L) )
        return V
    else:
        print("No key words present in arguments")


# Drag Equation
def drag_func(CD,Vel,rho,Sw):
    D = CD * 0.5 * rho * (Vel**2) * Sw
    return D


# Oswald Efficieny Factor
def oef(C_D2,Ra):
    e = 1/ (C_D2* np.pi* Ra)
    return e


# Stall Airspeed Equation V_min
def V_min_func(CL_max,W,Sw,rho):
    a = np.sqrt(2/CL_max)
    b = np.sqrt((W/Sw)/rho)
    V_min = a*b
    return V_min


# Minimum Drag Airspeed; Glide Ratio Optimizing Speed
def V_MD_func(e,Ra,C_D0,W,Sw,rho):
    a1 = np.sqrt(2)
    a2 = (np.pi* e* Ra* C_D0)**(0.25)
    b = np.sqrt((W/Sw)/rho)
    V_MD = (a1/a2) * b
    return V_MD


# Minimum Power Required Airspeed
def V_MDV_func(e,Ra,C_D1,C_D0,W,Sw,rho):
    a1 = (np.pi* e* Ra* C_D1)
    a2 = (np.pi* e* Ra* C_D1) **2
    a3 = 12* e* Ra* C_D0
    b = np.sqrt((W/Sw)/rho)
    V_MDV = (2*b) / (np.sqrt(a1 + (np.sqrt(a2 + a3))))
    return V_MDV



# Ballast Weight Equation
def W_b_func(rho_b,bdiam_inner,bdiam_outer, **kwargs):
    cs_area = ( (np.pi*(bdiam_outer**2))/4 ) - ( (np.pi*(bdiam_inner**2))/4)
    if kwargs is not None and 'W_b' in kwargs:
        W_b = kwargs['W_b']
        bthick = (W_b/rho_b) / cs_area
        return bthick
    elif kwargs is not None and 'bthick' in kwargs:
        print("are we going in here?")
        bthick = kwargs['bthick']
        W_b = (bthick*cs_area) * rho_b
        return W_b
    else:
        print("No key words present in arguments")


# Dowel Weight Equation
def W_d_func(rho_d,ddiam, **kwargs):
    cs_area = (np.pi*(ddiam**2))/4
    if kwargs is not None and 'W_d' in kwargs:
        W_d = kwargs['W_d']
        dowel = W_d / (cs_area*rho_d)
        return dowel
    elif kwargs is not None and 'dowel' in kwargs:
        dowel = kwargs['dowel']
        W_d = cs_area * dowel * rho_d
        return W_d
    else:
        print("No key words present in arguments")


#  Fuselage Weight Equation
def W_f_func(**kwargs):
    if 'W_b' in kwargs and 'W_d' in kwargs:
        W_b = kwargs['W_b']
        W_d = kwargs['W_d']
        W_f = W_b + W_d
        return W_f
    elif 'W_f' in kwargs and 'W_d' in kwargs:
        W_f = kwargs['W_f']
        W_d = kwargs['W_d']
        W_b = W_f - W_d
        return W_b
    elif 'W_f' in kwargs and 'W_b' in kwargs:
        W_f = kwargs['W_f']
        W_b = kwargs['W_b']
        W_d = W_f - W_b
        return W_d
    else:
        print("No key words present in arguments")


def area2412(chord):
    # Calculates the area of the NACA2412 airfoil given a chord length
    x = np.arange(0,1.00001,0.00001)
    p = 0.4
    m = 0.02
    yc1 = m/p**2*(2*p*x - x**2)*np.where(x<=p,x,0)
    yc2 = m/(1-p)**2*(1-2*p + 2*p*x - x**2)*np.where(x>p,x,0)
    yc = yc1 + yc2
    dyc1 = 2*m/p**2*(p-x)*np.where(x<=p,x,0)
    dyc2 = 2*m/(1-p)**2*(p-x)*np.where(x>p,x,0)
    dyc = dyc1+dyc2
    yt = 0.12/0.2*(0.2969*x**0.5 - 0.126*x - 0.3516*x**2 + 0.2843*x**3 -0.1015*x**4)
    theta = np.arctan(dyc)
    xu = x - yt*np.sin(theta)
    yu = yc + yt*np.cos(theta)
    xl = x + yt*np.sin(theta)
    yl = yc - yt*np.cos(theta)
    area1 = np.trapz(yu*chord,xu*chord)
    area2 = abs(np.trapz(yl*chord,xl*chord))
    return area1 + area2

def area0012(chord):
    # Calculates the area of the NACA0012 airfoil given a chord length
    x = np.arange(0,1.00001,0.00001)
    # yc1 = m/p**2*(2*p*x - x**2)*np.where(x<=p,x,0)
    # yc2 = m/(1-p)**2*(1-2*p + 2*p*x - x**2)*np.where(x>p,x,0)
    yc = 0
    # dyc1 = 2*m/p**2*(p-x)*np.where(x<=p,x,0)
    # dyc2 = 2*m/(1-p)**2*(p-x)*np.where(x>p,x,0)
    dyc = 0
    yt = 0.12/0.2*(0.2969*x**0.5 - 0.126*x - 0.3516*x**2 + 0.2843*x**3 -0.1015*x**4)
    theta = np.arctan(dyc)
    xu = x - yt*np.sin(theta)
    yu = yc + yt*np.cos(theta)
    xl = x + yt*np.sin(theta)
    yl = yc - yt*np.cos(theta)
    area1 = np.trapz(yu*chord,xu*chord)
    area2 = abs(np.trapz(yl*chord,xl*chord))
    return area1 + area2

# Wing Weight Function
# def wing_weight(fthick,SA,rho_f):
#     W_w = fthick * SA * rho_f
#     return W_w


# Generic Plotter Function
def gen_plotter(x,y, **kwargs):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    if "y2" not in kwargs:
        ax.plot(x,y)
    if "y2" in kwargs:
        ax.plot(x,y,label=kwargs["set1_label"])
        ax.plot(x,kwargs["y2"],label=kwargs["set2_label"])
        ax.legend()
    if 'xlabel' in kwargs:
        ax.set_xlabel(kwargs['xlabel'])
    if 'ylabel' in kwargs:
        ax.set_ylabel(kwargs['ylabel'])
    if 'title' in kwargs:
        ax.set_title(kwargs['title'])
    plt.show()
    if 'save' in kwargs:
        fig.savefig(kwargs['save']) 







def main():
    vel = kin_energy(g,KE, W = 12)
    weight = kin_energy(g,KE, V = 15.7)
    W_w = wing_weight(fthick, SA, rho_f)
    print(W_w)
    W = 1.0
    W_d = W_d_func(rho_d,ddiam, dowel = dowel_allow)
    print("The weight of the dowel is: ", W_d)
    W_f = W - W_w
    
    W_b = W_f_func(W_f = W_f, W_d = W_d)
    print("The weight of the ballast is: ", W_b)

    bthick = W_b_func(rho_b, bdiam_inner, bdiam_outer, W_b = W_b)
    print("The thickness of the ballast is: %0.6f [ft]" % (bthick))
    print("The thickness of the ballast is: %0.6f [in]" % (bthick*12))

    # if bthick prescribed:
    bthick2 = 0.5/12
    W_b2 = W_b_func(rho_b, bdiam_inner, bdiam_outer, bthick = bthick2)
    print("The new W_b value is: ", W_b2)
    W_f2 = W_f_func(W_b = W_b2, W_d = W_d)
    print("The new W_f value is: ", W_f2)
    W2 = W_f2 + W_w
    print("The new W value is: ", W2)
    vel2 = kin_energy(g,KE, W = W2)
    print("The new Velocity is: ", vel2)
    CL2 = C_L_func(rho,Sw,V = vel2, W = W2)
    print("The new C_L values is: ",CL2)


#main()



