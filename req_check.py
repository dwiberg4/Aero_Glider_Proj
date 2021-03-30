import machupX as MX
import numpy as np
import json
from scipy.optimize import minimize

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



def check_requirements(filename):
    with open(filename) as f:
        airplane = json.load(f)
    f.close()
    total_volume = 0
    total_span = 0
    main_section = 0
    main_chord_average = 0
    wing_weight = 0
    area = 0
    for wing in airplane['wings']:
        double = 0
        if airplane['wings'][wing]['side'] == 'both':
            double = 2
        else:
            double = 1
        airfoil = airplane['wings'][wing]['airfoil']
        if airfoil == 'NACA_2412':
            airfoil_func = area2412
        if airfoil == 'NACA_0012':
            airfoil_func = area0012
        chord_average = 0
        i = 0
        for chord in airplane['wings'][wing]['chord']:
            chord_average +=chord[1]
            i +=1
        chord_average /=i
        avg_area = airfoil_func(chord_average)
        volume = airplane['wings'][wing]['semispan']*avg_area*double
        if airplane['wings'][wing]['is_main']:
            total_span += double*airplane['wings'][wing]['semispan']
            main_chord_average +=chord_average
            main_section +=1
        wing_weight +=0.7*volume
        area += chord_average*double*airplane['wings'][wing]['semispan']
    main_chord = main_chord_average/main_section
    Wf = 0.018408*27+wing_weight
    W = airplane['weight']
    aspect_ratio = total_span**2/(total_span*main_chord)
    rhs = (16/3*234)/(1 + Wf*total_span/(W-Wf))
    structure = False
    if aspect_ratio<=rhs:
        print('This airplane is structurally sound')
        print('The aspect ratio can be increased by {0:.3f} %.'.format(rhs/aspect_ratio))
    if area<=9:
        print('This airplane is less than 9 square feet.')
        print('There is {0:.3f} square foot area left.'.format(9-area))


if __name__ == '__main__':
    check_requirements('glider_2.0.json')


        