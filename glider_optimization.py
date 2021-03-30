import machupX as MX
import numpy as np
import json
from scipy.optimize import minimize

def trim_func(x,Cl,static_margin, Cn_beta, Cl_beta):
    with open('glider_2.0.json') as f:
        airplane = json.load(f)
    f.close()
    with open('scene.json') as f:
        scene = json.load(f)
    f.close()
    airplane['CG'] = [x[0],0.0,0.0]
    airplane['wings']['main_wing']['twist'] = [[0.0, x[1],
                                                1.0, x[1]]]
    airplane['wings']['outside_wings']['twist'] =  [[0.0, x[1],
                                                    1.0, x[1]]]
    airplane['wings']['H_Stab']['connect_to']['dx'] = x[2]
    airplane['wings']['V_Stab']['connect_to']['dx'] = x[3]
    airplane['wings']['outside_wings']['dihedral'] = x[4]
    weight = airplane['weight']
    #scene['scene']['aircraft']['glider_2.0']['state']['velocity'] = np.sqrt(10/(0.5*weight/32.217))
    with open('glider_2.0.json', 'w') as fp:
        json.dump(airplane,fp,indent = 4)
    fp.close()
    with open('scene.json', 'w') as fp:
        json.dump(scene,fp,indent = 4)
    fp.close()
    my_scene = MX.Scene('C:/Users/Matheson Price/Glider_Project/scene.json')
    FM_results = my_scene.solve_forces(non_dimensional = True)
    derivs = my_scene.derivatives()['glider_2.0']
    residual = FM_results['glider_2.0']['total']['Cm']**2 + (derivs['stability']['%_static_margin']/100-static_margin)**2+\
        + (derivs['stability']['Cn,b'] - Cn_beta)**2 + (derivs['stability']['Cl,b']-Cl_beta)**2 + (FM_results['glider_2.0']['total']['CL']-Cl)**2
    print('The moment coefficient is: ',FM_results['glider_2.0']['total']['Cm'])
    print('The static margin is: ',derivs['stability']['%_static_margin']/100)
    print('Cn,b is: ',derivs['stability']['Cn,b'])
    print('Cl,b is: ',derivs['stability']['Cl,b'])
    print('Cl is: ',FM_results['glider_2.0']['total']['CL'])
    print('L/D is: ',FM_results['glider_2.0']['total']['CL']/FM_results['glider_2.0']['total']['CD'])
    return residual

def trim(Cl,static_margin, Cn_beta, Cl_beta):
    x = [
        -0.141,              # CGx                                   0
        3.2,                # Mounting Angle fo the main wing       1
        -2.825,             # Location of Horizontal Stabilizer     2
        -3.161,             # Location of Vertical Stabilizer       3
        5,               # Main Wing Dihedral                    4
    ]
    constraints = ({'type':'ineq','fun':lambda x: x[4]})
    args = (Cl,static_margin, Cn_beta, Cl_beta)
    res = minimize(trim_func, x, constraints = constraints,args = args)
    print(res)

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


def structural(x):
    wing_volume = area2412(x[2])*x[0]       # volume of wing
    emp_volume = area0012(x[1])*x[0]/3      # Span of horiz and vert are 1/3 span of wing
    foam_weight = (wing_volume + emp_volume) * 0.7
    aspect_ratio = x[0]**2/(x[0]*x[1])
    wf = x[4] - foam_weight
    rhs = 16/3*243/(1+wf*x[0]/(x[2]-wf))
    return rhs-aspect_ratio

def square9(x):
    # Check that all the foam elements are within the 9 square feet
    wing_area = x[0]*x[1]
    emp_area = x[0]*x[1]*3
    return 9-(wing_area+emp_area)

def weight(x):
    wing_volume = area2412(x[1])*x[0]       # volume of wing
    emp_volume = area0012(x[1])*x[0]/3      # Span of horiz and vert are 1/3 span of wing
    foam_weight = (wing_volume + emp_volume) * 0.7
    dowel_weight = (0.75/24)**2*np.pi*6*27

if __name__ == '__main__':
    trim(0.5, 0.35, 0.1, -0.05)
    # To use, define everything in the JSON file, then define your trim(static_margin, Cn_beta, Cl_beta)
    # Pressing run will update the CGx, H and V Stab location and dihedral inside of the JSON file.


