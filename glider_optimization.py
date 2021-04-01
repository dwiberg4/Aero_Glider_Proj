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
                                                1.0, x[1]-2.5]]
    airplane['wings']['outside_wings']['twist'] =  [[0.0, 0,
                                                    1.0, 0]]
    airplane['wings']['H_Stab']['connect_to']['dx'] = x[2]
    airplane['wings']['V_Stab']['connect_to']['dx'] = x[3]
    airplane['wings']['main_wing']['dihedral'] = x[4]
    weight = airplane['weight']
    #scene['scene']['aircraft']['glider_2.0']['state']['velocity'] = np.sqrt(10/(0.5*weight/32.217))
    with open('glider_2.0.json', 'w') as fp:
        json.dump(airplane,fp,indent = 4)
    fp.close()
    with open('scene.json', 'w') as fp:
        json.dump(scene,fp,indent = 4)
    fp.close()
    my_scene = MX.Scene('scene.json')
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

def trim(static_margin, Cn_beta, Cl_beta):
    x = [
        -0.141,              # CGx                                   0
        3.2,                # Mounting Angle fo the main wing       1
        -2.825,             # Location of Horizontal Stabilizer     2
        -3.161,             # Location of Vertical Stabilizer       3
        5,               # Main Wing Dihedral                    4
    ]
    with open('glider_2.0.json') as f:
        plane = json.load(f)
    f.close()
    weight = plane['weight']
    velocity = np.sqrt((32.17/weight)*2*10)
    CL = weight/(0.5*0.0020482*velocity**2*plane['reference']['area'])
    constraints = ({'type':'ineq','fun':lambda x: x[4]})
    args = (CL,static_margin, Cn_beta, Cl_beta)
    res = minimize(trim_func, x, constraints = constraints,args = args)
    print(res)

if __name__ == '__main__':
    trim(0.35, 0.1, -0.05)
    # To use, define everything in the JSON file, then define your trim(static_margin, Cn_beta, Cl_beta)
    # Pressing run will update the CGx, H and V Stab location and dihedral inside of the JSON file.


