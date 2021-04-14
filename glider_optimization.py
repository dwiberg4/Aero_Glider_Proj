import machupX as MX
import numpy as np
import json
from scipy.optimize import minimize
import matplotlib.pyplot as plt




def trim_func(x,Cl,static_margin, Cn_beta, Cl_beta):
    with open('glider_2.0.json') as f:
        airplane = json.load(f)
    f.close()
    with open('scene.json') as f:
        scene = json.load(f)
    f.close()
    wing_twist = 2.5        # angle to twist wing down from mounting angle.
    airplane['CG'][0] = x[0]
    airplane['wings']['main_wing']['twist'] = [[0.0, x[1],
                                                1.0, x[1]-wing_twist]]
    airplane['wings']['outside_wings']['twist'] =  [[0.0, 0,
                                                    1.0, 0]]
    airplane['wings']['H_Stab']['connect_to']['dx'] = x[2]
    airplane['wings']['V_Stab']['connect_to']['dx'] = x[3]
    airplane['wings']['main_wing']['dihedral'] = x[4]
    with open('glider_2.0.json', 'w') as fp:
        json.dump(airplane,fp,indent = 4)
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
    trim(0.35, 0.15, -0.03)
    import req_check
    req_check.check_requirements('glider_2.0.json')
    import a10_a13_wrapper
    a10_a13_wrapper
    import dyna_matheson
    dyna_matheson


    # import dynamic
    # import a10_a13_wrapper

    # statmar = np.linspace(0.1,0.5,20)
    # n_beta = np.linspace(-0.25,0.25,20)
    # l_beta = np.linspace(-0.25,0.25,20)
    # short_period_data = []
    # phugoid_data = []
    # roll_data = []
    # spiral_data = []
    # dutch_roll_data = []
    # for i in statmar:
    #     trim(i,0.1,-0.05)
    #     a10_a13_wrapper
    #     lat,long = dynamic.dyna()
    #     short_period_data.append(np.abs(long[0]))
    #     phugoid_data.append(np.abs(long[2]))
    #     roll_data.append(np.abs(lat[0]))
    #     spiral_data.append(np.abs(lat[1]))
    #     dutch_roll_data.append(np.abs(lat[4]))
    # plt.scatter(statmar,short_period_data,label = 'Short Period Magnitude')
    # plt.xlabel('Static Margin')
    # plt.ylabel('Magnitude of Eigenvalue')
    # plt.title('Static Margin versus Eigenvalues')
    # plt.legend()
    # plt.savefig('StaticMargin_SP.png',dpi=800)
    # plt.show()
    # plt.scatter(statmar,phugoid_data,label = 'Phugoid Period Magnitude')
    # plt.xlabel('Static Margin')
    # plt.ylabel('Magnitude of Eigenvalue')
    # plt.title('Static Margin versus Eigenvalues')
    # plt.legend()
    # plt.savefig('StaticMargin_Ph.png',dpi=800)
    # plt.show()
    # plt.scatter(statmar,roll_data,label = 'Roll Mode Magnitude')
    # plt.xlabel('Static Margin')
    # plt.ylabel('Magnitude of Eigenvalue')
    # plt.title('Static Margin versus Eigenvalues')
    # plt.legend()
    # plt.savefig('StaticMargin_Roll.png',dpi=800)
    # plt.show()
    # plt.scatter(statmar,spiral_data,label='Spiral Mode Magnitude')
    # plt.xlabel('Static Margin')
    # plt.ylabel('Magnitude of Eigenvalue')
    # plt.title('Static Margin versus Eigenvalues')
    # plt.legend()
    # plt.savefig('StaticMargin_Spiral.png',dpi=800)
    # plt.show()
    # plt.scatter(statmar,dutch_roll_data,label = 'Dutch Roll Magnitude')
    # plt.xlabel('Static Margin')
    # plt.ylabel('Magnitude of Eigenvalue')
    # plt.title('Static Margin versus Eigenvalues')
    # plt.legend()
    # plt.savefig('StaticMargin_Dutch.png',dpi=800)
    # plt.show()


    # short_period_data = []
    # phugoid_data = []
    # roll_data = []
    # spiral_data = []
    # dutch_roll_data = []
    # for i in n_beta:
    #     trim(0.35,i,-0.05)
    #     a10_a13_wrapper
    #     lat,long = dynamic.dyna()
    #     short_period_data.append(np.abs(long[0]))
    #     phugoid_data.append(np.abs(long[2]))
    #     roll_data.append(np.abs(lat[0]))
    #     spiral_data.append(np.abs(lat[1]))
    #     dutch_roll_data.append(np.abs(lat[4]))
    # plt.scatter(n_beta,short_period_data,label = 'Short Period Magnitude')
    # plt.scatter(n_beta,phugoid_data,label = 'Phugoid Period Magnitude')
    # plt.scatter(n_beta,roll_data,label = 'Roll Mode Magnitude')
    # plt.scatter(n_beta,spiral_data,label='Spiral Mode Magnitude')
    # plt.scatter(n_beta,dutch_roll_data,label = 'Dutch Roll Magnitude')
    # plt.xlabel('Cn_b')
    # plt.ylabel('Magnitude of Eigenvalue')
    # plt.title('Cn_b versus Eigenvalues')
    # plt.legend()
    # plt.savefig('Cnb.png',dpi=800)
    # plt.show()


    # short_period_data = []
    # phugoid_data = []
    # roll_data = []
    # spiral_data = []
    # dutch_roll_data = []
    # for i in l_beta:
    #     trim(0.35,0.1,i)
    #     a10_a13_wrapper
    #     lat,long = dynamic.dyna()
    #     short_period_data.append(np.abs(long[0]))
    #     phugoid_data.append(np.abs(long[2]))
    #     roll_data.append(np.abs(lat[0]))
    #     spiral_data.append(np.abs(lat[1]))
    #     dutch_roll_data.append(np.abs(lat[4]))
    # plt.scatter(n_beta,short_period_data,label = 'Short Period Magnitude')
    # plt.scatter(n_beta,phugoid_data,label = 'Phugoid Period Magnitude')
    # plt.scatter(n_beta,roll_data,label = 'Roll Mode Magnitude')
    # plt.scatter(n_beta,spiral_data,label='Spiral Mode Magnitude')
    # plt.scatter(n_beta,dutch_roll_data,label = 'Dutch Roll Magnitude')
    # plt.xlabel('Cl_b')
    # plt.ylabel('Magnitude of Eigenvalue')
    # plt.title('Cl_b versus Eigenvalues')
    # plt.legend()
    # plt.savefig('Clb.png',dpi=800)
    # plt.show()
    # # To use, define everything in the JSON file, then define your trim(static_margin, Cn_beta, Cl_beta)
    # # Pressing run will update the CGx, H and V Stab location and dihedral inside of the JSON file.


