# locus fixer prog



import numpy as np
import glider_funcs as gfs

    
    
####---------------------------------------------------------------------####
#                                         Create Data Matrix, cases, and fill
# Columns:
# 0: aoa,  1: CL,  2: CD,  3: Cm,  4: CL_a,  5: Cm_a,  6: SM,  7: Cl_b,  8: Cn_b
aoa = np.linspace(-14,14,15)
data_mat = np.zeros((aoa.size,9))

data_mat = np.array(([-14,-0.6211,	0.0309,	0.46],\
        [-12,-0.4325,0.0224,0.3942],\
        [-10,-0.2441,	0.0158,	0.3287],\
        [-8,-0.0553,	0.012,	0.2634],\
        [-6,0.1339,	0.011,	0.198],\
        [-4,0.3238,	0.0127,	0.1322],\
        [-2,0.5146,	0.0172,	0.0658],\
        [0,0.7059,	0.0246,	0],\
        [2,0.8993,	0.0348,	-0.0701],\
        [4,1.0936,	0.0481,	-0.1402],\
        [6,1.2893,	0.0643,	-0.212],\
        [8,1.4866,	0.0837,	-0.2857],\
        [10,1.6856,	0.1063,	-0.3615],\
        [12,1.8864,	0.1321,	-0.4395],\
        [14,2.0891,	0.1615,	-0.52]))


#            Then Create Functions Matrix, use data, call functions, and fill
# Columns:
# 0: aoa,  1: L/D,  2: Vel,  3: Drag,  4: Pwr,  5: GR,  6: Sr
func_mat = np.zeros((aoa.size,7))
#                            Then Create Aero Calc Matrix, use data, and fill
# Columns:
# 0: aoa,  1: CN,  2: CA,  3: Cm,a,  4: CN,a,  5: CA,a,
# 6: Cm,a,a,  7: CN,a,a,  8: CA,a,a,  9: x,  10: y
aero_mat = np.zeros((aoa.size,11))

for i in range(aoa.size):
    # state = {
    #     "type" : "aerodynamic",
    #     "velocity" : params["vel"],
    #     "alpha" : aoa[i],
    #     "beta" : 0.0
    # }
    # my_scene.set_aircraft_state(state = state, aircraft = "glider_2.0")
    # FM_results = my_scene.solve_forces(dimensional=False, non_dimensional=True)
    # derivs = my_scene.derivatives()
    
    # # Unpackage and organize the required values.
    # # 0: aoa,  1: CL,  2: CD,  3: Cm,  4: CL_a,  5: Cm_a,  6: SM,  7: Cl_b,  8: Cn_b
    # data_mat[i][0] = aoa[i]
    # data_mat[i][1] = FM_results["glider_2.0"]["total"]["CL"]
    # data_mat[i][2] = FM_results["glider_2.0"]["total"]["CD"]
    # data_mat[i][3] = FM_results["glider_2.0"]["total"]["Cm"]
    # data_mat[i][4] = derivs["glider_2.0"]["stability"]["CL,a"]
    # data_mat[i][5] = derivs["glider_2.0"]["stability"]["Cm,a"]
    # data_mat[i][6] = derivs["glider_2.0"]["stability"]["%_static_margin"]
    # data_mat[i][7] = derivs["glider_2.0"]["stability"]["Cl,b"]
    # data_mat[i][8] = derivs["glider_2.0"]["stability"]["Cn,b"]

    
    
    # Run a few calcs
    # 0: aoa,  1: L/D,  2: Vel,  3: Drag,  4: Pwr,  5: GR,  6: Sr
    # func_mat[i][0] = aoa[i]
    # func_mat[i][1] = data_mat[i][1] / data_mat[i][2]
    # func_mat[i][2] = gfs.C_L_func(params["rho"],params["Sw"],C_L = data_mat[i][1], W= params["W"])
    # func_mat[i][3] = gfs.drag_func(data_mat[i][2],func_mat[i][2],params["rho"],params["Sw"])
    # func_mat[i][4] = func_mat[i][2] * func_mat[i][3]
    # func_mat[i][5] = data_mat[i][1] / data_mat[i][2]
    # func_mat[i][6] = func_mat[i][4] / params["W"]

    # Aero Center Calcs
    # 0: aoa,  1: CN,  2: CA,  3: Cm,a,  4: CN,a,  5: CA,a,
    # 6: Cm,a,a,  7: CN,a,a,  8: CA,a,a,  9: x,  10: y
    aero_mat[i][0] = aoa[i]
    aero_mat[i][1] = data_mat[i][1]*np.cos(aoa[i]*(np.pi/180)) + data_mat[i][2]*np.sin(aoa[i]*(np.pi/180))
    aero_mat[i][2] = data_mat[i][2]*np.cos(aoa[i]*(np.pi/180)) - data_mat[i][1]*np.sin(aoa[i]*(np.pi/180))

for i in range(aoa.size):
    if ((i == 0) or (i == aoa.size-1)):
        aero_mat[i][3] = 0
        aero_mat[i][4] = 0
        aero_mat[i][5] = 0
        aero_mat[i][6] = 0
        aero_mat[i][7] = 0
        aero_mat[i][8] = 0
        aero_mat[i][9] = 0
        aero_mat[i][10] = 0
    else:
        aero_mat[i][3] = (data_mat[i+1][3]-data_mat[i-1][3]) /(2*2*np.pi/180)
        aero_mat[i][4] = (aero_mat[i+1][1]-aero_mat[i-1][1]) /(2*2*np.pi/180)
        aero_mat[i][5] = (aero_mat[i+1][2]-aero_mat[i-1][2]) /(2*2*np.pi/180)
        aero_mat[i][6] = (data_mat[i+1][3]-2*data_mat[i][3]+data_mat[i-1][3]) /((2*np.pi/180)**2)
        aero_mat[i][7] = (aero_mat[i+1][1]-2*aero_mat[i][1]+aero_mat[i-1][1]) /((2*np.pi/180)**2)
        aero_mat[i][8] = (aero_mat[i+1][2]-2*aero_mat[i][2]+aero_mat[i-1][2]) /((2*np.pi/180)**2)
        bot = (aero_mat[i][4]*aero_mat[i][8] - aero_mat[i][5]*aero_mat[i][7])
        aero_mat[i][9] = ( (aero_mat[i][5]*aero_mat[i][6] - aero_mat[i][3]*aero_mat[i][8]) /bot) * 0.75
        aero_mat[i][10] = ( (aero_mat[i][4]*aero_mat[i][6] - aero_mat[i][3]*aero_mat[i][7]) /bot) * 0.75
        #### =(I103*J103-G103*L103)/(H103*L103-I103*K103)*$H$99
        #### =(H105*J105-G105*K105)/(H105*L105-I105*K105)*$H$99
        #### I = CAa, J=Cmaa, G=Cma, L=CAaa, H=CNa, K=CNaa


print(data_mat)
print(func_mat)
print(aero_mat)

np.savetxt("locus.csv", aero_mat, delimiter=",")