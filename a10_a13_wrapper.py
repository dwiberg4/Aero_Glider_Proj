# The A10 through A13 Wrapper script for the Glider Project
# -Static Analysis at this point

# A10. Compute the following as a function of alpha from −14degrees to +14degrees
# C_L
# C_D
# C_m
# L/D
# Velocity [ft/s] required for the lift to be equal to the weight at the alpha of interest
# Drag [lbf]
# Required Power [ft*lbf/s]
# Static Margin
# No-Wind Glide Ratio
# Sink Rate [ft/s]
# Pitch, Roll, and Yaw stability derivatives

# A11. Plot the above as a function of velocity from the data computed

# A12.Compute the three drag coefficients C_Do, C_D1, and C_D2 that describe the 
# drag of the entire aircraft as a function of lift.
# Compute the Oswald Efficiency for your aircraft from C_D2

# A13. Compute the stall speed (V_min), minimum drag airspeed (V_MD),
# and minimum power airspeed (V_MDV).

import numpy as np 
import machupX as MX
import json
import glider_funcs as gfs
import bal_bal_tool as bbt
import matplotlib.pyplot as plt
import pprint

def main(input_file,plot):
    ####---------------------------------------------------------------------####
    #                                                            Initialize scene
    my_scene = MX.Scene(input_file)
    my_scene.display_wireframe()


    ####---------------------------------------------------------------------####
    #                     Bring in a few values from the scene and aircraft jsons
    # Create the Params Dictionary:
    params = {}
    params["h"] = 100                         # 100 ft above the ground
    params["alpha"] = 0                       # Launch aoa is 0 deg.
    params["gamma"] = 0                       # Launch elevation angle is 0 deg.
    params["KE"] = 10                         # Contest constant Kinetic Energy is 10 ft-lbf
    params["g"] = 32.1741                     # Gravity in ft/s^2
    params["rho"] = 0.0020482                 # [slugs/ft^2] Atmospheric density
    params["SA_allow"] = 9                    # 9 ft^2 of allowable EPS foam area
    params["dowel_allow"] = 6                 # 6 ft dowel length; allowable x-range
    params["ddiam"] = 0.75/12                 # [ft] dowel diameter
    params["bdiam_inner"] = 0.75/12           # [ft] inner diam of ballast
    params["bdiam_outer"] = 1.47/12           # [ft] outer diam of ballast
    params["J"] = 234                         # [ft] Foam Structural parameter
    params["rho_f"] = 0.7                     # [lbm/ft^3] Foam density
    params["rho_d"] = 27                      # [lbm/ft^3] Wooden Dowel density
    params["rho_b"] = 1186.13                 # [lbm/ft^3] Ballast density
    params["fthick"] = 0.061079               # [ft] Back-calc'd thickness of foam per D.H.
    with open(input_file, "r") as f:
        scene = json.load(f)
        params["rho"] = scene["scene"]["atmosphere"]["rho"]
        params["aircraft_file"] = scene["scene"]["aircraft"]["glider_2.0"]["file"]
        params["vel"] = scene["scene"]["aircraft"]["glider_2.0"]["state"]["velocity"]
        with open(params["aircraft_file"], "r") as g:
            aircraft = json.load(g)
            params["CG"] = aircraft["CG"]
            params["W"] = aircraft["weight"]
            params["Sw"] = aircraft["reference"]["area"]
            params["CL_max"] = aircraft["airfoils"]["NACA_2412"]["CL_max"]
            mainwing_b = (aircraft["wings"]["main_wing"]["semispan"] * 2)
            outerwing_b = (aircraft["wings"]["outside_wings"]["semispan"] * 2)
            params["b"] = mainwing_b + outerwing_b
            m_c_root = (aircraft["wings"]["main_wing"]["chord"][0][1])
            m_c_tip = (aircraft["wings"]["main_wing"]["chord"][1][1])
            mainwing_c = (m_c_root + m_c_tip) / 2
            o_c_root = (aircraft["wings"]["outside_wings"]["chord"][0][1])
            o_c_tip = (aircraft["wings"]["outside_wings"]["chord"][1][1])
            outerwing_c = (o_c_root + o_c_tip) / 2
            params["wing_avg_c"] = (mainwing_c + outerwing_c) / 2
            if params["Sw"] != ((mainwing_b * mainwing_c) + (outerwing_b * outerwing_c)):
                print("UPDATE THE REFERENCE AREA TO MATCH")
            params["Sw"] = (mainwing_b * mainwing_c) + (outerwing_b * outerwing_c)
            hstab_b = (aircraft["wings"]["H_Stab"]["semispan"] * 2)
            h_c_root = (aircraft["wings"]["H_Stab"]["chord"][0][1])
            h_c_tip = (aircraft["wings"]["H_Stab"]["chord"][1][1])
            hstab_c = (h_c_root + h_c_tip) / 2
            Sw_H = hstab_b * hstab_c
            params["Sw_H"] = Sw_H
            params["hstab_b"] = hstab_b
            params["hstab_c"] = hstab_c
            vstab_b = (aircraft["wings"]["V_Stab"]["semispan"])
            v_c_root = (aircraft["wings"]["V_Stab"]["chord"][0][1])
            v_c_tip = (aircraft["wings"]["V_Stab"]["chord"][1][1])
            vstab_c = (v_c_root + v_c_tip) / 2
            Sw_V = vstab_b * vstab_c
            params["Sw_V"] = Sw_V
            params["vstab_b"] = vstab_b
            params["vstab_c"] = vstab_c
            params["x_hstab"] = aircraft["wings"]["H_Stab"]["connect_to"]["dx"]
            params["x_vstab"] = aircraft["wings"]["V_Stab"]["connect_to"]["dx"]
            params["CLa_2412"] = aircraft["airfoils"]["NACA_2412"]["CLa"]
            params["CLa_0012"] = aircraft["airfoils"]["NACA_0012"]["CLa"]

    params["Ra"] = ((mainwing_b + outerwing_b)**2) / params["Sw"]
    print("\nThe aspect ratio is: ", params["Ra"])

    ####---------------------------------------------------------------------####
    #                                         Create Data Matrix, cases, and fill
    # Columns:
    # 0: aoa,  1: CL,  2: CD,  3: Cm,  4: CL_a,  5: Cm_a,  6: SM,  7: Cl_b,  8: Cn_b
    aoa = np.linspace(-5,5,29)
    data_mat = np.zeros((aoa.size,9))
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
        state = {
            "type" : "aerodynamic",
            "velocity" : params["vel"],
            "alpha" : aoa[i],
            "beta" : 0.0
        }
        my_scene.set_aircraft_state(state = state, aircraft = "glider_2.0")
        FM_results = my_scene.solve_forces(dimensional=False, non_dimensional=True)
        derivs = my_scene.derivatives()
        
        # Unpackage and organize the required values.
        # 0: aoa,  1: CL,  2: CD,  3: Cm,  4: CL_a,  5: Cm_a,  6: SM,  7: Cl_b,  8: Cn_b
        data_mat[i][0] = aoa[i]
        data_mat[i][1] = FM_results["glider_2.0"]["total"]["CL"]
        data_mat[i][2] = FM_results["glider_2.0"]["total"]["CD"]
        data_mat[i][3] = FM_results["glider_2.0"]["total"]["Cm"]
        data_mat[i][4] = derivs["glider_2.0"]["stability"]["CL,a"]
        data_mat[i][5] = derivs["glider_2.0"]["stability"]["Cm,a"]
        data_mat[i][6] = derivs["glider_2.0"]["stability"]["%_static_margin"]
        data_mat[i][7] = derivs["glider_2.0"]["stability"]["Cl,b"]
        data_mat[i][8] = derivs["glider_2.0"]["stability"]["Cn,b"]
        
        # Run a few calcs
        # 0: aoa,  1: L/D,  2: Vel,  3: Drag,  4: Pwr,  5: GR,  6: Sr
        func_mat[i][0] = aoa[i]
        func_mat[i][1] = data_mat[i][1] / data_mat[i][2]
        func_mat[i][2] = gfs.C_L_func(params["rho"],params["Sw"],C_L = data_mat[i][1], W= params["W"])
        func_mat[i][3] = gfs.drag_func(data_mat[i][2],func_mat[i][2],params["rho"],params["Sw"])
        func_mat[i][4] = func_mat[i][2] * func_mat[i][3]
        func_mat[i][5] = data_mat[i][1] / data_mat[i][2]
        func_mat[i][6] = func_mat[i][4] / params["W"]

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
            aero_mat[i][9] = ((aero_mat[i][5]*aero_mat[i][6] - aero_mat[i][3]*aero_mat[i][8])/bot) * params["wing_avg_c"]
            aero_mat[i][10] = ((aero_mat[i][4]*aero_mat[i][6] - aero_mat[i][3]*aero_mat[i][7])/bot) * params["wing_avg_c"]
            np.savetxt("locus.csv", aero_mat, delimiter=",")
            #### =(I103*J103-G103*L103)/(H103*L103-I103*K103)*$H$99
            #### =(H105*J105-G105*K105)/(H105*L105-I105*K105)*$H$99
            #### I = CAa, J=Cmaa, G=Cma, L=CAaa, H=CNa, K=CNaa




    if plot == True:
        # Generate Plots
        title = "C_L as a Function of Velocity"
        gfs.gen_plotter(func_mat[:,2],data_mat[:,1],xlabel='Velocity',ylabel='CL',title=title,save='fig_CL.png')
        title = "C_D as a Function of Velocity"
        gfs.gen_plotter(func_mat[:,2],data_mat[:,2],xlabel='Velocity',ylabel='CD',title=title,save='fig_CD.png')
        title = "C_m as a Function of Velocity"
        gfs.gen_plotter(func_mat[:,2],data_mat[:,3],xlabel='Velocity',ylabel='Cm',title=title,save='fig_Cm.png')
        title = "L/D as a Function of Velocity"
        gfs.gen_plotter(func_mat[:,2],func_mat[:,1],xlabel='Velocity',ylabel='L/D',title=title,save='fig_L_D.png')
        title = "Drag as a Function of Velocity"
        gfs.gen_plotter(func_mat[:,2],func_mat[:,3],xlabel='Velocity',ylabel='Drag',title=title,save='fig_Drag.png')
        title = "Power Required as a Function of Velocity"
        gfs.gen_plotter(func_mat[:,2],func_mat[:,4],xlabel='Velocity',ylabel='Pwr_req',title=title,save='fig_Pwr.png')
        title = "Static Margin as a Function of Velocity"
        gfs.gen_plotter(func_mat[:,2],data_mat[:,6],xlabel='Velocity',ylabel='Static Margin',title=title,save='fig_SM.png')
        title = "Glide Ratio as a Function of Velocity"
        gfs.gen_plotter(func_mat[:,2],func_mat[:,5],xlabel='Velocity',ylabel='Glide Ratio',title=title,save='fig_GR.png')
        title = "Sink Rate as a Function of Velocity"
        gfs.gen_plotter(func_mat[:,2],func_mat[:,6],xlabel='Velocity',ylabel='Sink Rate',title=title,save='fig_SR.png')
        title = "C_m_a as a Function of Velocity"
        gfs.gen_plotter(func_mat[:,2],data_mat[:,5],xlabel='Velocity',ylabel='C_m_a',title=title,save='fig_Cma.png')
        title = "C_l_b as a Function of Velocity"
        gfs.gen_plotter(func_mat[:,2],data_mat[:,7],xlabel='Velocity',ylabel='C_l_b',title=title,save='fig_Clb.png')
        title = "C_n_b as a Function of Velocity"
        gfs.gen_plotter(func_mat[:,2],data_mat[:,8],xlabel='Velocity',ylabel='C_n_b',title=title,save='fig_Cnb.png')
        title = "Drag as a Function of Lift"
        gfs.gen_plotter(data_mat[:,1],data_mat[:,2],xlabel='LIFT',ylabel='Drag',title=title,save='fig_D(lift).png')

        title = "Locus of Aerodynamics Centers as a Function of Angle of Attack"
        set1 = "x_ac"
        set2 = "y_ac"
        gfs.gen_plotter(aero_mat[:,0],aero_mat[:,9], y2=aero_mat[:,10],xlabel='Angle of Attack',ylabel='AC Location',title=title,set1_label=set1,set2_label=set2,save='fig_ac_locus.png')

        # Create the x vs. y ac plot:
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.scatter(aero_mat[:,9],aero_mat[:,10])
        ax.set_xlabel("x_ac")
        ax.set_ylabel("y_ac")
        ax.set_title("Aerodynamic Centers; x vs y location")
        plt.show()
        fig.savefig("x_vs_y_ac.png")

    print(data_mat[:,1])
    print(data_mat[:,2])
    # print(data_mat)
    # print(func_mat)
    #print(aero_mat)



    ####---------------------------------------------------------------------####
    #                         Perform Polynomial Regression to determine D_coeffs
    lift = data_mat[:,1]
    drag = data_mat[:,2]
    D_coeffs = np.polyfit(lift,drag,2)
    #print(D_coeffs)
    C_D2 = D_coeffs[0]
    C_D1 = D_coeffs[1]
    C_D0 = D_coeffs[2]

    li = np.linspace(-0.25,1.5,30)
    dr = np.zeros((li.size))
    for i in range(li.size):
        dr[i] = C_D0 + (C_D1 * li[i]) + (C_D2 * (li[i]**2))
    
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot(li,dr,label="polyfit",color = "orange")
    ax.scatter(data_mat[:,1],data_mat[:,2],label="DATA")
    ax.legend()
    # if 'xlabel' in kwargs:
    #     ax.set_xlabel(kwargs['xlabel'])
    # if 'ylabel' in kwargs:
    #     ax.set_ylabel(kwargs['ylabel'])
    # if 'title' in kwargs:
    #     ax.set_title(kwargs['title'])
    plt.show()
    fig.savefig("Comparison.png") 
    # title = "Drag as a Function of Lift"
    #     gfs.gen_plotter(data_mat[:,1],data_mat[:,2],xlabel='LIFT',ylabel='Drag',title=title,save='fig_D(lift).png')

    ####---------------------------------------------------------------------####
    #                                               Perform Airspeed Calculations
    e = gfs.oef(C_D2,params["Ra"])
    V_min = gfs.V_min_func(params["CL_max"],params["W"],params["Sw"],params["rho"])
    V_MD = gfs.V_MD_func(e,params["Ra"],C_D0,params["W"],params["Sw"],params["rho"])
    V_MDV = gfs.V_MDV_func(e,params["Ra"],C_D1,C_D0,params["W"],params["Sw"],params["rho"])

    print("\n\nThe Stall speed is %.6f ft/s" %V_min)
    print("The Minimum Drag airspeed is %.6f ft/s" %V_MD)
    print("The Minimum Power airspeed is %.6f ft/s" %V_MDV)

    new_weight = gfs.kin_energy(params["g"],params["KE"], V= V_MD)
    print("\n\nThe new weight would be: ",new_weight)

    new_CL_desired = gfs.C_L_func(params["rho"],params["Sw"], W = params["W"], V = V_MD)
    print("\n\nThe new desired CL would be: ",new_CL_desired)


    ####---------------------------------------------------------------------####       
    #                                    Generate and Fill the Results Dictionary

    state = {
        "type" : "aerodynamic",
        "velocity" : params["vel"],
        "alpha" : 0.0,
        "beta" : 0.0
    }
    my_scene.set_aircraft_state(state = state, aircraft = "glider_2.0")
    # 0.440613
    target_CL = my_scene.target_CL(CL = 0.55, filename="target_CL.txt", set_state = True)
    FM_results = my_scene.solve_forces(filename="Output_FM",report_by_segment=True, dimensional=False, non_dimensional=True,)
    derivs = my_scene.derivatives(filename="Output_Derivs",report_by_segment=True)
    results = {}
    results["CY_p"] = derivs["glider_2.0"]["damping"]["Cy,pbar"]
    results["Cl_p"] = derivs["glider_2.0"]["damping"]["Cl,pbar"]
    results["Cn_p"] = derivs["glider_2.0"]["damping"]["Cn,pbar"]
    results["CL_q"] = derivs["glider_2.0"]["damping"]["CL,qbar"]
    results["CD_q"] = derivs["glider_2.0"]["damping"]["CD,qbar"]
    results["Cm_q"] = derivs["glider_2.0"]["damping"]["Cm,qbar"]
    results["CY_r"] = derivs["glider_2.0"]["damping"]["Cy,rbar"]
    results["Cl_r"] = derivs["glider_2.0"]["damping"]["Cl,rbar"]
    results["Cn_r"] = derivs["glider_2.0"]["damping"]["Cn,rbar"]
    
    # 0: aoa,  1: CL,  2: CD,  3: Cm,  4: CL_a,  5: Cm_a,  6: SM,  7: Cl_b,  8: Cn_b
    results["CL"] = FM_results["glider_2.0"]["total"]["CL"]
    results["CD"] = FM_results["glider_2.0"]["total"]["CD"]
    results["Cm"] = FM_results["glider_2.0"]["total"]["Cm"]
    results["CL_a"] = derivs["glider_2.0"]["stability"]["CL,a"]
    results["Cm_a"] = derivs["glider_2.0"]["stability"]["Cm,a"]
    results["SM"] = derivs["glider_2.0"]["stability"]["%_static_margin"]
    results["Cl_b"] = derivs["glider_2.0"]["stability"]["Cl,b"]
    results["Cn_b"] = derivs["glider_2.0"]["stability"]["Cn,b"]


    # 0: aoa,  1: L/D,  2: Vel,  3: Drag,  4: Pwr,  5: GR,  6: Sr
    results["L/D"] = results["CL"] / results["CD"] 
    results["V"] = gfs.C_L_func(params["rho"],params["Sw"],C_L = results["CL"], W= params["W"])
    results["Drag"] = gfs.drag_func(results["CD"],results["V"],params["rho"],params["Sw"])
    results["Pwr"] = results["V"] * results["Drag"]
    results["GR"] = results["CL"] / results["CD"] 
    results["Sr"] = results["Pwr"] / params["W"]


    results["C_D2"]= C_D2
    results["C_D1"]= C_D1
    results["C_D0"]= C_D0
    results["e"]= e
    results["V_stall"]= V_min
    results["V_min_drag"]= V_MD
    results["V_min_pwr"]= V_MDV
    results["New_weight"]= new_weight
    results["New_CL_Desired"]= new_CL_desired


    ####---------------------------------------------------------------------####
    #                                             Perform Additional Calculations
    # Determine the MAX allowable arm:
    max_arm = 6 + (params["x_vstab"])
    params["max_arm"] = max_arm

    area_w = gfs.area2412(params["wing_avg_c"])
    area_h = gfs.area0012(params["hstab_c"])
    area_v = gfs.area0012(params["vstab_c"])
    wing_volume = (area_w * params["b"]) + (area_h * params["hstab_b"]) + (area_v * params["vstab_b"])
    W_w = wing_volume * params["rho_f"]
    W_d = gfs.W_d_func(params["rho_d"],params["ddiam"], dowel = params["dowel_allow"])
    W_a = W_w + W_d
    current_b_mass = params["W"] - W_a
    current_b_thick = gfs.W_b_func(params["rho_b"],params["bdiam_inner"],params["bdiam_outer"],W_b = current_b_mass)
    
    params["W_w"] = W_w
    params["W_d"] = W_d
    params["W_a"] = W_a
    params["current_b_mass"] = current_b_mass
    params["current_b_thick"] = current_b_thick


    target_arm = bbt.balbal(params, results,False, b_mass = current_b_mass)
    params["target_arm"] = target_arm

    # If we were to max the moment are (which simultaneously reduces the weight):
    reduced_mass = bbt.balbal(params,results,False, b_arm = max_arm, mass = True)
    params["reduced_mass"] = reduced_mass
   

    print("\nThe current_b_mass is: ", current_b_mass)
    print("The target arm, therefore, is: ", target_arm)
    print("\nIf we were to max the mmt arm at: ", max_arm)
    print("Then the required reduced mass would be: ",reduced_mass)



    ####---------------------------------------------------------------------####
    #                                                 Prepare the simulator .JSON

    l_wt = 1.1 * ((params["wing_avg_c"]*-0.75) - params["x_hstab"])
    first = (4 * params["Sw_H"] * l_wt) / (np.pi * (params["b"]**2) * params["wing_avg_c"])
    CL_ahat = first * params["CLa_2412"] * params["CLa_0012"]
    Cm_ahat = CL_ahat * params["x_hstab"]

    with open("0000.json", "r") as f:
        sim_json = json.load(f)

        sim_json["aircraft"]["name"] = "OUR_Glider"
        sim_json["aircraft"]["wing_area[ft^2]"] = params["Sw"]
        sim_json["aircraft"]["wing_span[ft]"] = params["b"]
        sim_json["operating"]["weight[lbf]"] = params["W"]
        sim_json["operating"]["climb[deg]"] = 0.0
        sim_json["operating"]["density[slugs/ft^3]"] = 0.0023769
        
        sim_json["reference"]["CL"] = results["CL"] # data_mat[15][1]

        sim_json["reference"]["Ixx[slugs*ft^2]"] = 676.9938/144
        sim_json["reference"]["Iyy[slugs*ft^2]"] = 400.91910371/144
        sim_json["reference"]["Izz[slugs*ft^2]"] = 1076.76979/144
        sim_json["reference"]["Ixy[slugs*ft^2]"] = -0.00000252/144
        sim_json["reference"]["Ixz[slugs*ft^2]"] = 0.47384/144
        sim_json["reference"]["Iyz[slugs*ft^2]"] = 0.00001588/144
        sim_json["reference"]["CD0"] = C_D0
        sim_json["reference"]["CD1"] = C_D1
        sim_json["reference"]["CD2"] = C_D2

        sim_json["reference"]["CL,a"] = derivs["glider_2.0"]["stability"]["CL,a"]

        sim_json["reference"]["Cm,a"] = derivs["glider_2.0"]["stability"]["Cm,a"]

        sim_json["reference"]["CL,ahat"] = CL_ahat
        sim_json["reference"]["Cm,ahat"] = Cm_ahat

        sim_json["reference"]["CY,b"] = derivs["glider_2.0"]["stability"]["Cy,b"]
        sim_json["reference"]["Cl,b"] = derivs["glider_2.0"]["stability"]["Cl,b"]
        sim_json["reference"]["Cm,b"] = derivs["glider_2.0"]["stability"]["Cm,b"]
        sim_json["reference"]["Cn,b"] = derivs["glider_2.0"]["stability"]["Cn,b"]

        sim_json["reference"]["CY,p"] = results["CY_p"]
        sim_json["reference"]["Cl,p"] = results["Cl_p"]
        sim_json["reference"]["Cn,p"] = results["Cn_p"]
        sim_json["reference"]["CL,q"] = results["CL_q"]
        sim_json["reference"]["CD,q"] = results["CD_q"]
        sim_json["reference"]["Cm,q"] = results["Cm_q"]
        sim_json["reference"]["CY,r"] = results["CY_r"]
        sim_json["reference"]["Cl,r"] = results["Cl_r"]
        sim_json["reference"]["Cn,r"] = results["Cn_r"]


    with open("edited_0000.json", "w") as edited:
        json.dump(sim_json, edited,indent = 4)



    return results, params


input_file = "scene.json"
(results,params) = main(input_file,False)

print("\n\nThe PARAMETERS Dictionary is as follows: ")
pprint.pprint(params)
print("\n\nThe RESULTS Dictionary is as follows: ")
pprint.pprint(results)

launch_v = gfs.kin_energy(params["g"],params["KE"],W = params["W"])
print(launch_v)
