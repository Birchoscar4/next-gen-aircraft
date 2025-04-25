#Weight & Balance Analysis Function 
#Jude Doherty 25/04/2025
import pandas as pd
import math
import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
pi = float(math.pi)
plt.ion()

#############################################################
#                                                           #
#               WEIGHT & BALANCE ANALYSIS                   #
#                                                           #
#############################################################

# Load the component masses CSV from structures output
def WeightBalance():

    W_LH2, W_tank, r, L = LH2_tank_sizing_code()
    tank_location_x = input("Input X coordinate of LH2 tank location")
    tank_location_y
    tank_location_z
    
    ld3_location = input("Input X coordinate of LD3 containers location")

    
    masses_df = pd.read_csv('aircraft_masses.csv')
    

    # Check columns for spaces
    masses_df.columns = [col.strip() for col in masses_df.columns]

    # Now set 'Component' as the index
    masses_df.set_index('Component', inplace=True)

    #Weight and Balance Array
    wb_data = {
        'Component': ['Fuselage', 'Aft Fuselage', 'Wings', 'Vertical Tails', 'LH2 Tanks', 'Engines', 'Passengers and Crew', 'LD3 Containers', 'Electrical Systems & Actuators', 'MLG', 'NLG'],
        'Weight (kg)': [ float(masses_df.loc['Fuselage', 'mass']), float(masses_df.loc['Aft Body', 'mass']), 2*float(masses_df.loc['Wing', 'mass']), 2*float(masses_df.loc['vstab', 'mass']), (W_LH2 + W_Tank), 3902, 12580, 6352, 5836, 2*float(masses_df.loc['Main Landing Gear', 'mass']), float(masses_df.loc['Nose Landing Gear', 'mass'])],
        'Xi (m)': [float(masses_df.loc['Fuselage', 'x_cg']), float(masses_df.loc['Aft Body', 'x_cg']), float(masses_df.loc['Wing', 'x_cg']), float(masses_df.loc['vstab', 'x_cg']), 6.2295, 15.775, 8.43, 4.2265, 8.23, float(masses_df.loc['Main Landing Gear', 'x_cg']), float(masses_df.loc['Nose Landing Gear', 'x_cg'])], #x distance from RP
        'Yi (m)': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], #z distance from RP
        'Zi (m)': [float(masses_df.loc['Fuselage', 'z_cg']), float(masses_df.loc['Aft Body', 'z_cg']), float(masses_df.loc['Wing', 'z_cg']), float(masses_df.loc['vstab', 'z_cg']), -1.118, 2, 0.4, -0.955, -1, float(masses_df.loc['Main Landing Gear', 'z_cg']), float(masses_df.loc['Nose Landing Gear', 'z_cg'])], #y distance from RP
        'Ixx Personal': [float(masses_df.loc['Fuselage', 'Ixx']), float(masses_df.loc['Aft Body', 'Ixx']), 2*float(masses_df.loc['Wing', 'Ixx']), 2*float(masses_df.loc['vstab', 'Ixx']), 2*(0.5*((W_LH2+W_Tank)/2)*(r**2)), 36072, 0, 0, 0, 2*float(masses_df.loc['Main Landing Gear', 'Ixx']), float(masses_df.loc['Nose Landing Gear', 'Ixx']),],
        'Iyy Personal': [float(masses_df.loc['Fuselage', 'Iyy']), float(masses_df.loc['Aft Body', 'Iyy']), 2*float(masses_df.loc['Wing', 'Iyy']), 2*float(masses_df.loc['vstab', 'Iyy']), 2*((1/12)*((W_LH2+W_Tank/2)*((3*(r**2)+(l**2)))), 19735, 0, 3377, 0, 2*float(masses_df.loc['Main Landing Gear', 'Iyy']), float(masses_df.loc['Nose Landing Gear', 'Iyy']),], 
        'Izz Personal': [float(masses_df.loc['Fuselage', 'Izz']), float(masses_df.loc['Aft Body', 'Izz']), 2*float(masses_df.loc['Wing', 'Izz']), 2*float(masses_df.loc['vstab', 'Iyy']), 2*((1/12)*((W_LH2+W_Tank/2)*((3*(r**2)+(l**2)))), 19735, 0, 2490, 0, 2*float(masses_df.loc['Main Landing Gear', 'Izz']), float(masses_df.loc['Nose Landing Gear', 'Izz']), ]    
    }
    print(wb_data)
    df = pd.DataFrame(wb_data)

    # Calculate moments
    df['X Moment'] = df['Weight (kg)'] * df['Xi (m)'] 
    df['Y Moment'] = df['Weight (kg)'] * df['Yi (m)']
    df['Z Moment'] = df['Weight (kg)'] * df['Zi (m)']


    # Calculate total weight and total moments
    total_weight = df['Weight (kg)'].sum()
    total_moment_x = df['X Moment'].sum()
    total_moment_y = df['Y Moment'].sum()
    total_moment_z = df['Z Moment'].sum()

    # Calculate CG positions
    cg_x = total_moment_x / total_weight
    cg_y = total_moment_y / total_weight
    cg_z = total_moment_z / total_weight

    print ('X CoG: ', (cg_x))
    print ('Y CoG: ', (cg_y))
    print ('Z CoG: ', (cg_z))
    

    # Calculate moments of inertia 
    #Global Intertia Calculations
    df['Ixx'] = df['Weight (kg)'] * (((df['Yi (m)'] - cg_y)**2) + ((df['Zi (m)'] - cg_z)**2)); 
    df['Iyy'] = df['Weight (kg)'] * (((df['Zi (m)'] - cg_z)**2) + ((df['Xi (m)'] - cg_x)**2)); 
    df['Izz'] = df['Weight (kg)'] * (((df['Xi (m)'] - cg_x)**2) + ((df['Yi (m)'] - cg_y)**2)); 
    df['Ixy'] = df['Weight (kg)'] * ((df['Xi (m)'] - cg_x) * (df['Yi (m)'] - cg_y)); 
    df['Iyz'] = df['Weight (kg)'] * ((df['Yi (m)'] - cg_y) * (df['Zi (m)'] - cg_z)); 
    df['Izx'] = df['Weight (kg)'] * ((df['Zi (m)'] - cg_z) * (df['Xi (m)'] - cg_x)); 

    Ixx_total = df['Ixx'].sum() + df['Ixx Personal'].sum(); 
    Iyy_total = df['Iyy'].sum() + df['Iyy Personal'].sum(); 
    Izz_total = df['Izz'].sum() + df['Izz Personal'].sum(); 
    Ixy_total = df['Ixy'].sum(); 
    Iyz_total = df['Iyz'].sum(); 
    Izx_total = df['Izx'].sum(); 

    print ('Ixx total: ', (Ixx_total))
    print ('Iyy total: ', (Iyy_total))
    print ('Izz total: ', (Izz_total))

    return masses_df, total_weight, Ixx_total, Iyy_total, Izz_total, Izx_total

#############################################################
#                                                           #
#               STATIC STABILITY ANALYSIS                   #
#                                                           #
#############################################################
def StaticStability():

 
    # Path to your .stab file
    stab_file_path = 'T_area_12_angle_60.stab'

    # Empty dictionary to store the values
    flight_data = {}

    # Read file
    with open(stab_file_path, 'r') as file:
        for line in file:
            if line.strip() and not line.strip().startswith("#"):
                parts = line.split()
                if len(parts) >= 2:
                    name = parts[0]
                    value = float(parts[1])
                    unit = parts[2] if len(parts) > 2 else ""
                    
                    # Store into dictionary
                    flight_data[name] = {'Value': value, 'Unit': unit}
                    
                    # Check if we reached Yaw___Rate
                    if name == "Yaw___Rate":
                        break  # ⛔ Stop reading further!

    # Convert to a DataFrame for easy viewing
    flight_data_df = pd.DataFrame.from_dict(flight_data, orient='index')

    mach_number = flight_data_df.loc['Mach_', 'Value']; 
    velocity = flight_data_df.loc['Vinf_', 'Value'];  
    rho = flight_data_df.loc['Rho_', 'Value']; 
    q_dash0 = 0.5*rho*(velocity**2); 
    cd_0 = 0.021; 
    wing_area = flight_data_df.loc['Sref_', 'Value']; 
    thrust_mach_number = -12590; 
    mac = flight_data_df.loc['Cref_', 'Value']; 
    wingspan = flight_data_df.loc['Bref_', 'Value']; 
    gravity = 9.81

    # Function to read OpenVSPAero .stab file
    def read_stab_file(filename):
        with open(filename, 'r') as file:
            lines = file.readlines()
        
        # Find the header row that contains 'Coef' (marks the start of the table)
        for i, line in enumerate(lines):
            if "Coef" in line:
                header_row = i
                break

        # Extract column names
        columns = lines[header_row].split()
        
        # Extract data rows (skipping header and comment lines)
        data = [line.split() for line in lines[header_row + 3:] if not line.startswith("#")]

        # Separate stability derivatives from SM and x_np
        full_data = []
        extra_data = {}

        for row in data:
            if len(row) == len(columns):  # Full row with all derivatives
                full_data.append(row)
            elif len(row) == 2:  # SM and x_np (only two values)
                extra_data[row[0]] = float(row[1])

        # Convert full stability derivative table to a DataFrame
        df_stab = pd.DataFrame(full_data, columns=columns)
        
        # Convert numeric columns to float
        df_stab.iloc[:, 1:] = df_stab.iloc[:, 1:].astype(float)

        # Convert SM and x_np to a separate DataFrame
        df_extra = pd.DataFrame(extra_data.items(), columns=["Coef", "Value"])

        return df_stab, df_extra

    def extract_sm_xnp(filename):
        with open(filename, 'r') as file:
            lines = file.readlines()

        extra_data = {}
        found_section = False

        for line in lines:
            # Check if we've reached the extra data section
            if "Result" in line and "Value" in line:  
                found_section = True
                continue  # Skip to next line

            if found_section:
                line = line.strip()
                if line == "" or line.startswith("#"):  # Skip empty/comment lines
                    continue
                
                values = line.split()
                if len(values) >= 2:  # Ensure at least "SM <value>"
                    key = values[0]  # First item is the coefficient name
                    try:
                        value = float(values[1])  # Convert second item to float
                        extra_data[key] = value
                    except ValueError:
                        print(f"Warning: Could not convert {values} to float.")

        # Convert extracted data into a DataFrame
        df_extra = pd.DataFrame(extra_data.items(), columns=["Coef", "Value"])

        if df_extra.empty:
            print("No SM or X_np values were found. Check the file format!")

        return df_extra

    #Extracting Cm Alpha Dot
    def extract_cm_q_alpha_dot(filepath):
        with open(filepath, 'r') as file:
            for line in file:
                if "CMm_(q + alpha_dot)" in line:
                    # Split the line and take the last value
                    value = float(line.strip().split()[-1])
                    return value
        return None  # if not found

    #Load .qstab file
    qstab_file = 'Poster_Showcase_Model 1_DegenGeom_Q.qstab'  # Replace with your actual file path
    cm_q_plus_alpha_dot = extract_cm_q_alpha_dot(qstab_file)

    # Load .stab file
    stab_data, extra_data = read_stab_file(stab_file_path)
    extra_data = extract_sm_xnp(stab_file_path)

    # Display DataFrames
    print("🔹 Stability Derivative Data:")
    print(stab_data)

    print("\n🔹 Additional Stability Data (SM & x_np):")
    print(extra_data)

    #Forward Speed Disturbance Calculations 
    cd_u = stab_data.loc[stab_data['Coef'] == 'CD', 'U'].values[0]; 
    Ct_xu = ((mach_number/(q_dash0*wing_area))*thrust_mach_number)-(2*cd_0); 
    print("Forward Speed Disturbance Coefficient:")
    print (Ct_xu - cd_u); 
    if (Ct_xu - cd_u) < 0: 
        print ("Aircraft Derivative is Stable")
    else: 
        print ("Aircraft Derivative is not Stable, optimise configuration")

    print("")

    #Side Speed Disturbance Calculations
    cy_beta = stab_data.loc[stab_data['Coef'] == 'CFy', 'Beta'].values[0]; 
    print("Side Speed Disturbance Coefficient:")
    print (cy_beta) 
    if cy_beta < 0: 
        print ("Aircraft Derivative is Stable")
    else: 
        print ("Aircraft Derivative is not Stable, optimise configuration")

    print("")

    #Vertical Speed Disturbance Calculations 
    cl_alpha = stab_data.loc[stab_data['Coef'] == 'CL', 'Alpha'].values[0]; 
    print("Vertical Speed Disturbance Coefficient:")
    print (cl_alpha) 
    if cl_alpha > 0: 
        print ("Aircraft Derivative is Stable")
    else: 
        print ("Aircraft Derivative is not Stable, optimise configuration")

    print("")

    #Dihedral Effect Disturbance Calculations
    c_i_beta = stab_data.loc[stab_data['Coef'] == 'CMl', 'Beta'].values[0]; 
    print("Dihedral Effect Disturbance Coefficient:")
    print (c_i_beta) 
    if c_i_beta < 0: 
        print ("Aircraft Derivative is Stable")
    else: 
        print ("Aircraft Derivative is not Stable, optimise configuration")

    print("")

    #Roll Rate Disturbance Calculations 
    c_i_p = stab_data.loc[stab_data['Coef'] == 'CMl', 'p'].values[0]; 
    print("Roll Rate Disturbance Coefficient:")
    print (c_i_p) 
    if c_i_p < 0: 
        print ("Aircraft Derivative is Stable")
    else: 
        print ("Aircraft Derivative is not Stable, optimise configuration")

    print("")

    #Static Longitudinal Stability Calculations 
    cm_alpha = stab_data.loc[stab_data['Coef'] == 'CMm', 'Alpha'].values[0]; 
    print("Static Longitudinal Stability Coefficient:")
    print (cm_alpha) 
    if cm_alpha < 0: 
        print ("Aircraft Derivative is Stable")
    else: 
        print ("Aircraft Derivative is not Stable, optimise configuration")

    print("")

    #Pitch Rate Disturbance Calculations
    cm_q = stab_data.loc[stab_data['Coef'] == 'CMm', 'q'].values[0]; 
    print("Pitch Rate Disturbance Coefficient:")
    print (cm_q) 
    if cm_q < 0: 
        print ("Aircraft Derivative is Stable")
    else: 
        print ("Aircraft Derivative is not Stable, optimise configuration")

    print("")

    #Static Directional Stability Calculations 
    cn_beta = stab_data.loc[stab_data['Coef'] == 'CMn', 'Beta'].values[0]; 
    print("Static Directional Stability Coefficient:")
    print (cn_beta) 
    if cn_beta > 0:
        print ("Aircraft Derivative is Stable ")
    else: 
        print ("Aircraft Derivative is not Stable, optimise configuration")

    print("")

    #Yaw Rate Disturbance Calculations 
    cn_r = stab_data.loc[stab_data['Coef'] == 'CMn', 'r'].values[0]; 
    print("Yaw Rate Disturbance Coefficient:")
    print (cn_r) 
    if cn_r < 0: 
        print ("Aircraft Derivative is Stable")
    else: 
        print ("Aircraft Derivative is not Stable, optimise configuration")

    print("")

    #Static Margin 
    static_margin = extra_data.loc[extra_data['Coef'] == 'SM', 'Value'].values[0]; 
    print("Static Margin:")
    print (static_margin) 
    if static_margin > 0: 
        print ("Aircraft Derivative is Stable")
    else: 
        print ("Aircraft Derivative is not Stable, optimise configuration")

    print("")

    return velocity, gravity, q_dash0, wing_area, wingspan, mac, rho, stab_data, cd_0, cn_beta, c_i_beta, cn_r, c_i_p, cm_q_plus_alpha_dot, cy_beta, Ct_xu, cm_q, cl_alpha, cm_alpha
 
#############################################################
#                                                           #
#               DYNAMIC STABILITY ANALYSIS                  #
#                                                           #
#############################################################
def DynamicStabilityControl():

    # Obtaining values from Weight & Balance Function
    masses_df, total_weight, Ixx_total, Iyy_total, Izz_total, Izx_total = WeightBalance()
    
    # Obtaining Values from Static Stability function
    (velocity, gravity, q_dash0, wing_area, wingspan, mac, rho, 
     stab_data, cd_0, cn_beta, c_i_beta, cn_r, c_i_p, 
     cm_q_plus_alpha_dot, cy_beta, Ct_xu, cm_q, cl_alpha, cm_alpha) = StaticStability()

    ########## Required Stability Derivatives##############

    cl_r = stab_data.loc[stab_data['Coef'] == 'CMl', 'r'].values[0]; 
    cy_r = stab_data.loc[stab_data['Coef'] == 'CFy', 'r'].values[0]; 
    cl_u = stab_data.loc[stab_data['Coef'] == 'CL', 'U'].values[0]; 
    cl_0 = (total_weight*gravity)/(q_dash0*wing_area)

    #Control Surface Coeff. 
    ae_Z_delta_e = stab_data.loc[stab_data['Coef'] == 'CL', 'ConGrp_1'].values[0]; #Outboard Elevon
    ie_Z_delta_e = stab_data.loc[stab_data['Coef'] == 'CL', 'ConGrp_2'].values[0]; #Inboard Elevon

    ae_M_delta_e = stab_data.loc[stab_data['Coef'] == 'CMm', 'ConGrp_1'].values[0]; #Outboard Elevon
    ie_M_delta_e = stab_data.loc[stab_data['Coef'] == 'CMm', 'ConGrp_2'].values[0]; #Inboard Elevon

    ae_L_delta_a = stab_data.loc[stab_data['Coef'] == 'CMl', 'ConGrp_1'].values[0]; #Outboard Elevon

    C_n_delta_r = stab_data.loc[stab_data['Coef'] == 'CMn', 'ConGrp_3'].values[0]; #Rudder

    #Elevator, Aileron and Rudder Power Control Derivative Calculations 
    Z_deltaE = (-q_dash0*wing_area*(ae_Z_delta_e + ie_Z_delta_e))/total_weight
    M_deltaE = (q_dash0*wing_area*mac*(ae_M_delta_e + ie_M_delta_e))/(Iyy_total)
    L_delta_a = (q_dash0*wing_area*wingspan*ae_L_delta_a)/Ixx_total
    N_delta_r = (q_dash0*wing_area*wingspan*C_n_delta_r)/Izz_total

    #################DYNAMIC MODE CALCS######################################

    #Rolling Mode Approximation Calculations 
    l_p = (q_dash0*wing_area*wingspan*wingspan*c_i_p)/(2*Ixx_total*velocity)
    t_roll = -1/l_p; 

    print("Roll Time Period (s)")
    print(t_roll)

    if t_roll < 1.4:
        print("Aircraft Meets Level 1 Flying Requirements for Rolling Mode")
    elif t_roll < 3:
        print("Aircraft Meets Level 2 Flying Requirements for Rolling Mode")
    elif t_roll < 10:
        print("Aircraft Meets Level 3 Flying Requirements for Rolling Mode")
    else: 
        print("Aircraft is unacceptable for Rolling Mode, optimise configuration")

    print("")
    ######Roll Mode Plotting########

    # Numerator and denominator coefficients of the transfer function
    s_roll, t_plot_roll = sp.symbols('s t'); 

    roll_step = (-10*pi)/180 #5 Degree Input

    phi_roll_s = ((roll_step*L_delta_a)/((s_roll * (s_roll - l_p))))
    phi_roll_t = sp.inverse_laplace_transform(phi_roll_s, s_roll, t_plot_roll)

    #Plot Function
    phi_func_roll = sp.lambdify(t_plot_roll, phi_roll_t, modules=['numpy'])
    t_values_roll = np.linspace(0, 1, 500)
    phi_values_roll = phi_func_roll(t_values_roll)

    plt.figure(figsize=(10, 6))
    plt.plot(t_values_roll, np.rad2deg(phi_values_roll), label='Roll Angle (deg)', color='r')  # Convert to degrees
    plt.xlabel('Time (s)')
    plt.ylabel('Roll Angle, phi(t) (deg)')
    plt.title('Roll Response to 10° Aileron Step Input')
    plt.grid()
    plt.legend()
    plt.show()


    #Spiral Mode Approximation Calculations 
    l_beta = (q_dash0*wing_area*wingspan*c_i_beta)/(Ixx_total)
    n_r = (q_dash0*wing_area*wingspan*wingspan*cn_r)/(2*Izz_total*velocity)
    l_r = (q_dash0*wing_area*wingspan*wingspan*cl_r)/(2*Ixx_total*velocity)
    n_beta = (q_dash0*wing_area*wingspan*cn_beta)/(Izz_total)

    spiral_time_constant = (l_beta*0.6931471806)/((n_beta*l_r)-(l_beta*n_r))

    print("Spiral to Double Amplitude Constant (s)")
    print(spiral_time_constant)

    print("")

    if spiral_time_constant > 20:
        print("Aircraft Meets Level 1 Flying Requirements for Spiral Mode")
    elif spiral_time_constant > 12:
        print("Aircraft Meets Level 2 Flying Requirements for Spiral Mode")
    elif spiral_time_constant > 4:
        print("Aircraft Meets Level 3 Flying Requirements for Spiral Mode")
    else: 
        print("Aircraft is unacceptable for Spiral Mode, optimise configuration")

    print("")

    #Dutch-Roll Mode Approximation Calculations 
    y_beta = (q_dash0*wing_area*cy_beta)/total_weight
    y_r = (q_dash0*wing_area*wingspan*cy_r)/(2*total_weight*velocity)

    omega_ndt = math.sqrt((1/velocity)*((y_beta*n_r)+(n_beta*velocity)-(n_beta*y_r)))
    zeta_dt = (-1/(2*omega_ndt))*(n_r+(y_beta/velocity))

    print("Dutch Roll Zeta")
    print(zeta_dt)

    print("")

    print("Dutch Roll Natural Frequency")
    print(omega_ndt)

    print("")

    dt_real_part = -omega_ndt*zeta_dt
    dt_imaginary_part = omega_ndt*math.sqrt(abs((zeta_dt**2)-1))
    print(f"The Dutch Roll Transfer Function solutions are: {dt_real_part} + {dt_imaginary_part}i and {dt_real_part} - {dt_imaginary_part}i")

    print("")
    if zeta_dt > 0.08 and (-dt_real_part) > 0.35 and omega_ndt > 0.4:
        print("Aircraft Meets Level 1 Flying Requirements for Dutch Roll Mode")
    elif zeta_dt > 0.02 and (-dt_real_part) > 0.05 and omega_ndt > 0.4:
        print("Aircraft Meets Level 2 Flying Requirements for Dutch Roll Mode")
    elif zeta_dt > 0.02 and omega_ndt > 0.4: 
        print("Aircraft Meets Level 3 Flying Requirements for Dutch Roll Mode")
    else: 
        print("Aircraft is unacceptable for Dutch Roll, optimise configuration")

    print("")
    #Phugoid Approximation Calculations

    z_u = (-q_dash0*wing_area*(cl_u+(2*cl_0)))/(total_weight*velocity)
    x_u = (-q_dash0*wing_area*(2*cd_0))/(total_weight*velocity)

    omega_np = math.sqrt((-z_u*9.81)/velocity); 
    zeta_p = -x_u/(2*omega_np)

    p_real_part = -omega_np*zeta_p
    p_imaginary_part = omega_np*math.sqrt(abs((zeta_p**2)-1))

    print("")

    print(f"The Phugoid Transfer Function solutions are: {p_real_part} + {p_imaginary_part}i and {p_real_part} - {p_imaginary_part}i") 

    print("")

    print("Phugoid Damping Ratio")
    print(zeta_p)

    print("")

    if zeta_p > 0.04:
        print("Aircraft Meets Level 1 Flying Requirements for Phugoid")
    elif zeta_p > 0:
        print("Aircraft Meets Level 2 Flying Requirements for Phugoid")
    else: 
        print("Aircraft has Level 3 or worse Flying Requirements for Phugoid") 

    print("")
    ########Phugoid Plotting###########
    x_tu = ((q_dash0*wing_area)*(Ct_xu + (2*cd_0)))/(total_weight*velocity)

    # Numerator and denominator coefficients of the transfer function
    s_phugoid, t_phugoid = sp.symbols('s t'); 

    elevator_step_phugoid = (15*pi)/180 #5 Degree Input

    #Pitch Transfer Function
    theta_phugoid_s = (elevator_step_phugoid * Z_deltaE * (s_phugoid - x_u - x_tu))/( (-velocity * ((s_phugoid**2) - ((x_u*x_tu)*s_phugoid) - ((z_u*gravity)/velocity)))) 

    theta_phugoid_t = sp.inverse_laplace_transform(theta_phugoid_s, s_phugoid, t_phugoid)

    #Plot Function
    theta_pugoid_func = sp.lambdify(t_phugoid, theta_phugoid_t, modules=['numpy'])
    t_values_theta_phugoid = np.linspace(0, 360, 500)
    theta_values_phugoid = theta_pugoid_func(t_values_theta_phugoid)

    plt.figure(figsize=(10, 6))
    plt.plot(t_values_theta_phugoid, np.rad2deg(theta_values_phugoid), label='Pitch Attitude (deg)', color='r')  # Convert to degrees
    plt.xlabel('Time (s)')
    plt.ylabel('Pitch Attitude Angle, theta(t) (deg)')
    plt.title('Phugoid Response to 30° Elevon Step Input')
    plt.grid()
    plt.legend()
    plt.show()


    #Short Period Approximation Calculations 
    cm_alpha_dot = cm_q_plus_alpha_dot - cm_q; 

    Z_alpha = (-q_dash0*wing_area*(cl_alpha+cd_0))/total_weight
    m_q = (q_dash0*wing_area*(mac**2)*cm_q)/(2*Iyy_total*velocity)
    m_alpha = (q_dash0*wing_area*mac*cm_alpha)/(Iyy_total)
    m_alpha_dot = (q_dash0*wing_area*(mac**2)*cm_alpha_dot)/(2*Iyy_total*velocity)

    omega_nsp = math.sqrt(((Z_alpha*m_q)/velocity)-m_alpha); 
    zeta_sp = (-(m_q+(Z_alpha/velocity)+m_alpha_dot))/(2*omega_nsp)


    sp_real_part = -omega_nsp*zeta_sp
    sp_imaginary_part = omega_nsp*(math.sqrt(abs((zeta_p**2)-1)))

    print(f"The solutions are: {sp_real_part} + {sp_imaginary_part}i and {sp_real_part} - {sp_imaginary_part}i")

    print("")
    print("Short Period Damping Ratio")
    print(zeta_sp)

    print("")
    print("Short Period Natural Freq.")
    print(omega_nsp)


    print("")
    if zeta_sp > 0.3 and zeta_sp < 2.0:
        print("Aircraft Meets Level 1 Flying Requirements for Short Period")
    elif zeta_sp > 0.2 and zeta_sp < 2.0:
        print("Aircraft Meets Level 2 Flying Requirements for Short Period")
    elif zeta_sp > 0.15: 
        print("Aircraft Meets Level 3 Flying Requirements for Short Period")
    else: 
        print("Aircraft is unacceptable for Short Period, optimise configuration")

    print("")

    ######Short Period Plotting########
    # Numerator and denominator coefficients of the transfer function
    s_sp, t_sp = sp.symbols('s t'); 

    elevator_step_sp = (10*pi)/180 #5 Degree Input

    #Pitch Transfer Function
    alpha_sp_s = (elevator_step_sp * ((Z_deltaE * s_sp) + ((M_deltaE*velocity) - (m_q*Z_deltaE))))/(velocity * ((s_sp**2) - (s_sp * (m_q + (Z_alpha/velocity) + m_alpha_dot)) + (((Z_alpha*m_q)/velocity) - m_alpha)))

    alpha_sp_t = sp.inverse_laplace_transform(alpha_sp_s, s_sp, t_sp)

    #Plot Function
    alpha_sp_func = sp.lambdify(t_sp, alpha_sp_t, modules=['numpy'])
    t_values_alpha_sp = np.linspace(0, 8, 500)
    alpha_values_sp = alpha_sp_func(t_values_alpha_sp)

    plt.figure(figsize=(10, 6))
    plt.plot(t_values_alpha_sp, np.rad2deg(alpha_values_sp), label='AoA (deg)', color='r')  # Convert to degrees
    plt.xlabel('Time (s)')
    plt.ylabel('AoA, alpha(t) (deg)')
    plt.title('Short Period Response to 30° Elevon Step Input')
    plt.grid()
    plt.legend()
    plt.show()

    theta_sp_s = (elevator_step_sp * ((s_sp * ((velocity * M_deltaE) + (Z_deltaE * m_alpha_dot))))  + ((m_alpha * Z_deltaE) - (Z_alpha * M_deltaE))) / (s_sp * velocity * ((s_sp**2) - (s_sp * (m_q + (Z_alpha/velocity) + m_alpha_dot)) + (((Z_alpha*m_q)/velocity) - m_alpha)))
    theta_sp_t = sp.inverse_laplace_transform(theta_sp_s, s_sp, t_sp)

    theta_sp_func = sp.lambdify(t_sp, theta_sp_t, modules=['numpy'])
    t_values_theta_sp = np.linspace(0, 8, 500)
    theta_values_sp = theta_sp_func(t_values_theta_sp)

    plt.figure(figsize=(10, 6))
    plt.plot(t_values_theta_sp, np.rad2deg(theta_values_sp), label='Theta (deg)', color='r')  # Convert to degrees
    plt.xlabel('Time (s)')
    plt.ylabel('Pitch Attitude, theta(t) (deg)')
    plt.title('Short Period Response to 30° Elevon Step Input')
    plt.grid()
    plt.legend()
    plt.show()



    ###############Controlability Calculations#############
    #Directional Control in event of Engine Out Yawing Moment
    total_thrust = 51857.087041393

    max_rudder_deflection = (15*pi)/180 #Assumed to be maximum Rudder deflection
    v_mc = math.sqrt(((total_thrust/2)*3.5)/(0.5*rho*wing_area*wingspan*C_n_delta_r))

    print("Minimum Control Speed at Altitude"); 
    print(v_mc)

    print("")

    exit_code = input("Press Enter to Exit")


