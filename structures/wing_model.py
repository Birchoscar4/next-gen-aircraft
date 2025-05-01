def compute_wing_mass(spar_cap_radius, num_stringers, skin_thickness, spar_thickness, max_q_velocity, max_q_density, load_factors, LE_spar, TE_spar, file_path):
    import numpy as np
    import re
    import pandas as pd

    buckled = 0
    load_factor = load_factors

    #file_path = "pressure.slc"
    aerofoil_file = "C:\\0714sc.dat"
    rho = max_q_density #0.36518 # kg/m3 - air
    v = max_q_velocity #210 # m/s
    q_inf = 0.5 * rho * v **2
    #LE_spar = 0.1
    #TE_spar = 0.7
    wingbox_fraction = TE_spar - LE_spar
    safetyfactor = 1.5

    #Material Properties
    E = 6.8947573e10 # Pa
    v = 0.4
    G = E / (2*(1+v))
    density = 1550 # kg/m3

    #Structural Parameters
    I_spar_cap=np.pi*(2*spar_cap_radius)**4/64
    A_spar_cap=np.pi*(spar_cap_radius**2)
    stringer_thickness = 0.006
    stringer_depth = 0.07
    stringer_width = 0.02
    A_stringer = stringer_thickness*stringer_depth + 2*(stringer_thickness*stringer_width)
    I_stringer = stringer_depth**3*stringer_thickness/12 + 2*(stringer_thickness**3*stringer_width/12 + stringer_thickness*stringer_width*(stringer_thickness+stringer_depth)**2/4)

    # Aerodynamics functions
    def compute_lift_forces(file_path, q_inf):
        slices = {}
        current_slice = None

        with open(file_path, "r", errors="ignore") as file:
            for line in file:
                # block headers
                match = re.match(r"BLOCK Cut_\d+_at_Y:_([\d\.\-]+)", line)
                if match:
                    current_slice = float(match.group(1))
                    slices[current_slice] = []
                    continue

                if current_slice is None:
                    continue

                if "Mach" in line or "Alpha" in line or "Beta" in line:
                    continue

                if "x" in line and "y" in line and "z" in line and "dCp" in line:
                    continue

                values = line.split()
                if len(values) == 4:
                    try:
                        x, y, z, dCp = map(float, values)
                        slices[current_slice].append((x, dCp))
                    except ValueError:
                        continue 

        # find cl and lift per unit span
        lift_slices = []
        y_values = sorted(slices.keys())

        for y_slice, data in slices.items():
            if len(data) < 2:
                continue

            data.sort(key=lambda point: point[0])
            x_vals, dCp_vals = np.array(data)[:, 0], np.array(data)[:, 1]

            dx_vals = np.diff(x_vals)

            chord_length = np.max(x_vals) - np.min(x_vals)
            if chord_length == 0:
                continue 

            # compute lift coefficient
            Cl = np.sum(-dCp_vals[:-1] * dx_vals) / chord_length

            # compute lift per unit span
            L_prime = Cl * q_inf * chord_length

            lift_slices.append((L_prime))
        return lift_slices

    def compute_x_cop_percentage(file_path):
        slices = {}
        current_slice = None

        with open(file_path, "r", errors="ignore") as file:
            for line in file:
                match = re.match(r"BLOCK Cut_\d+_at_Y:_([\d\.\-]+)", line)
                if match:
                    current_slice = float(match.group(1))
                    slices[current_slice] = []
                    continue

                if current_slice is None:
                    continue

                if "Mach" in line or "Alpha" in line or "Beta" in line:
                    continue

                if "x" in line and "y" in line and "z" in line and "dCp" in line:
                    continue

                values = line.split()
                if len(values) == 4:
                    try:
                        x, y, z, dCp = map(float, values)
                        slices[current_slice].append((x, dCp))
                    except ValueError:
                        continue

        results = []

        for y_slice, data in slices.items():
            if not data:
                continue 

            data = np.array(data) 
            x_vals, dCp_vals = data[:, 0], data[:, 1]

            x_min, x_max = np.min(x_vals), np.max(x_vals)
            chord_length = x_max - x_min

            if chord_length == 0:
                continue 

            # compute x position of cp as a percentage
            sum_dCp = np.sum(dCp_vals)
            if sum_dCp != 0:
                x_CoP = np.sum(x_vals * dCp_vals) / sum_dCp
                x_CoP_percentage = (x_CoP - x_min) / chord_length
            else:
                x_CoP_percentage = np.nan 

            results.append(x_CoP_percentage)

        return results

    def extract_wing_geometry(file_path):
        slices = {}
        current_slice = None

        with open(file_path, "r", errors="ignore") as file:
            for line in file:
                match = re.match(r"BLOCK Cut_\d+_at_Y:_([\d\.\-]+)", line)
                if match:
                    current_slice = float(match.group(1))
                    slices[current_slice] = []
                    continue

                if current_slice is None:
                    continue

                if "Mach" in line or "Alpha" in line or "Beta" in line or "x" in line:
                    continue

                values = line.split()
                if len(values) == 4:
                    try:
                        x, y, z, dCp = map(float, values)
                        slices[current_slice].append(x)
                    except ValueError:
                        continue

        if not slices:
            return None

        y_values = sorted(slices.keys())
        y_root, y_tip = y_values[0], y_values[-1]
        outerboard_span = y_tip-y_root

        #find root and tip x coords
        x_root = np.array(slices[y_root])
        x_tip = np.array(slices[y_tip])

        # find le and te x positions
        x_LE_root, x_TE_root = np.min(x_root), np.max(x_root)
        x_LE_tip, x_TE_tip = np.min(x_tip), np.max(x_tip)

        c_root = x_TE_root - x_LE_root # root chord
        c_tip = x_TE_tip - x_LE_tip # tip chord
        taper_ratio = c_tip / c_root if c_root != 0 else np.nan

        x_quarter_root = x_LE_root + 0.25 * c_root
        x_quarter_tip = x_LE_tip + 0.25 * c_tip

        delta_x_quarter = x_quarter_tip - x_quarter_root
        delta_y = y_tip - y_root
        sweep_quarter_chord = np.arctan(delta_x_quarter / delta_y) if delta_y != 0 else np.nan

        return {
            "Root Chord": c_root,
            "Tip Chord": c_tip,
            "Taper Ratio": taper_ratio,
            "1/4 Chord Sweep Angle": sweep_quarter_chord,
            "Outerboard Span": outerboard_span,
            "x_root": x_quarter_root,
            "x_tip": x_quarter_tip,
            "y_root": y_root,
            "y_tip": y_tip
        }

    # Aerofoil Functions

    def load_aerofoil(aerofoil_file):
        # loads airfoil coordinates from a .dat file
        with open(aerofoil_file, 'r') as file:
            lines = file.readlines()

        coordinates = []
        for line in lines:
            parts = line.strip().split()
            if len(parts) == 2:
                try:
                    x, y = float(parts[0]), float(parts[1])
                    coordinates.append((x, y))
                except ValueError:
                    continue 

        coordinates = np.array(coordinates)
        coordinates = coordinates[np.argsort(coordinates[:, 0])] # sort by x coordinates
        coordinates = coordinates[:-1]
        return coordinates

    def find_spar_heights(aerofoil, LE_spar, TE_spar):
        front_spar_x = LE_spar
        rear_spar_x = TE_spar

        front_spar_y = abs(aerofoil[np.abs(aerofoil[:, 0] - front_spar_x).argmin()+1][1] - aerofoil[np.abs(aerofoil[:, 0] - front_spar_x).argmin()][1])
        rear_spar_y = abs(aerofoil[np.abs(aerofoil[:, 0] - rear_spar_x).argmin()+1][1] - aerofoil[np.abs(aerofoil[:, 0] - rear_spar_x).argmin()][1])

        t_over_c = max(row[1] for row in aerofoil) - min(row[1] for row in aerofoil)
        max_index = max(range(len(aerofoil)), key=lambda i: aerofoil[i][1])
        max_thickness_location = aerofoil[max_index,0]
        avg_spar_height = (front_spar_y + rear_spar_y)/2
        return front_spar_y, rear_spar_y, avg_spar_height, t_over_c, max_thickness_location

    # Structural Functions

    def compute_shear_center(f_height, r_height, c_height, t_spar, t_skin, wingbox_chord, max_thickness_location_wingbox, A_spar_cap, I_spar_cap):
        sy = 1
        f_halfchord = max_thickness_location_wingbox
        r_halfchord = wingbox_chord - max_thickness_location_wingbox
        A_cap = A_spar_cap

        # lengths of sections
        L12 = f_height
        L23 = np.sqrt(f_halfchord**2+(0.5*(c_height-f_height))**2)
        L34 = np.sqrt(r_halfchord**2+(0.5*(c_height-r_height))**2)
        L45 = r_height
        L56 = L34
        L61 = L23

        #y-pos of points
        y1 = -f_height/2
        y2 = f_height/2
        y3 = c_height/2
        y4 = r_height/2
        y5 = -r_height/2
        y6 = -c_height/2

        # areas of booms
        A1 = abs(L12*t_spar*(2+y2/y1))/6 + abs(L61*t_skin*(2+y6/y1))/6 + A_cap
        A2 = abs(L12*t_spar*(2+y1/y2))/6 + abs(L23*t_skin*(2+y3/y2))/6 + A_cap
        A3 = abs(L23*t_skin*(2+y2/y3))/6 + abs(L34*t_skin*(2+y4/y3))/6
        A4 = abs(L34*t_skin*(2+y3/y4))/6 + abs(L45*t_spar*(2+y5/y4))/6 + A_cap
        A5 = abs(L45*t_spar*(2+y4/y5))/6 + abs(L56*t_skin*(2+y6/y5))/6 + A_cap
        A6 = abs(L56*t_skin*(2+y5/y6))/6 + abs(L61*t_skin*(2+y1/y6))/6

        A_tot = A1+A2+A3+A4+A5+A6

        Ycg=(A1*y1+A2*y2+A3*y3+A4*y4+A5*y5+A6*y6)/A_tot

        #x-pos of points
        x1 = -f_halfchord
        x2 = -f_halfchord
        x3 = 0
        x4 = r_halfchord
        x5 = r_halfchord
        x6 = 0
        Xcg=(A1*x1+A2*x2+A3*x3+A4*x4+A5*x5+A6*x6)/A_tot

        # contribution of the spars 
        I_f_spar = (t_spar * (f_height ** 3)) / 12
        I_r_spar = (t_spar * (r_height ** 3)) / 12

        # contribution of the top and bottom skins 
        f_theta = np.arctan((c_height-f_height)/(2*f_halfchord))
        r_theta = np.arctan((c_height-r_height)/(2*r_halfchord))
        I_f_skin = ((((np.sin(f_theta) ** 2) * t_skin * (L23 ** 3)) / 12 + L23*t_skin) * ((c_height + f_height) / 4) ** 2) * 2
        I_r_skin = ((((np.sin(r_theta) ** 2) * t_skin * (L34 ** 3)) / 12 + L34*t_skin) * ((c_height + r_height) / 4) ** 2) * 2

        #spar caps at corners of wingbox
        I_f_spar_caps = I_spar_cap + A_spar_cap * (f_height / 2) ** 2
        I_r_spar_caps = I_spar_cap + A_spar_cap * (r_height / 2) ** 2

        # total second moment of area
        Ixx = I_f_spar + I_r_spar + I_f_skin + I_r_skin + I_f_spar_caps + I_r_spar_caps

        # open shear flows
        q1openx= 0
        q2openx=q1openx-(sy/Ixx)*A2*(y2-Ycg)
        q3openx=q2openx-(sy/Ixx)*A3*(y3-Ycg)
        q4openx=q3openx-(sy/Ixx)*A4*(y4-Ycg)
        q5openx=q4openx-(sy/Ixx)*A5*(y5-Ycg)
        q6openx=q5openx-(sy/Ixx)*A6*(y6-Ycg)

        # calculate p
        a1 = - y2 + y1
        b1 = x2 - x1
        c1 = x1*-a1 - y1*b1
        p1 = abs((a1*0+b1*0+c1)/np.sqrt(a1*a1 + b1*b1))
        a2 = - y3 + y2
        b2 = x3 - x2
        c2 = x2*-a2 - y2*b2
        p2 = abs((a2*0+b2*0+c2)/np.sqrt(a2*a2 + b2*b2))
        a3 = - y4 + y3
        b3 = x4 - x3
        c3 = x3*-a3 - y3*b3
        p3 = abs((a3*0+b3*0+c3)/np.sqrt(a3*a3 + b3*b3))
        a4 = - y5 + y4
        b4 = x5 - x4
        c4 = x4*-a4 - y4*b4
        p4 = abs((a4*0+b4*0+c4)/np.sqrt(a4*a4 + b4*b4))
        a5 = - y6 + y5
        b5 = x6 - x5
        c5 = x5*-a5 - y5*b5
        p5 = abs((a5*0+b5*0+c5)/np.sqrt(a5*a5 + b5*b5))
        a6 = - y1 + y6
        b6 = x1 - x6
        c6 = x6*-a6 - y6*b6
        p6 = abs((a6*0+b6*0+c6)/np.sqrt(a6*a6 + b6*b6))

        enclosed_area = (f_height+c_height)*f_halfchord*0.5 + (r_height+c_height)*r_halfchord*0.5

        q_s0 = -(q1openx*L12*p1+q2openx*L23*p2+q3openx*L34*p3+q4openx*L45*p4+q5openx*L56*p5+q6openx*L61*p6)/(2*enclosed_area)

        q1 = q1openx + q_s0
        q2 = q2openx + q_s0
        q3 = q3openx + q_s0
        q4 = q4openx + q_s0
        q5 = q5openx + q_s0
        q6 = q6openx + q_s0

        # take moments about 0,0
        x_shear = ((q1*p1*L12+q2*p2*L23+q3*p3*L23+q4*p4*L45+q5*p5*L56+q6*p6*L61)/sy + f_halfchord) / (f_halfchord+r_halfchord)
        return x_shear
    
    def compute_bending_moment(total_loads, rib_pitch):
        n = len(total_loads)
        bending_moment = np.zeros(n)

        for i in range(n):
            moment = 0
            for j in range(i, n):
                moment += total_loads[j] * ((j - i) * rib_pitch)
            bending_moment[i] = moment

        return bending_moment

    def compute_twisting_moment(aerodynamic_loads, CoPs, wingbox_chord, LE_spar, wingbox_fraction, xshear):
        n = len(aerodynamic_loads)
        twisting_moment = np.zeros(n)
        wing_chord = np.zeros(n)

        for i in range(n-2, -1, -1):
            wing_chord[i] = wingbox_chord[i] / wingbox_fraction
            dx = CoPs[i]*wing_chord[i] - LE_spar * wing_chord[i] - wingbox_chord[i] * xshear
            twisting_moment[i] = aerodynamic_loads[i] * dx + twisting_moment[i+1]

        return twisting_moment

    def compute_second_moment(spar_thickness, wingbox_height, skin_thickness, wingbox_chord, I_spar_cap, I_stringer, A_spar_cap, A_stringer, num_stringers, spar_height):
        # spars
        I_spar = (spar_thickness * (wingbox_height ** 3)) / 12

        # skins - angled thin wall approximation
        theta_skin = np.arctan((wingbox_height - spar_height) / wingbox_chord)
        I_skin = (((np.sin(theta_skin) ** 2) * skin_thickness * ((np.cos(theta_skin) * wingbox_chord / 2) ** 3)) / 12  + ((np.cos(theta_skin) * wingbox_chord)*skin_thickness*0.5) * ((wingbox_height + spar_height) / 4) ** 2) * 2 #Hexagonal wingbox

        # stringers
        if num_stringers > 0:
            I_stringers =  I_stringer + A_stringer * ((wingbox_height + spar_height)/ 4) ** 2
        else:
            I_stringers = 0

        # spar caps
        I_spar_caps = I_spar_cap + A_spar_cap * (wingbox_height / 2) ** 2

        # total
        I_total = 2 * I_spar + 4 * I_spar_caps + 2 * I_skin + 2 * num_stringers * I_stringers
        return I_total

    def compute_torsional_constant(spar_thickness, wingbox_height, skin_thickness, wingbox_chord):

        # lengths of each section of the closed wingbox
        a = wingbox_chord # top and bottom skins
        b = wingbox_height # front and rear spars
        t1 = skin_thickness
        t2 = spar_thickness

        J = 2*(2*t1*t2*(a-t2)**2*(b-t1)**2)/(a*t2+b*t1-t1**2-t2**2)
        return J

    def skin_second_moment(skin_thickness, wingbox_chord, I_stringer, num_stringers, spar_height, wingbox_height):
        # angled thin wall approximation
        theta_skin = np.arctan((wingbox_height - spar_height) / wingbox_chord)
        I_skin = (((np.sin(theta_skin) ** 2) * skin_thickness * ((np.cos(theta_skin) * wingbox_chord / 2) ** 3)) / 12) * 2

        # Contribution of stringers
        if num_stringers > 0:
            I_stringer =  I_stringer 
        else:
            I_stringer = 0

        # total
        I_skin_total = I_skin + num_stringers * I_stringer
        return I_skin_total

    # Buckling Functions
    def calculate_pcrit(skinEI, rib_pitch, wingbox_chord):
        k = 0.53*np.exp(-0.505 * rib_pitch / wingbox_chord)
        pcrit = (3.1416 ** 2 * skinEI) / ((rib_pitch * k) ** 2)
        return pcrit

    # Mass Functions
    def wingbox_mass_inertia(density, root_chord, taper_ratio, y_root, y_tip, t_over_c, skin_thickness, spar_thickness, A_spar_cap, A_stringer, num_stringers, sweep_angle, LE_spar, TE_spar, wingbox_fraction):
        def chord(y):
            return root_chord * (1 - (1 - taper_ratio) * (y - y_root) / (y_tip - y_root))
        
        def section_properties(y):
            c = chord(y)
            h = t_over_c * c
            spar_height = avg_spar_height_fraction * h
            theta_skin = np.arctan((h - spar_height) / c)
            A_skin = 2 * (c*np.cos(theta_skin)) * skin_thickness
            A_spars = 2 * spar_height * spar_thickness / np.cos(sweep_angle)
            A_LE = skin_thickness*np.pi*LE_spar*c/wingbox_fraction
            A_TE = skin_thickness*2*(1-TE_spar)*c/wingbox_fraction
            A_total = A_skin + A_spars + 4 * A_spar_cap / np.cos(sweep_angle) + num_stringers * A_stringer * 2 / np.cos(sweep_angle) + A_TE + A_LE
            return A_total
        
        mass = 0
        I_x = 0
        I_y = 0
        I_z = 0
        y_moment = 0

        i = 0
        dy = (y_tip - y_root) / elements
        # integration along span for mass and moment
        for y in np.linspace(y_root, y_tip, elements):
            A = section_properties(y)
            dmass[i] = density * A * dy
            mass += dmass[i]
            y_moment += dmass[i] * y
            i += 1
        
        y_cg = y_moment / mass
        x_cg = -(x_root * (y_cg-y_root)/(y_tip-y_root) + x_tip * (1 - (y_cg-y_root)/(y_tip-y_root))) # approximation
        z_cg = 0 # wing height above nose - can be defined elsewhere in the system

        i = 0
        # integration along span for inertias
        for y in np.linspace(y_root, y_tip, elements):
            A = section_properties(y)
            dmass[i] = density * A * dy
            x = x_root * (y-y_root)/(y_tip-y_root) + x_tip * (1 - (y-y_root)/(y_tip-y_root))
            I_x += dmass[i] * (y-y_cg) ** 2
            I_y += dmass[i] * (-x-x_cg) ** 2
            I_z = I_x
            i += 1
        return mass, I_x, I_y, I_z, x_cg, y_cg, z_cg, dmass

    #extract wing geometry from pressure data file
    wing_geometry = extract_wing_geometry(file_path)
    sweep_angle = wing_geometry["1/4 Chord Sweep Angle"]
    root_chord = wing_geometry["Root Chord"] * wingbox_fraction * np.cos(sweep_angle)
    taper_ratio = wing_geometry["Taper Ratio"]
    outerboardspan = wing_geometry["Outerboard Span"]

    tip_deflection = np.zeros(2)
    tip_twist = np.zeros(2)
    buckled = np.zeros(2)

    for load_case in range(2):
        #Define aerodynamic force arrays
        aerodynamic_loads = compute_lift_forces(file_path, q_inf)
        elements = len(aerodynamic_loads)
        projected_rib_pitch = outerboardspan / (elements- 1)
        aerodynamic_loads = [load_factor[load_case] * num * projected_rib_pitch for num in aerodynamic_loads]
        rib_pitch = projected_rib_pitch / np.cos(sweep_angle)

        # initialise arrays
        beam_length= outerboardspan / np.cos(sweep_angle)
        wingbox_height=np.zeros(elements)
        wingbox_chord = np.zeros(elements)
        spar_height = np.zeros(elements)
        skin_force = np.zeros(elements)
        dmass = np.zeros(elements)
        EI = np.zeros(elements)
        I = np.zeros(elements)
        J = np.zeros(elements)
        GJ = np.zeros(elements)
        z = np.linspace(0, beam_length, elements)

        aerofoil = load_aerofoil(aerofoil_file)
        spar_heights = find_spar_heights(aerofoil, LE_spar, TE_spar)
        t_over_c = spar_heights[3]
        avg_spar_height_fraction = spar_heights[2] / t_over_c
        max_thickness_location = spar_heights[4]

        #compute EI and GJ
        for i in range(0, elements):
            wingbox_chord[i] = root_chord * (1 - ((1 - taper_ratio) * (z[i] / (beam_length))))
            wingbox_height[i] = t_over_c * (wingbox_chord[i] / wingbox_fraction)
            spar_height[i] = wingbox_height[i] * avg_spar_height_fraction
            I[i] = compute_second_moment(spar_thickness=spar_thickness, wingbox_height=wingbox_height[i], skin_thickness=skin_thickness, wingbox_chord=wingbox_chord[i], I_spar_cap=I_spar_cap, A_spar_cap=A_spar_cap, A_stringer=A_stringer, I_stringer=I_stringer, num_stringers=num_stringers, spar_height=spar_height[i])
            EI[i] = E * I[i]
            J[i] = compute_torsional_constant(spar_thickness=spar_thickness, wingbox_height=wingbox_height[i],skin_thickness=skin_thickness,wingbox_chord=wingbox_chord[i])
            GJ[i] = G * J[i]

        #Compute Mass Properties
        x_root = wing_geometry["x_root"]
        x_tip = wing_geometry["x_tip"]
        y_root = wing_geometry["y_root"]
        y_tip = wing_geometry["y_tip"]
        mass, I_x, I_y, I_z, x_cg, y_cg, z_cg, dmass = wingbox_mass_inertia(density, root_chord, taper_ratio, y_root, y_tip, t_over_c, skin_thickness, spar_thickness, A_spar_cap, A_stringer, num_stringers, sweep_angle, LE_spar, TE_spar, wingbox_fraction)

        # Define mass force arrays
        secondary_masses = np.zeros(elements)
        inertial_loads = -(dmass + secondary_masses) * 9.81 * load_factor[load_case]
        total_loads = inertial_loads + aerodynamic_loads
        M = -compute_bending_moment(total_loads, rib_pitch)

        dz = z[1] - z[0] #spanwise step size

        # initialise variables
        Ve = np.zeros_like(z)
        zdef = np.zeros_like(z)
        ztotaldef = np.zeros_like(z)
        length = 0
        theta = np.zeros_like(z)
        upper_skin_stress = np.zeros_like(z)

        #Central difference deflection
        for i in range(0, len(z) - 1):
            d2Ve_dz2 = -M[i] / EI[i] #bending equation
            Ve[i + 1] = 2 * Ve[i] - Ve[i-1] + d2Ve_dz2 * (dz ** 2)

        CoPs = compute_x_cop_percentage(file_path)

        CoPs = np.array(CoPs, dtype=float)
        wingbox_chord = np.array(wingbox_chord, dtype=float)

        # Compute shear centre as a percentage of wingbox chord from the front spar
        xshear = compute_shear_center(spar_heights[0], spar_heights[1], spar_heights[3], spar_thickness, skin_thickness, TE_spar-LE_spar, max_thickness_location-LE_spar,A_spar_cap, I_spar_cap)
        Q = compute_twisting_moment(aerodynamic_loads, CoPs, wingbox_chord, LE_spar, wingbox_fraction, xshear) 

        # Initialize twist angle with forward euler
        alpha_e = np.zeros_like(z)
        dAlpha_dz = Q[0] / GJ[0]
        alpha_e[0] = 0
        alpha_e[1] = 0 + dz * dAlpha_dz

        for i in range(1, len(z) - 1):
            dAlpha_dz = Q[i] / GJ[i] # Apply torsion equation
            alpha_e[i + 1] = alpha_e[i - 1] + dAlpha_dz * dz * 2 # Integrate to get twist angle with central difference

        
        for k in range(1, len(z)):
            theta[k] = (Ve[k] - Ve[k - 1]) / dz #Backward Euler to find slope
            ztotaldef[k] = (beam_length/len(z))-((beam_length/len(z))*np.cos(theta[k])) + ztotaldef[k-1]
            zdef[k] = z[k] - ztotaldef[k] # find spanwise deflection
            length = length+np.sqrt((Ve[k]-Ve[k-1])*(Ve[k]-Ve[k-1])+(zdef[k]-zdef[k-1])*(zdef[k]-zdef[k-1]))
        
        #BUCKLING
        #Determine forces in skin segments, compare against critical buckling force
        for i in range(1, len(z)):
            upper_skin_stress[i] = -M[i-1] * wingbox_height[i] * 0.5 / I[i-1]
            skin_force = upper_skin_stress[i] * (num_stringers * A_stringer + (wingbox_chord[i] + wingbox_chord[i-1]) * 0.5 * skin_thickness)
            skinEI = E * skin_second_moment(skin_thickness, (wingbox_chord[i] + wingbox_chord[i-1]) * 0.5, I_stringer=I_stringer, num_stringers=num_stringers, spar_height=spar_height[i], wingbox_height=wingbox_height[i])
            pcrit = calculate_pcrit(skinEI, rib_pitch, (wingbox_chord[i] + wingbox_chord[i-1]) * 0.5)
            if skin_force > pcrit / safetyfactor:
                print("panel ", i, " buckled")
                buckled[load_case] = 1

        import csv

        if load_case == 0:
            with open("def_twist_pos.csv", "w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(["Deflection", "Twist"])
                writer.writerows(zip(Ve, alpha_e))

        if load_case == 1:
            with open("def_twist_neg.csv", "w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(["Deflection", "Twist"])
                writer.writerows(zip(Ve, alpha_e))

        tip_deflection[load_case] = Ve[elements-1]
        tip_twist[load_case] = alpha_e[elements-1]

    return [x_cg,y_cg,z_cg,I_x,I_y,I_z,mass,np.max(buckled), tip_deflection[0], tip_deflection[1], np.max(abs(tip_twist)), TE_spar]

#compute_wing_mass(0.08, 8, 0.01, 0.02, 210, 0.38, [3.5,-1], 0.1, 0.7)
