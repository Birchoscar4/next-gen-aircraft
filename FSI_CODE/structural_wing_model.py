def structural_wing_model():
    import numpy as np
    import re
    #import matplotlib.pyplot as plt

    file_path = "pressure.slc"
    aerofoil_file = "0714sc.dat"
    rho = 1.2256 # kg/m3 - air
    v = 272 # m/s
    q_inf = 0.5 * rho * v **2
    LE_spar = 0.08
    TE_spar = 0.7
    wingbox_fraction = TE_spar - LE_spar
    safetyfactor=1.5
    #t_over_c = 0.14

    #Material Properties
    E = 6.8947573e10 # Pa
    v = 0.4
    G = E / (2*(1+v))
    density = 1550 # kg/m3

    #Structural Parameters
    spar_cap_radius = 0.03
    I_spar_cap=np.pi*(2*spar_cap_radius)**4/64
    A_spar_cap=np.pi*(spar_cap_radius**2)
    num_stringers = 8
    stringer_thickness = 0.005
    stringer_depth = 0.04
    stringer_width = 0.02
    A_stringer = stringer_thickness*stringer_depth + 2*(stringer_thickness*stringer_width)
    I_stringer = stringer_depth**3*stringer_thickness/12 + 2*(stringer_thickness**3*stringer_width/12 + stringer_thickness*stringer_width*(stringer_thickness+stringer_depth)**2/4)
    skin_thickness = 0.003
    spar_thickness = 0.005
    rib_thickness = 0.002

    MTOW = 60000
    load_factor = [2.1 + (24000 / (MTOW * 2.2 + 10000)), -1.0]
    if load_factor[0] < 2.5:
        load_factor[0] = 2.5

    # Aerodynamics functions
    def compute_lift_forces(file_path, q_inf):
        slices = {}
        current_slice = None

        with open(file_path, "r", errors="ignore") as file:
            for line in file:
                # Detect block headers
                match = re.match(r"BLOCK Cut_\d+_at_Y:_([\d\.\-]+)", line)
                if match:
                    current_slice = float(match.group(1))  # Extract Y position
                    slices[current_slice] = []  # Initialize slice data
                    continue

                # Ensure a valid slice is initialized before reading data
                if current_slice is None:
                    continue

                # Ignore metadata lines
                if "Mach" in line or "Alpha" in line or "Beta" in line:
                    continue

                # Ignore column headers
                if "x" in line and "y" in line and "z" in line and "dCp" in line:
                    continue

                # Parse data rows (x, y, z, dCp)
                values = line.split()
                if len(values) == 4:
                    try:
                        x, y, z, dCp = map(float, values)
                        slices[current_slice].append((x, dCp))  # Store (x, dCp)
                    except ValueError:
                        continue  # Skip lines that cannot be parsed correctly

        # Compute Cl and lift per unit span (L') for each slice
        lift_slices = []
        y_values = sorted(slices.keys())  # Sorted spanwise positions

        for y_slice, data in slices.items():
            if len(data) < 2:
                continue  # Need at least two points to compute lift

            # Sort data by x (in case they are unordered)
            data.sort(key=lambda point: point[0])
            x_vals, dCp_vals = np.array(data)[:, 0], np.array(data)[:, 1]

            # Compute delta x (spacing between adjacent points)
            dx_vals = np.diff(x_vals)

            # Compute chord length (max difference in x) at each slice
            chord_length = np.max(x_vals) - np.min(x_vals)
            if chord_length == 0:
                continue  # Avoid division by zero

            # Compute sectional lift coefficient Cl
            Cl = np.sum(-dCp_vals[:-1] * dx_vals) / chord_length

            # Compute lift per unit span (L')
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
                continue  # Skip empty slices

            data = np.array(data)  # Convert to numpy array
            x_vals, dCp_vals = data[:, 0], data[:, 1]

            # Compute chord length dynamically
            x_min, x_max = np.min(x_vals), np.max(x_vals)
            chord_length = x_max - x_min

            if chord_length == 0:
                continue  # Avoid division by zero

            # Compute x_CoP using weighted averages
            sum_dCp = np.sum(dCp_vals)
            if sum_dCp != 0:
                x_CoP = np.sum(x_vals * dCp_vals) / sum_dCp
                x_CoP_percentage = (x_CoP - x_min) / chord_length   # Convert to percentage
            else:
                x_CoP_percentage = np.nan  # Handle cases with zero sum

            results.append(x_CoP_percentage)

        return results

    def extract_wing_geometry(file_path):
        slices = {}
        current_slice = None

        with open(file_path, "r", errors="ignore") as file:
            for line in file:
                # Detect block headers for each slice
                match = re.match(r"BLOCK Cut_\d+_at_Y:_([\d\.\-]+)", line)
                if match:
                    current_slice = float(match.group(1))  # Extract Y position
                    slices[current_slice] = []  # Initialize slice data
                    continue

                # Ensure a valid slice is initialized before reading data
                if current_slice is None:
                    continue

                # Ignore metadata lines and column headers
                if "Mach" in line or "Alpha" in line or "Beta" in line or "x" in line:
                    continue

                # Parse data rows (x, y, z, dCp)
                values = line.split()
                if len(values) == 4:
                    try:
                        x, y, z, dCp = map(float, values)
                        slices[current_slice].append(x)  # Store x-coordinates only
                    except ValueError:
                        continue  # Skip lines that cannot be parsed correctly

        # Identify root and tip slices
        if not slices:
            return None

        y_values = sorted(slices.keys())  # Sorted y-values
        y_root, y_tip = y_values[0], y_values[-1]
        outerboard_span = y_tip-y_root

        # Get x-coordinates at root and tip
        x_root = np.array(slices[y_root])
        x_tip = np.array(slices[y_tip])

        # Compute leading and trailing edge locations
        x_LE_root, x_TE_root = np.min(x_root), np.max(x_root)
        x_LE_tip, x_TE_tip = np.min(x_tip), np.max(x_tip)

        # Compute geometric parameters
        c_root = x_TE_root - x_LE_root  # Root chord
        c_tip = x_TE_tip - x_LE_tip  # Tip chord
        taper_ratio = c_tip / c_root if c_root != 0 else np.nan

        # Compute 1/4 chord positions
        x_quarter_root = x_LE_root + 0.25 * c_root
        x_quarter_tip = x_LE_tip + 0.25 * c_tip

        # Compute 1/4 chord sweep angle
        delta_x_quarter = x_quarter_tip - x_quarter_root
        delta_y = y_tip - y_root
        sweep_quarter_chord = np.arctan(delta_x_quarter / delta_y) if delta_y != 0 else np.nan

        return {
            "Root Chord": c_root,
            "Tip Chord": c_tip,
            "Taper Ratio": taper_ratio,
            "1/4 Chord Sweep Angle": sweep_quarter_chord,
            "Outerboard Span": outerboard_span,
            "y_root": y_root,
            "y_tip": y_tip
        }

    # Aerofoil Functions
    def load_aerofoil(aerofoil_file):
        """Loads airfoil coordinates from a .dat file"""
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
                    continue  # Skip invalid lines

        coordinates = np.array(coordinates)
        coordinates = coordinates[np.argsort(coordinates[:, 0])]  # Sort by x-coordinates
        coordinates = coordinates[:-1]
        return coordinates

    def find_spar_heights(aerofoil, LE_spar, TE_spar): #NOT WORKING RIGHT
        front_spar_x = LE_spar
        rear_spar_x = TE_spar

        front_spar_y = aerofoil[np.abs(aerofoil[:, 0] - front_spar_x).argmin()+1][1] - aerofoil[np.abs(aerofoil[:, 0] - front_spar_x).argmin()][1]
        rear_spar_y = aerofoil[np.abs(aerofoil[:, 0] - front_spar_x).argmin()+1][1] - aerofoil[np.abs(aerofoil[:, 0] - rear_spar_x).argmin()][1]

        t_over_c = max(row[1] for row in aerofoil) - min(row[1] for row in aerofoil)
        max_index = max(range(len(aerofoil)), key=lambda i: aerofoil[i][1])
        max_thickness_location = aerofoil[max_index,0]
        avg_spar_height = (front_spar_y + rear_spar_y)/2
        return front_spar_y, rear_spar_y, avg_spar_height, t_over_c, max_thickness_location

    def compute_rib_area(aerofoil, wingbox_chord, wingbox_fraction):
        perim_scale = wingbox_chord / wingbox_fraction
        area_scale = perim_scale * perim_scale
        for i in range(0, len(aerofoil)):
            dx = aerofoil[i*2+2, 0] -aerofoil[i*2, 0]
            #print (dx)


        area =1
        return area

    # Structural Functions
    def compute_shear_center(f_height, r_height, c_height, t_spar, t_skin, wingbox_chord, max_thickness_location_wingbox, A_spar_cap, I_spar_cap): #ALMOST WORKING
        sy = 1
        f_halfchord = 21#max_thickness_location_wingbox
        r_halfchord = 1#wingbox_chord - max_thickness_location_wingbox
        A_cap = A_spar_cap
        f_height = 1#
        r_height = 1#
        c_height =1#

        # lengths of sections
        L12 = f_height
        L23 = np.sqrt(f_halfchord**2+(0.5*(c_height-f_height))**2)
        L34 = np.sqrt(r_halfchord**2+(0.5*(c_height-r_height))**2)
        L45 = r_height
        L56 = L34
        L61 = L23

        # Areas of sections
        A1 = (t_skin*L61)/2 + (t_spar*L12)/2 + A_cap
        A2 = (t_spar*L12)/2 + (t_skin*L23)/2 + A_cap
        A3 = (t_skin*L23)/2 + (t_skin*L34)/2
        A4 = (t_skin*L34)/2 + (t_spar*L45)/2 + A_cap
        A5 = (t_spar*L45)/2 + (t_skin*L56)/2 + A_cap
        A6 = (t_skin*L56)/2 + (t_skin*L61)/2

        A_tot = A1+A2+A3+A4+A5+A6

        #y-pos of points
        y1 = -f_height/2
        y2 = f_height/2
        y3 = c_height/2
        y4 = r_height/2
        y5 = -r_height/2
        y6 = -c_height/2
        Ycg=(A1*y1+A2*y2+A3*y3+A4*y4+A5*y5+A6*y6)/A_tot

        #x-pos of points
        x1 = -f_halfchord
        x2 = -f_halfchord
        x3 = 0
        x4 = r_halfchord
        x5 = r_halfchord
        x6 = 0
        Xcg=(A1*x1+A2*x2+A3*x3+A4*x4+A5*x5+A6*x6)/A_tot

        # Contribution of the spars 
        I_f_spar = (t_spar * (f_height ** 3)) / 12
        I_r_spar = (t_spar * (r_height ** 3)) / 12

        # Contribution of the top and bottom skins 
        f_theta = np.arctan((c_height-f_height)/(2*f_halfchord))
        r_theta = np.arctan((c_height-r_height)/(2*r_halfchord))
        I_f_skin = ((((np.sin(f_theta) ** 2) * t_skin * (L23 ** 3)) / 12  + L23*t_skin) * ((c_height + f_height) / 4) ** 2) * 2
        I_r_skin = ((((np.sin(r_theta) ** 2) * t_skin * (L34 ** 3)) / 12  + L34*t_skin) * ((c_height + r_height) / 4) ** 2) * 2

        #Spar caps at corners of wingbox
        I_f_spar_caps = I_spar_cap + A_spar_cap * (f_height / 2) ** 2
        I_r_spar_caps = I_spar_cap + A_spar_cap * (r_height / 2) ** 2

        # Total second moment of area
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
        x_shearabs =(q1*p1*L12+q2*p2*L23+q3*p3*L23+q4*p4*L45+q5*p5*L56+q6*p6*L61)/sy
        #print(x_shearabs)
        x_shear = ((q1*p1*L12+q2*p2*L23+q3*p3*L23+q4*p4*L45+q5*p5*L56+q6*p6*L61)/sy + f_halfchord) / (f_halfchord+r_halfchord)
        return x_shear
    
    def calculate_bending_moment(total_loads, rib_pitch):
        n = len(total_loads)
        bending_moment = np.zeros(n)

        # Compute bending moment at each force application point
        for i in range(n):
            moment = 0
            for j in range(i, n):
                moment += total_loads[j] * ((j - i) * rib_pitch)
            bending_moment[i] = moment

        return bending_moment

    def calculate_twisting_moment(aerodynamic_loads, CoPs, wingbox_chord, LE_spar, wingbox_fraction):
        n = len(aerodynamic_loads)
        twisting_moment = np.zeros(n)
        wing_chord = np.zeros(n)

        # Compute twisting moment at each point
        for i in range(n-2, -1, -1):
            wing_chord[i] = wingbox_chord[i] / wingbox_fraction
            dx = CoPs[i]*wing_chord[i] - LE_spar * wing_chord[i] - wingbox_chord[i] * 0.5
            twisting_moment[i] = aerodynamic_loads[i] * dx + twisting_moment[i+1]

        return twisting_moment

    def wingbox_second_moment(spar_thickness, wingbox_height, skin_thickness, wingbox_chord, I_spar_cap, I_stringer, A_spar_cap, A_stringer, num_stringers, spar_height):
        # Contribution of the spars 
        I_spar = (spar_thickness * (wingbox_height ** 3)) / 12

        # Contribution of the top and bottom skins 
        #I_skin = 0 * (wingbox_chord * (skin_thickness ** 3)) / 12  + (wingbox_chord*skin_thickness) * (wingbox_height / 2) ** 2 # rectangular wingbox
        theta_skin = np.arctan((wingbox_height - spar_height) / wingbox_chord)
        I_skin = (((np.sin(theta_skin) ** 2) * skin_thickness * ((np.cos(theta_skin) * wingbox_chord / 2) ** 3)) / 12  + ((np.cos(theta_skin) * wingbox_chord)*skin_thickness*0.5) * ((wingbox_height + spar_height) / 4) ** 2) * 2 #Hexagonal wingbox

        # Contribution of stringers
        if num_stringers > 0:
            I_stringers =  I_stringer + A_stringer * ((wingbox_height + spar_height)/ 4) ** 2
        else:
            I_stringers = 0

        #Spar caps at corners of wingbox
        I_spar_caps = I_spar_cap + A_spar_cap * (wingbox_height / 2) ** 2

        # Total second moment of area
        I_total = 2 * I_spar + 4 * I_spar_caps + 2 * I_skin + 2 * num_stringers * I_stringers
        return I_total

    def wingbox_polar_moment(spar_thickness, wingbox_height, skin_thickness, wingbox_chord, I_spar_cap, I_stringer, A_spar_cap, A_stringer, num_stringers, spar_height):
        # Contribution of the spars to Ix (bending about z-axis)
        I_spar_x = (spar_thickness * (wingbox_height ** 3)) / 12
        
        # Contribution of the spars to Iy (bending about y-axis)
        I_spar_y = (wingbox_height * (spar_thickness ** 3)) / 12 + (wingbox_height * spar_thickness) * (wingbox_chord / 2) ** 2
        
        # Contribution of the top and bottom skins to Ix
        theta_skin = np.arctan((wingbox_height - spar_height) / wingbox_chord)
        I_skin_x = (((np.sin(theta_skin) ** 2) * skin_thickness * ((np.cos(theta_skin) * wingbox_chord / 2) ** 3)) / 12  + ((np.cos(theta_skin) * wingbox_chord)*skin_thickness*0.5) * ((wingbox_height + spar_height) / 4) ** 2) * 2 #Hexagonal wingbox
        
        # Contribution of the top and bottom skins to Iy
        I_skin_y = (skin_thickness * (wingbox_chord ** 3)) / 12 
        
        # Spar caps at the corners of the wingbox
        I_spar_caps_x = I_spar_cap + A_spar_cap * (wingbox_height / 2) ** 2
        I_spar_caps_y = I_spar_cap + A_spar_cap * (wingbox_chord / 2) ** 2
        
        # Total second moment of area about x and y axes
        I_x_total = 2 * I_spar_x + 4 * I_spar_caps_x + 2 * I_skin_x
        I_y_total = 2 * I_spar_y + 4 * I_spar_caps_y + 2 * I_skin_y 
        
        # Polar moment of inertia (J = Ix + Iy)
        J_total = I_x_total + I_y_total
        
        return J_total

    def skin_second_moment(skin_thickness, wingbox_chord, I_stringer, num_stringers, spar_height, wingbox_height):
        theta_skin = np.arctan((wingbox_height - spar_height) / wingbox_chord)
        I_skin = (((np.sin(theta_skin) ** 2) * skin_thickness * ((np.cos(theta_skin) * wingbox_chord / 2) ** 3)) / 12) * 2
        #I_skin = (wingbox_chord * skin_thickness**3) / 12
        # Contribution of stringers
        if num_stringers > 0:
            I_stringer =  I_stringer 
        else:
            I_stringer = 0

        # Total second moment of area
        I_skin_total =  I_skin + num_stringers * I_stringer
        #print("iskin",I_skin_total)
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
            A_TE = skin_thickness*2*TE_spar*c/wingbox_fraction
            A_total = A_skin + A_spars + 4 * A_spar_cap / np.cos(sweep_angle) + num_stringers * A_stringer * 2 / np.cos(sweep_angle) + A_TE + A_LE
            return A_total
        
        mass = 0
        I_x = 0
        i = 0
        dy = (y_tip - y_root) / elements  # Integration step
        for y in np.linspace(y_root, y_tip, elements):
            A = section_properties(y)
            dmass[i] = density * A * dy
            mass += dmass[i]
            I_x += dmass[i] * y ** 2
            i += 1

        return mass, I_x, dmass

    #Extract wing geometry from pressure data file
    wing_geometry = extract_wing_geometry(file_path)
    sweep_angle = wing_geometry["1/4 Chord Sweep Angle"]
    root_chord = wing_geometry["Root Chord"] / wingbox_fraction * np.cos(sweep_angle)
    taper_ratio = wing_geometry["Taper Ratio"]
    outerboardspan = wing_geometry["Outerboard Span"]

    #Define aerodynamic force arrays
    aerodynamic_loads = compute_lift_forces(file_path, q_inf)
    elements = len(aerodynamic_loads)
    projected_rib_pitch = outerboardspan / (elements- 1)
    aerodynamic_loads = [load_factor[0] * num * projected_rib_pitch for num in aerodynamic_loads]
    rib_pitch = projected_rib_pitch / np.cos(sweep_angle)
    centres_of_pressure = np.zeros(elements)

    # Define given discrete data along the span
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
    z = np.linspace(0, beam_length, elements)  # Beam spanwise positions from root to tip

    aerofoil = load_aerofoil(aerofoil_file)
    spar_heights = find_spar_heights(aerofoil, LE_spar, TE_spar)
    t_over_c = spar_heights[3]
    avg_spar_height_fraction = spar_heights[2] / t_over_c
    max_thickness_location = spar_heights[4]
    #print(max_thickness_location)


    #Compute EI and GJ with geometric parameters along span
    for i in range(0, elements):
        wingbox_chord[i] = root_chord * (1 - ((1 - taper_ratio) * (z[i] / (beam_length))))
        wingbox_height[i] = t_over_c * (wingbox_chord[i] / wingbox_fraction)
        spar_height[i] = wingbox_height[i] * avg_spar_height_fraction
        I[i] = wingbox_second_moment(spar_thickness=spar_thickness, wingbox_height=wingbox_height[i], skin_thickness=skin_thickness, wingbox_chord=wingbox_chord[i], I_spar_cap=I_spar_cap, A_spar_cap=A_spar_cap, A_stringer=A_stringer, I_stringer=I_stringer, num_stringers=num_stringers, spar_height=spar_height[i])
        EI[i] = E * I[i]
        J[i] = wingbox_polar_moment(spar_thickness=spar_thickness, wingbox_height=wingbox_height[i], skin_thickness=skin_thickness, wingbox_chord=wingbox_chord[i], I_spar_cap=I_spar_cap, A_spar_cap=A_spar_cap, A_stringer=A_stringer, I_stringer=I_stringer, num_stringers=num_stringers, spar_height=spar_height[i])
        GJ[i] = G * J[i]

    #Compute Mass Properties
    y_root = wing_geometry["y_root"]
    y_tip = wing_geometry["y_tip"]
    #for i in range(1, elements):
        #rib_area = compute_rib_area(aerofoil, wingbox_chord[i], wingbox_fraction)
    mass, I_x, dmass = wingbox_mass_inertia(density, root_chord, taper_ratio, y_root, y_tip, t_over_c, skin_thickness, spar_thickness, A_spar_cap, A_stringer, num_stringers, sweep_angle, LE_spar, TE_spar, wingbox_fraction)
    # Define mass force arrays
    secondary_masses = np.zeros(elements)
    inertial_loads = -(dmass + secondary_masses) * 9.81 * load_factor[0]
    total_loads = inertial_loads + aerodynamic_loads
    M = -calculate_bending_moment(total_loads, rib_pitch)

    dz = z[1] - z[0]  # Spanwise step size

    # Solve for deflection using finite difference integration
    Ve = np.zeros_like(z)  # Initialize deflection
    zdef = np.zeros_like(z)  # Initialize deflection
    ztotaldef = np.zeros_like(z)  # Initialize deflection
    length = 0
    theta = np.zeros_like(z)  # Initialize slope (dVe/dz)
    upper_skin_stress = np.zeros_like(z)  # Initialize skin stress

    #Central difference deflection
    for i in range(0, len(z) - 1):
        d2Ve_dz2 = -M[i] / EI[i]  # Apply bending equation
        #theta[i + 1] = theta[i] + d2Ve_dz2 * dz  # Integrate to get slope - forward euler
        #Ve[i + 1] = Ve[i] + theta[i] * dz  # Integrate slope to get deflection - forward euler
        Ve[i + 1] = 2 * Ve[i] - Ve[i-1] + d2Ve_dz2 * (dz ** 2)

    CoPs = compute_x_cop_percentage(file_path)
    # Ensure variables are NumPy arrays for element-wise operations
    CoPs = np.array(CoPs, dtype=float)
    wingbox_chord = np.array(wingbox_chord, dtype=float)
    Q = calculate_twisting_moment(aerodynamic_loads, CoPs, wingbox_chord, LE_spar, wingbox_fraction)  # Example torque distribution
    #plt.plot(z,Q)

    # Solve for twist angle using numerical integration
    # Initialize twist angle with forward euler
    alpha_e = np.zeros_like(z)  
    dAlpha_dz = Q[0] / GJ[0]
    alpha_e[0] = 0
    alpha_e[1] = 0 + dz * dAlpha_dz

    for i in range(1, len(z) - 1):
        dAlpha_dz = Q[i] / GJ[i]  # Apply torsion equation
        #alpha_e[i + 1] = alpha_e[i] + dAlpha_dz * dz #Forward euler
        alpha_e[i + 1] = alpha_e[i - 1] + dAlpha_dz * dz * 2 # Integrate to get twist angle with central difference

    #Find spanwise deflection
    for k in range(1, len(z)):
        theta[k] = (Ve[k] - Ve[k - 1]) / dz #Backward Euler to find slope
        ztotaldef[k] = (beam_length/len(z))-((beam_length/len(z))*np.cos(theta[k])) + ztotaldef[k-1]
        zdef[k] = z[k] - ztotaldef[k]   # find spanwise deflection
        length = length+np.sqrt((Ve[k]-Ve[k-1])*(Ve[k]-Ve[k-1])+(zdef[k]-zdef[k-1])*(zdef[k]-zdef[k-1]))
        
    #BUCKLING
    #Determine forces in skin segments, compare against critical buckling force
    for i in range(1, len(z)):
        upper_skin_stress[i] = -M[i-1] * wingbox_height[i] * 0.5 / I[i-1]
        skin_force = upper_skin_stress[i] * (num_stringers * A_stringer + (wingbox_chord[i] + wingbox_chord[i-1]) * 0.5 * skin_thickness)
        #print("Skin force:   ", skin_force)
        skinEI = E * skin_second_moment(skin_thickness, (wingbox_chord[i] + wingbox_chord[i-1]) * 0.5, I_stringer=I_stringer, num_stringers=num_stringers, spar_height=spar_height[i], wingbox_height=wingbox_height[i])
        pcrit = calculate_pcrit(skinEI, rib_pitch, (wingbox_chord[i] + wingbox_chord[i-1]) * 0.5)
        #print("Critical Load:", pcrit)
        if skin_force > pcrit / safetyfactor:
            print("panel ", i, " buckled")


    print("Mass:", mass,"kg")
    print("Mass Moment of Inertia:", I_x,"kgm^2")
    print("Tip deflection (m):", Ve[elements-1])

    # Plot results
    """
    plt.figure(figsize=(12, 5))

    plt.subplot(1, 2, 1)
    plt.plot(zdef, Ve)
    plt.xlim(0,beam_length)
    plt.ylim(0,2)
    plt.xlabel("Spanwise Position z")
    plt.ylabel("Deflection (m)")
    plt.title("Deflection Along Span")
    plt.grid()

    plt.subplot(1, 2, 2)
    plt.plot(z, alpha_e, color='r')
    plt.xlabel("Spanwise Position z")
    plt.ylabel("Twist Angle (rad)")
    plt.title("Twist Along Span")
    plt.grid()

    plt.show()
    """

    import csv

    with open("def_twist.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(zip(Ve, alpha_e))

    #xshear = compute_shear_center(spar_heights[0], spar_heights[1], spar_heights[3], spar_thickness, skin_thickness, TE_spar-LE_spar, max_thickness_location-LE_spar,A_spar_cap, I_spar_cap)
    #xshear = compute_shear_center(1, 1, 1, spar_thickness, skin_thickness, 2, 0.5,A_spar_cap, I_spar_cap,)
    #print("xshear",xshear)
    return Ve[elements-1]
