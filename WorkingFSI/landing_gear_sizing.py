def landing_gear_sizing(ground_clearance, MTOW, nlg_distance_from_cg, mlg_distance_from_cg, sink_speed, z_cg):
    
    #Imports

    import numpy as np
    import math

    # SI Units for all parameters

    #Define material properties and safety factor e.g 300M steel
    safety_factor = 1.5
    E = 205e9
    density = 7870
    yield_strength = 1.931e9

    #Define gear and wheel parameters
    #ground_clearance = 1.98
    wheel_diameter = 0.7 
    mlg_length = ground_clearance - wheel_diameter/2
    nlg_length = ground_clearance - wheel_diameter/2
    gear_efficient_factor = 0.7 # Empirical
    lm = mlg_distance_from_cg #1.5 # Distance between CG and main gear
    ln = nlg_distance_from_cg #8 # Distance between CG and nose gear
    weight_distribution_factor = ln / (lm + ln) # Fraction of total weight supported by main gear
    num_wheels = 2 # wheels per gear for single bogey
    tire_stroke = 0.1 * wheel_diameter
    shock_strut_mass_factor = 0.7 # Fraction of total structure mass for the shock strut, empirical NASA paper. Other masses include braces and fittings

    #Define landing parameters
    sink_speed = 3.8
    touchdown_aoa = 3 #degrees
    touchdown_aoa = np.radians(touchdown_aoa)
    MTOW = MTOW * 9.81
    lift = MTOW
    design_reaction_factor = 1.2 # Widely used
    ns = 0.47 #shock efficient factor
    nt = 0.8 #tire efficiecny factor

    #Define braking parameters
    ax_over_g = 0.5 # Braking acceleration as a fraction of g

    #Compute stroke lengths
    mlg_stroke = (MTOW*sink_speed**2/(2*9.81)+(MTOW - lift)*tire_stroke)/(ns*design_reaction_factor*MTOW+nt*design_reaction_factor*MTOW - MTOW + lift)
    nlg_stroke = 1

    #Define load cases
    #load case 0 - level landing conditions (25.479), load case 1 - side load conditions (25.485), load case 2 - one gear landing conditions (25.483)
    afm = [0.75, 0.5, 1] #Axial force multiplier
    sfm = [0.25, 0.8, 0] #Side force multiplier
    dfm = [0.4, 0, 0.25] #Drag force multiplier

    mlg_mass = np.zeros(3)
    nlg_mass = np.zeros(3)
    mlg_cg = np.zeros(3)
    nlg_cg = np.zeros(3)

    for load_case in range(3):
        #Compute main gear forces (axial and transverse from FAR 25) impact forces
        def compute_mlg_forces():
            mlg_axial_force = weight_distribution_factor * MTOW * (sink_speed * sink_speed / 9.81 + mlg_stroke * math.cos(touchdown_aoa)) / (gear_efficient_factor * mlg_stroke * math.cos(touchdown_aoa))
            mlg_drag_force = dfm[load_case] * mlg_axial_force 
            mlg_side_force = sfm[load_case] * mlg_axial_force 
            mlg_axial_force = mlg_axial_force * afm[load_case]
            mlg_transverse_force = np.sqrt(mlg_side_force**2 + mlg_drag_force**2)
            return mlg_axial_force, mlg_transverse_force
        mlg_axial_force, mlg_transverse_force = compute_mlg_forces()

        #Compute nose gear forces (axial and transverse from FAR 25) low speed braking acceleration forces
        def compute_nlg_forces():
            nlg_axial_force = MTOW * (lm + z_cg * ax_over_g)/(lm + ln) # axial force for the nose isn't affected by load cases as there is only one nose gear
            nlg_drag_force = dfm[load_case] * nlg_axial_force 
            nlg_side_force = sfm[load_case] * nlg_axial_force 
            nlg_transverse_force = np.sqrt(nlg_side_force**2 + nlg_drag_force**2)
            return nlg_axial_force, nlg_transverse_force
        nlg_axial_force, nlg_transverse_force = compute_nlg_forces()

        #Compute moment due to transverse forces
        mlg_transverse_moment = mlg_transverse_force * mlg_length * 0.5
        nlg_transverse_moment = nlg_transverse_force * nlg_length * 0.5

        #Buckling with calibrated k value
        def compute_radius(axial_force, strut_length):
            K = 2
            I = safety_factor * axial_force * pow(strut_length * K, 2) / (E * 9.8696)
            radius = np.sqrt(np.sqrt(64 * I / 3.1416)) / 2
            return radius, I

        mlg_radius, mlg_I = compute_radius(mlg_axial_force, mlg_length)
        nlg_radius, nlg_I = compute_radius(nlg_axial_force, nlg_length)

        #Stress calculation
        def compute_stress(axial_force, radius, transverse_moment, transverse_force, I):
            area = np.pi * radius ** 2
            normal_stress = axial_force / area + transverse_moment/I * radius
            shear_stress = transverse_force / area
            von_mises = np.sqrt(normal_stress**2 + 3*shear_stress**2)
            return von_mises
        mlg_stress = compute_stress(mlg_axial_force, mlg_radius, mlg_transverse_moment, mlg_transverse_force, mlg_I)
        nlg_stress = compute_stress(nlg_axial_force, nlg_radius, nlg_transverse_moment, nlg_transverse_force, nlg_I)

        while mlg_stress*safety_factor > yield_strength:
            mlg_radius += 0.005
            I = (mlg_radius*2)**4 * np.pi/64
            mlg_stress = compute_stress(mlg_axial_force, mlg_radius, mlg_transverse_moment, mlg_transverse_force, I)

        while nlg_stress*safety_factor > yield_strength:
            nlg_radius += 0.005
            I = (nlg_radius*2)**4 * np.pi/64
            nlg_stress = compute_stress(nlg_axial_force, nlg_radius, nlg_transverse_moment, nlg_transverse_force, I)

        #Wheel assembly calculation - empirical Currey
        def compute_wheel_assembly_mass(axial_force): 
            wheel_assembly_mass = (0.230893*(axial_force*2.2*wheel_diameter/(1000*3.28))**0.8482) / 2.205
            return wheel_assembly_mass
        mlg_wa_mass = compute_wheel_assembly_mass(mlg_axial_force)
        nlg_wa_mass = compute_wheel_assembly_mass(nlg_axial_force)

        #Mass and CoG calculation
        def compute_mass_properties(radius, length, wheel_assembly_mass):
            area = np.pi * radius ** 2
            mass = area * length * density / shock_strut_mass_factor + wheel_assembly_mass
            cg_pos = (0.5 * length * (area * length * density / shock_strut_mass_factor) + wheel_assembly_mass * length) / mass # postion of cg from hinge
            return mass, cg_pos

        mlg_mass[load_case], mlg_cg[load_case] = compute_mass_properties(mlg_radius, mlg_length, mlg_wa_mass)
        nlg_mass[load_case], nlg_cg[load_case] = compute_mass_properties(nlg_radius, nlg_length, nlg_wa_mass)

    mlg_mass_max = max(mlg_mass)
    nlg_mass_max = max(nlg_mass)
    mlg_max_index = np.argmax(mlg_mass)
    nlg_max_index = np.argmax(nlg_mass)
    mlg_cg_max = mlg_cg[mlg_max_index]
    nlg_cg_max = nlg_cg[nlg_max_index]
    return [mlg_mass_max, nlg_mass_max, mlg_cg_max, nlg_cg_max]