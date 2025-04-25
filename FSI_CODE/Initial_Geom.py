def build_initial_geometry():
    import openvsp as vsp
    import os
    import subprocess
    import platform
    import numpy as np

    # Initialize OpenVSP
    def initialize_vsp():
        vsp.ClearVSPModel()
        print("VSP model initialized.")

    # Add a wing geometry with specified parameters
    def add_wing(wing_name, sections):
        wing_id = vsp.AddGeom("WING", "")
        vsp.SetGeomName(wing_id, wing_name)

        # Add multi-sections and set parameters
        for i, section in enumerate(sections):
            if i < 3 + num_ribs:
                vsp.InsertXSec(wing_id, i+1, vsp.XS_FOUR_SERIES)
                vsp.SetParmVal(wing_id, "Span", f"XSec_{i+1}", section["span"])
                vsp.SetParmVal(wing_id, "Root_Chord", f"XSec_{i+1}", section["root_chord"])
                vsp.SetParmVal(wing_id, "Tip_Chord", f"XSec_{i+1}", section["tip_chord"])
                vsp.SetParmVal(wing_id, "Sweep", f"XSec_{i+1}", section["sweep"])
                vsp.SetParmVal(wing_id, "Sweep_Location", f"XSec_{i+1}", section["sweep_loc"])
                vsp.Update()
        vsp.CutXSec(wing_id, num_ribs+1)
        return wing_id

    """
    FSI Structural Coupling with Geometry"
    """

    def find_tip_chord(root_chord, taper_ratio, span, rib_pitch, i):
        tip_chord = root_chord * (1 - ((1 - taper_ratio) * (rib_pitch * i) / span))
        return tip_chord

    # Blend wing sections
    def blend_wing_sections(wing_id, blend_params):
        for section, params in blend_params.items():
            for param, value in params.items():
                vsp.SetParmVal(wing_id, param, section, value)
            vsp.Update()

    def update_aerofoils(wing_id, aerofoil_file_1, aerofoil_file_2):

        # Get the cross-section surface (XSecSurf) for the wing
        xsec_surf = vsp.GetXSecSurf(wing_id, 0)  # Assuming the main XSecSurf is at index 0

        # Loop through the specified cross-sections and update the airfoil
        for xsec_index in range(11):  # Inclusive of end_index
            if xsec_index <= 2:
            # Set the cross-section shape to XS_FILE_AIRFOIL
                vsp.ChangeXSecShape(xsec_surf, xsec_index, vsp.XS_FILE_AIRFOIL)

                # Get the cross-section (XSec) at the current index
                xsec_id = vsp.GetXSec(xsec_surf, xsec_index)
        
                vsp.ReadFileAirfoil(xsec_id, aerofoil_file_1)
            else:
                vsp.ChangeXSecShape(xsec_surf, xsec_index, vsp.XS_FILE_AIRFOIL)

                xsec_id = vsp.GetXSec(xsec_surf, xsec_index)

                vsp.ReadFileAirfoil(xsec_id, aerofoil_file_2)
            
            # Verify the cross-section type
            xsec_type = vsp.GetXSecShape(xsec_id)
            if xsec_type != vsp.XS_FILE_AIRFOIL:
                print(f"Error: XSec {xsec_index} is not of type XS_FILE_AIRFOIL. Skipping...")
                continue
            print(f"Updated airfoil for XSec {xsec_index}")

        # Update the model to reflect the changes
        vsp.Update()


    # Add vertical tail
    def add_vertical_tail(tail_name, position, params):
        tail_id = vsp.AddGeom("WING", "")
        vsp.SetGeomName(tail_id, tail_name)

        # Set position and parameters
        for param, value in position.items():
            vsp.SetParmVal(tail_id, param, "XForm", value)
        for param, value in params.items():
            vsp.SetParmVal(tail_id, param, "XSec_1", value)
        vsp.Update()

        return tail_id

    """
    MATHEMATICAL RELATION FOR SLAT AND ELEVON PLACEMENT" --- STILL NEED TO FIGURE OUT A MATHEMATICAL FORMULA FOR EVEN POSITIONING  
    """

    # Define the chordwise extent of each slat and elevon
    slat_chord_extent = 0.2  # Slats cover 20% of the chord
    elevon_chord_extent = 0.3  # Elevons cover 30% of the chord

    # Add elevons to the wing                           
    def add_elevons(wing_id, num_elevons=2, start_u=0.05, end_u=0.95):
        print("Adding Elevons")
        spacing = (end_u - start_u) / (num_elevons + 1)
        global elevon_parm
        elevon_parm = {}

        for i in range(num_elevons):
            # Calculate UStart and UEnd for elevons (TE)
            elevon_start_u = start_u + (i + 1) * spacing
            elevon_end_u = elevon_start_u + spacing

            # Ensure elevon does not exceed the maximum limit spanwise.
            if elevon_end_u > end_u:
                elevon_end_u = end_u

            elevon_id = vsp.AddSubSurf(wing_id, vsp.SS_CONTROL)  # Add a control surface
            vsp.SetParmVal(elevon_id, "LE_Flag", "SS_Control", 0)
            vsp.SetParmVal(elevon_id, "Surf_Type", "SS_Control", 2)  # Set as Elevon (type 2)
            vsp.SetParmVal(elevon_id, "UStart", "SS_Control", elevon_start_u)  # Start in chordwise direction
            vsp.SetParmVal(elevon_id, "UEnd", "SS_Control", elevon_end_u)  # End in chordwise direction
            vsp.SetParmVal(elevon_id, "Length_C_Start", "SS_Control", elevon_chord_extent)
            
            print(f"Added Elevon {i+1} with ID: {elevon_id}")
            elevon_parm[f'elevon_name{i+1}'] = elevon_id
            print(elevon_parm)
            vsp.Update()

            #cs_group_container_id = vsp.FindContainer("VSPAEROSettings", 0)

            #vsp.SetParmVal(vsp.FindParm(cs_group_container_id, "Surf_" + elevon_id + "_0_Gain", "ControlSurfaceGroup_0"), 1)
            #vsp.SetParmVal(vsp.FindParm(cs_group_container_id, "Surf_" + elevon_id + "_1_Gain", "ControlSurfaceGroup_0"), -1)
            #vsp.SetParmVal(vsp.FindParm(cs_group_container_id, "DeflectionAngle", "ControlSurfaceGroup_0"), 10.0)
            #vsp.Update()


    def add_slats(wing_id, num_slats=3, start_u=0.17, end_u=0.9):
        spacing = (end_u - start_u) / (num_slats + 1)
        global slat_parm
        slat_parm = {}
        for i in range(num_slats):
                # Calculate UStart and UEnd for elevons (TE)
                slat_start_u = start_u + (i + 1) * spacing
                slat_end_u = slat_start_u + spacing
            
                # Ensure that Slat does not exceed limit of 0.9 spanwise
                if slat_end_u > end_u:
                    slat_end_u = end_u
                
                # Add control surface
                global slat_id
                slat_id = vsp.AddSubSurf(wing_id, vsp.SS_CONTROL)
                vsp.SetParmVal(slat_id, "LE_Flag", "SS_Control", 2)
                vsp.SetParmVal(slat_id, "Surf_Type", "SS_Control", 2)  # Set as Elevon (type 2)
                vsp.SetParmVal(slat_id, "UStart", "SS_Control", slat_start_u)  # Start in chordwise direction
                vsp.SetParmVal(slat_id, "UEnd", "SS_Control", slat_end_u)  # End in chordwise direction
                vsp.SetParmVal(slat_id, "Length_C_Start", "SS_Control", slat_chord_extent)
                print(f"Added Slats {i+1} with ID: {slat_id}")
                slat_parm[f'slat_name{i+1}'] = slat_id
                print(slat_parm)
        return slat_id, slat_parm
    
    def add_rudder(tail_id):
        rudder_id = vsp.AddSubSurf(tail_id, vsp.SS_CONTROL)
        global rudder_parm
        rudder_parm = {}
        u_start = 0.35
        u_end = 0.65

        vsp.SetParmVal(rudder_id, "UStart", "SS_Control", u_start)
        vsp.SetParmVal(rudder_id, "UEnd", "SS_Control", u_end)
        vsp.SetParmVal(rudder_id, "Surf_Type", "SS_Control", 2)
        vsp.SetParmVal(rudder_id, "Length_C_Start", "SS_Control", 0.25)
        print(f"Added Rudder with ID: {rudder_id}")
        rudder_parm = rudder_id
        print (rudder_parm)
        return rudder_id
    
    # Add engines with propellers
    def add_engines_with_propellers():
        print("Adding Engines with Propeller Fans (Manually Symmetrical)")
        y_locations_eng = [4.0, -4.0]  # Y-locations for the engines (positive and negative for symmetry)
        y_locations_prop = [4.0, -4.0]  # Y-locations for the propellers (symmetry)

        for i, y_loc in enumerate(y_locations_eng):
            engine_id = vsp.AddGeom("POD", "")  # Add a POD geometry (engine nacelle)
            vsp.SetParmVal(engine_id, "X_Rel_Location", "XForm", 13.5)  # Position along the wing span
            vsp.SetParmVal(engine_id, "Y_Rel_Location", "XForm", y_loc)  # Position vertically (mirrored)
            vsp.SetParmVal(engine_id, "Z_Rel_Location", "XForm", 0.2)  # Position horizontally
            vsp.SetParmVal(engine_id, "FineRatio", "Design", 12)  # Set engine diameter
            vsp.SetParmVal(engine_id, "Length", "Design", 4.0)  # Set engine length

            # Add a Propeller Fan
            prop_id = vsp.AddGeom("PROP", engine_id)  # Add a propeller to the engine
            vsp.SetParmVal(prop_id, "NumBlade", "Design", 10.0)
            vsp.SetParmVal(prop_id, "Diameter", "Design", 4.0)  # Set propeller diameter
            vsp.SetParmVal(prop_id, "X_Rel_Location", "XForm", 16)  # Position relative to the engine
            vsp.SetParmVal(prop_id, "Y_Rel_Location", "XForm", y_locations_prop[i])  # Position relative to the engine
            vsp.SetParmVal(prop_id, "Z_Rel_Location", "XForm", 0.0)  # Position relative to the engine
            print(f"Added Engine {i+1} with Propeller Fan, Engine ID: {engine_id}, Propeller ID: {prop_id}")

    # Save and export the model
    def save_and_export_model(filename):
        vsp.Update()
        vsp.WriteVSPFile(f"{filename}.vsp3")
        vsp.ExportFile(f"{filename}.stl", vsp.SET_ALL, vsp.EXPORT_STL)
        print(f"Model saved as '{filename}.vsp3' and exported as '{filename}.stl'.")

    # Open the exported file in the default viewer
    def open_file_in_viewer(file_path):
        system = platform.system()
        if system == 'Windows':
            os.startfile(file_path)
        elif system == 'Darwin':
            subprocess.run(['open', file_path])
        elif system == 'Linux':
            subprocess.run(['xdg-open', file_path])
        else:
            print(f"Unsupported OS: {system}")

    outerwing_span = 12.5
    num_ribs = 10
    outer_root_chord = 5.576
    outer_tip_chord = 2.51
    taper_ratio = outer_tip_chord / outer_root_chord
    root_chord_array = np.zeros(num_ribs)
    root_chord_array[1] = outer_root_chord
    rib_pitch = (outerwing_span / (num_ribs - 1))

    for k in range(2, num_ribs):
        root_chord_array[k] = find_tip_chord(root_chord=outer_root_chord, taper_ratio=taper_ratio, span=outerwing_span, rib_pitch=rib_pitch, i=k)

    def group_cs():
    
        e_index = vsp.CreateVSPAEROControlSurfaceGroup()
        elevon_group = "Elevon_Group"
        vsp.SetVSPAEROControlGroupName( elevon_group, e_index)
        available_cs_vec = vsp.GetAvailableCSNameVec( e_index)
        print (available_cs_vec)
        a_ind_vec = [1, 2, 3, 4]
        vsp.AddSelectedToCSGroup(a_ind_vec, e_index)
        active_aileron_vec = vsp.GetActiveCSNameVec(e_index)
        aileron_added_set = {available_cs_vec[i-1] for i in a_ind_vec}
        active_aileron_set = set(active_aileron_vec)
        if active_aileron_set != aileron_added_set:
            print("Incorrect Elevon set")
        print("The elevon group contains:")
        print(active_aileron_set)

        vsp.Update()

        s_index = vsp.CreateVSPAEROControlSurfaceGroup()
        slat_group = "Slat_Group"
        vsp.SetVSPAEROControlGroupName( slat_group, s_index)
        slat_ind_vec = [5, 6, 7, 8, 9, 10]
        vsp.AddSelectedToCSGroup(slat_ind_vec, s_index)
        active_slat_vec = vsp.GetActiveCSNameVec(s_index)
        slat_added_set = {available_cs_vec[i-1] for i in slat_ind_vec}
        active_slat_set = set(active_slat_vec)
        if active_slat_set != slat_added_set:
            print("Incorrect Slat set")
        print("The slat group contains:")
        print(active_slat_set)

        vsp.Update()

        r_index = vsp.CreateVSPAEROControlSurfaceGroup()
        rudder_group = "Rudder_Group"
        vsp.SetVSPAEROControlGroupName( rudder_group, r_index)
        rudder_ind_vec = [11, 12]
        vsp.AddSelectedToCSGroup(rudder_ind_vec, r_index)
        active_rudder_vec = vsp.GetActiveCSNameVec(r_index)
        rudder_added_set = {available_cs_vec[i-1] for i in rudder_ind_vec}
        active_rudder_set = set(active_rudder_vec)
        if active_rudder_set != rudder_added_set:
            print("Incorrect rudder set")
        print("The rudder group contains:")
        print(active_rudder_set)

        cs_group_container_id = vsp.FindContainer("VSPAEROSettings", 0)

        vsp.Update()

        # Add slats
       
        print(type(slat_parm))
        slat1 = slat_parm['slat_name1']
        slat2 = slat_parm['slat_name2']
        slat3 = slat_parm['slat_name3']
        print(type(slat_id))

        elevon1 = elevon_parm['elevon_name1']

        cs_group_container_id = vsp.FindContainer("VSPAEROSettings", 0)

        vsp.Update()

        #Slat Deflection

        vsp.SetParmVal(vsp.FindParm(cs_group_container_id, "Surf_" + slat1 + "_0_Gain", "ControlSurfaceGroup_1"), 1 )
        vsp.SetParmVal(vsp.FindParm(cs_group_container_id, "Surf_" + slat1 + "_1_Gain", "ControlSurfaceGroup_1"), -1 )
        vsp.SetParmVal(vsp.FindParm(cs_group_container_id, "DeflectionAngle", "ControlSurfaceGroup_1"), 10.0)
        slat_deflection_angle = vsp.GetParmVal(cs_group_container_id, "DeflectionAngle", "ControlSurfaceGroup_1")
        print(f"Slat Deflection Angle: {slat_deflection_angle} degrees")

        #Aileron deflection
        vsp.SetParmVal(vsp.FindParm(cs_group_container_id, "Surf_" + elevon1 + "_0_Gain", "ControlSurfaceGroup_0"), 1 )
        vsp.SetParmVal(vsp.FindParm(cs_group_container_id, "Surf_" + elevon1 + "_1_Gain", "ControlSurfaceGroup_0"), -1 )
        vsp.SetParmVal(vsp.FindParm(cs_group_container_id, "DeflectionAngle", "ControlSurfaceGroup_0"), 10.0)
        elevon_deflection_angle = vsp.GetParmVal(cs_group_container_id, "DeflectionAngle", "ControlSurfaceGroup_0")
        print(f"Elevon Deflection Angle: {elevon_deflection_angle} degrees")

        #Rudder deflection
        vsp.SetParmVal(vsp.FindParm(cs_group_container_id, "Surf_" + rudder_parm + "_0_Gain", "ControlSurfaceGroup_2"), 1 )
        vsp.SetParmVal(vsp.FindParm(cs_group_container_id, "Surf_" + rudder_parm + "_1_Gain", "ControlSurfaceGroup_2"), -1 )
        vsp.SetParmVal(vsp.FindParm(cs_group_container_id, "DeflectionAngle", "ControlSurfaceGroup_2"), 0.5)
        rudder_deflection_angle = vsp.GetParmVal(cs_group_container_id, "DeflectionAngle", "ControlSurfaceGroup_2")
        print(f"Rudder Deflection Angle: {rudder_deflection_angle} degrees")

        



    # Main function to build the aircraft
    def build_aircraft(outerwing_span, root_chord_array, rib_pitch, num_ribs):
        initialize_vsp()

        # Define wing sections
        wing_sections = [
            {"span": 7.5, "root_chord": 17.684, "tip_chord": 6.97, "sweep": 45, "sweep_loc": 0},
            {"span": 2, "root_chord": 6.97, "tip_chord": 5.576, "sweep": 40, "sweep_loc": 0},
        ]

        for x in range(2, num_ribs):
            wing_sections.append({
                "span": rib_pitch,
                "root_chord": root_chord_array[x - 1],
                "tip_chord": root_chord_array[x],
                "sweep": 25,
                "sweep_loc": 0
            })
        wing_id = add_wing("WING", wing_sections)
        print(wing_sections)

        # Define wing blending parameters
        blend_params = {
            "XSec_0": {"OutLEMode": 5, "OutTEMode": 5, "OutLEStrength": 0.3, "OutTEStrength": 0.0},
            "XSec_1": {"InLEMode": 2, "OutLEMode": 4, "InLEStrength": 1.07, "OutLEStrength": 1.00,
                    "InTEMode": 3, "OutTEMode": 5, "InTEStrength": 1.00, "OutTEStrength": 0.00},
            "XSec_2": {"InLEMode": 2, "OutLEMode": 4, "InLEStrength": 3.30, "OutLEStrength": 0.00,
                    "InTEMode": 3, "OutTEMode": 5, "InTEStrength": 3.00, "OutTEStrength": 1.34},
            "XSec_3": {"InLEMode": 2, "OutLEMode": 4, "InLEStrength": 1.00, "OutLEStrength": 1.00,
                    "InTEMode": 3, "OutTEMode": 5, "InTEStrength": 1.00, "OutTEStrength": 1.00}
        }
        blend_wing_sections(wing_id, blend_params)

        aerofoil_file_1 = "C:\\23012.dat"
        aerofoil_file_2 = "C:\\0714sc.dat"
        
        update_aerofoils(wing_id, aerofoil_file_1, aerofoil_file_2)


        # Add vertical tail
        tail_position = {"X_Rel_Location": 14.00, "Y_Rel_Location": -8.60656, "X_Rel_Rotation": 90}
        tail_params = {"Span": 6.0, "Sweep": 30.0, "Sweep_Location": 0.0, "Root_Chord": 6.0, "Tip_Chord": 1.5, "ThickChord": 0.1}
        tail_id = add_vertical_tail("VerticalTail", tail_position, tail_params)
       
        #slat_parm = add_slats(wing_id)
        #print(slat_parm)

        #Control surface deflections
        #slat1 = (slat_parm[f'slat_name{1}'])
        #slat2 = (slat_parm[f'slat_name{2}'])
        #slat3 = (slat_parm[f'slat_name{3}'])

        cs_group_container_id = vsp.FindContainer("VSPAEROSettings", 0)

        vsp.Update()

        #Slat1 Deflection

        #vsp.SetParmVal(vsp.FindParm(cs_group_container_id, "Surf_" + slat1 + "_0_Gain", "ControlSurfaceGroup_1"), 1 )
        #vsp.SetParmVal(vsp.FindParm(cs_group_container_id, "Surf_" + slat1 + "_1_Gain", "ControlSurfaceGroup_1"), -1 )
        #vsp.SetParmVal(vsp.FindParm(cs_group_container_id, "DeflectionAngle", "ControlSurfaceGroup_1"), 10.0)
        #slat_deflection_angle = vsp.GetParmVal(cs_group_container_id, "DeflectionAngle", "ControlSurfaceGroup_1")
        #print(f"Slat Deflection Angle: {slat_deflection_angle} degrees")

        # Add elevons
        add_elevons(wing_id)

        # Add slats
       
        add_slats(wing_id)
        
        # Add Rudder
        add_rudder(tail_id)


        # Add engines with propellers
        add_engines_with_propellers()
        group_cs()
        #deflections()
        # Save and export the model
        save_and_export_model("Prelim_3")

        # Open the STL file in the default viewer
        # open_file_in_viewer("Prelim_1.stl")
        return wing_id, wing_sections #slat_parm
    
    #slat_parm = add_slats(wing_id)
    

    #def group_cs():
    
        e_index = vsp.CreateVSPAEROControlSurfaceGroup()
        elevon_group = "Elevon_Group"
        vsp.SetVSPAEROControlGroupName( elevon_group, e_index)
        available_cs_vec = vsp.GetAvailableCSNameVec( e_index)
        print (available_cs_vec)
        a_ind_vec = [1, 2, 3, 4]
        vsp.AddSelectedToCSGroup(a_ind_vec, e_index)
        active_aileron_vec = vsp.GetActiveCSNameVec(e_index)
        aileron_added_set = {available_cs_vec[i-1] for i in a_ind_vec}
        active_aileron_set = set(active_aileron_vec)
        if active_aileron_set != aileron_added_set:
            print("Incorrect Elevon set")
        print("The elevon group contains:")
        print(active_aileron_set)

        vsp.Update()

        s_index = vsp.CreateVSPAEROControlSurfaceGroup()
        slat_group = "Slat_Group"
        vsp.SetVSPAEROControlGroupName( slat_group, s_index)
        slat_ind_vec = [5, 6, 7, 8, 9, 10]
        vsp.AddSelectedToCSGroup(slat_ind_vec, s_index)
        active_slat_vec = vsp.GetActiveCSNameVec(s_index)
        slat_added_set = {available_cs_vec[i-1] for i in slat_ind_vec}
        active_slat_set = set(active_slat_vec)
        if active_slat_set != slat_added_set:
            print("Incorrect Slat set")
        print("The slat group contains:")
        print(active_slat_set)

        vsp.Update()

        r_index = vsp.CreateVSPAEROControlSurfaceGroup()
        rudder_group = "Rudder_Group"
        vsp.SetVSPAEROControlGroupName( rudder_group, r_index)
        rudder_ind_vec = [11, 12]
        vsp.AddSelectedToCSGroup(rudder_ind_vec, r_index)
        active_rudder_vec = vsp.GetActiveCSNameVec(r_index)
        rudder_added_set = {available_cs_vec[i-1] for i in rudder_ind_vec}
        active_rudder_set = set(active_rudder_vec)
        if active_rudder_set != rudder_added_set:
            print("Incorrect rudder set")
        print("The rudder group contains:")
        print(active_rudder_set)
    

    #def deflections():
        #Control surface deflections

        cs_group_container_id = vsp.FindContainer("VSPAEROSettings", 0)

        vsp.Update()

        #Slat1 Deflection

        #vsp.SetParmVal(vsp.FindParm(cs_group_container_id, "Surf_" + slat1 + "_0_Gain", "ControlSurfaceGroup_1"), 1 )
        #vsp.SetParmVal(vsp.FindParm(cs_group_container_id, "Surf_" + slat1 + "_1_Gain", "ControlSurfaceGroup_1"), -1 )
        #vsp.SetParmVal(vsp.FindParm(cs_group_container_id, "DeflectionAngle", "ControlSurfaceGroup_1"), 10.0)
        #slat_deflection_angle = vsp.GetParmVal(cs_group_container_id, "DeflectionAngle", "ControlSurfaceGroup_1")
        #print(f"Slat Deflection Angle: {slat_deflection_angle} degrees")
    

    # Run the main function
    if __name__ == "__main__":
        build_aircraft(outerwing_span=outerwing_span, root_chord_array=root_chord_array, rib_pitch=rib_pitch, num_ribs=num_ribs)
    
    vsp.ClearVSPModel()
    vsp.Update()
    vsp.ReadVSPFile("Prelim_3.vsp3")
    vsp.Update()
    print("--> Computing Geometry")
    #print("")
    # Set up analysis for VSPAero
    analysis_name = "VSPAEROComputeGeometry"
    #print(analysis_name)

    # Set default inputs
    vsp.SetAnalysisInputDefaults(analysis_name)
    # Set to Vortex Lattice Method
    analysis_method = list(vsp.GetIntAnalysisInput(analysis_name, "AnalysisMethod"))
    analysis_method[0] = vsp.VORTEX_LATTICE
    vsp.SetIntAnalysisInput( analysis_name, "AnalysisMethod", analysis_method )  # Vortex Lattice is typically 1

    # List inputs
    #vsp.PrintAnalysisInputs(analysis_name)
    #print("")

    # Execute
    print("\tExecuting...")
    rid = vsp.ExecAnalysis(analysis_name)
    print("COMPLETE")

    # Get and display results
    #vsp.PrintResults(rid)
    #print("")


    print("--> Computing VSPAERO")
    #print("")

    # Set up VSPAero Sweep analysis
    analysis_name = "VSPAEROSweep"
    #print(analysis_name)

    # Set default inputs
    vsp.SetAnalysisInputDefaults(analysis_name)

    # Set to Vortex Lattice Method
    analysis_method = list(vsp.GetIntAnalysisInput(analysis_name, "AnalysisMethod"))
    analysis_method[0] = vsp.VORTEX_LATTICE
    vsp.SetIntAnalysisInput(analysis_name, "AnalysisMethod", analysis_method )  # Vortex Lattice is typically 1
    vsp.SetDoubleAnalysisInput(analysis_name, "MachStart", (0.3,), 0)
    vsp.SetDoubleAnalysisInput(analysis_name, "MachNpts", (1,), 0)
    vsp.SetDoubleAnalysisInput(analysis_name, "ReCref", (6e+07,), 0)
    vsp.SetDoubleAnalysisInput(analysis_name, "Sref", (296.529,), 0)
    vsp.SetDoubleAnalysisInput(analysis_name, "bref", (41.22,), 0)
    vsp.SetDoubleAnalysisInput(analysis_name, "cref", (4.975,), 0)
    vsp.SetIntAnalysisInput(analysis_name, "AlphaNpts", (3,), 0)
    vsp.SetDoubleAnalysisInput(analysis_name, "AlphaStart", (0,), 0)
    vsp.SetDoubleAnalysisInput(analysis_name, "AlphaEnd", (1,), 0)
    vsp.SetDoubleAnalysisInput(analysis_name, "Clmax" , (1.5,), 0)
    vsp.SetDoubleAnalysisInput(analysis_name, "Vinf", (98.0,), 0)
    vsp.SetIntAnalysisInput(analysis_name, "NCPU", [4], 0)
    vsp.SetIntAnalysisInput(analysis_name, "WakeNumIter", [3], 0)
    vsp.SetDoubleAnalysisInput(analysis_name, "Rho", [0.8], 0)
    # Set AlphaNpts (assuming 1 point)
    vsp.SetIntAnalysisInput(analysis_name, "AlphaNpts", (1,), 0)

    # List inputs
    vsp.PrintAnalysisInputs(analysis_name)
    print("")

    # Execute
    print("\tExecuting...")
    rid = vsp.ExecAnalysis(analysis_name)
    print("COMPLETE")

    # Get and display results
    #vsp.PrintResults(rid)
    #print("")


    print("--> Generating Cp Slices")
    #print("")

    # Set up Cp Slicer analysis
    analysis_name = "CpSlicer"
    #print(analysis_name)

    # Set default inputs
    vsp.SetAnalysisInputDefaults(analysis_name)

    # Set to Vortex Lattice Method
    analysis_method = list(vsp.GetIntAnalysisInput(analysis_name, "AnalysisMethod"))
    analysis_method[0] = vsp.VORTEX_LATTICE
    vsp.SetIntAnalysisInput( analysis_name, "AnalysisMethod", analysis_method )

    # Set up YSlice positions
    ycuts = [0.0, 1.5, 3.0, 4.5, 5.0, 8.0]
    vsp.SetDoubleAnalysisInput(analysis_name, "YSlicePosVec", ycuts)
    zcuts = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 3.9]
    vsp.SetDoubleAnalysisInput(analysis_name, "ZSlicePosVec", zcuts)

    # List inputs
    #vsp.PrintAnalysisInputs(analysis_name)
    #print("")

    # Execute
    print("\tExecuting...")
    rid = vsp.ExecAnalysis(analysis_name)
    print("COMPLETE")
    
build_initial_geometry()




        #cs_group_container_id = vsp.FindContainer("VSPAEROSettings", 0)

        #vsp.SetParmVal(vsp.FindParm(cs_group_container_id, "Surf_" + rudder_id + "_0_Gain", "ControlSurfaceGroup_2"), 1 )
        #vsp.SetParmVal(vsp.FindParm(cs_group_container_id, "Surf_" + rudder_id + "_1_Gain", "ControlSurfaceGroup_2"), -1 )
        #vsp.SetParmVal(vsp.FindParm(cs_group_container_id, "DeflectionAngle", "ControlSurfaceGroup_2"), 10.0)

