def aero_analysis():
    import openvsp as vsp
    import os
    import subprocess
    import platform
    print ("Beginning CP analysis")
    stdout = vsp.cvar.cstdout
    errorMgr = vsp.ErrorMgrSingleton.getInstance()

    vsp.VSPCheckSetup()
    vsp.VSPRenew()
    vsp.ClearVSPModel()
    vsp.DeleteAllResults()
    # read the file again 
    vsp.ReadVSPFile("jetzero_edited_airfoil.vsp3")

    # Export the geometry that used for anaylsis
    compGeom = "VSPAEROComputeGeometry"
    vsp.SetAnalysisInputDefaults(compGeom)
    compGeom_results = vsp.ExecAnalysis(compGeom)

    VLMAnalysis = "VSPAEROSweep"

    vsp.SetAnalysisInputDefaults(VLMAnalysis)

    # Set the analysis type to vlm
    analysis_method = [vsp.VORTEX_LATTICE]
    vsp.SetIntAnalysisInput(VLMAnalysis, "AnalysisMethod", analysis_method)
    vsp.SetDoubleAnalysisInput(VLMAnalysis, "MachStart", (0.8,), 0)
    vsp.SetDoubleAnalysisInput(VLMAnalysis, "MachNpts", (1,), 0)
    vsp.SetDoubleAnalysisInput(VLMAnalysis, "ReCref", (6e+07,), 0)
    vsp.SetDoubleAnalysisInput(VLMAnalysis, "Sref", (339.545,), 0)
    vsp.SetDoubleAnalysisInput(VLMAnalysis, "bref", (46.228,), 0)
    vsp.SetDoubleAnalysisInput(VLMAnalysis, "cref", (9.399,), 0)
    vsp.SetDoubleAnalysisInput(VLMAnalysis, "AlphaStart", (0,), 0)
    vsp.SetDoubleAnalysisInput(VLMAnalysis, "AlphaEnd", (1,), 0)
    vsp.SetIntAnalysisInput(VLMAnalysis, "AlphaNpts", (1,), 0)
    vsp.SetDoubleAnalysisInput(VLMAnalysis, "ClMax" , (1.5,), 0)
    vsp.SetIntAnalysisInput(VLMAnalysis, "NCPU", [4], 0)
    vsp.SetIntAnalysisInput(VLMAnalysis, "WakeNumIter", [5], 0)
    analysis_cp = "CpSlicer"
    print(analysis_cp)

    # Set default inputs
    vsp.SetAnalysisInputDefaults(analysis_cp)

    # Set to Vortex Lattice Method
    analysis_method = list(vsp.GetIntAnalysisInput(analysis_cp, "AnalysisMethod"))
    analysis_method[0] = vsp.VORTEX_LATTICE
    vsp.SetIntAnalysisInput( analysis_cp, "AnalysisMethod", analysis_method )
    # Set up YSlice positions
    ycuts = [2.0, 4.5, 8.0]
    vsp.SetDoubleAnalysisInput(analysis_cp, "YSlicePosVec", ycuts)

    # List inputs
    vsp.PrintAnalysisInputs(analysis_cp)
    print("")
    vsp.PrintAnalysisInputs(VLMAnalysis)

    vsp.Update()

    #aileron control group
    a_index = vsp.CreateVSPAEROControlSurfaceGroup()
    aileron_group = "Aileron_Group"
    vsp.SetVSPAEROControlGroupName( aileron_group, a_index)
    available_cs_vec = vsp.GetAvailableCSNameVec( a_index)
    print (available_cs_vec)
    a_ind_vec = [1, 2]
    vsp.AddSelectedToCSGroup(a_ind_vec, a_index)
    active_aileron_vec = vsp.GetActiveCSNameVec(a_index)
    aileron_added_set = {available_cs_vec[i-1] for i in a_ind_vec}
    active_aileron_set = set(active_aileron_vec)
    if active_aileron_set != aileron_added_set:
        print("Incorrect Aileron set")
    print("The aileron group contains:")
    print(active_aileron_set)

    #Flap control group
    f_index = vsp.CreateVSPAEROControlSurfaceGroup()
    flap_group = "Flap_Group"
    vsp.SetVSPAEROControlGroupName(flap_group, f_index)
    f_ind_vec = [3, 4]
    vsp.AddSelectedToCSGroup(f_ind_vec, f_index)
    active_flap_vec = vsp.GetActiveCSNameVec(f_index)
    flap_added_set = {available_cs_vec[i-1] for i in f_ind_vec}
    active_flap_set = set(active_flap_vec)
    if active_flap_set != flap_added_set:
        print("Incorrect")
    print("The flap group contains:")
    print(active_flap_set)

    #Rudder control group
    r_index = vsp.CreateVSPAEROControlSurfaceGroup()
    rudder_group = "Rudder_Group"
    vsp.SetVSPAEROControlGroupName(rudder_group, r_index)
    r_ind_vec = [5, 6]
    vsp.AddSelectedToCSGroup(r_ind_vec, r_index)
    active_rudder_vec = vsp.GetActiveCSNameVec(r_index)
    rudder_added_set = {available_cs_vec[i-1] for i in r_ind_vec}
    active_rudder_set = set(active_rudder_vec)
    if active_rudder_set != rudder_added_set:
        print("Incorrect")
    print("The rudder group contains:")
    print(active_rudder_set)

    vsp.Update()

    #Control surface deflections

    cs_group_container_id = vsp.FindContainer("VSPAEROSettings", 0)

    vsp.Update()

    #Aileron Deflection

    vsp.SetParmVal(vsp.FindParm(cs_group_container_id, "Surf_HWIDFJVLZN_0_Gain", "ControlSurfaceGroup_0"), 1.0 )
    vsp.SetParmVal(vsp.FindParm(cs_group_container_id, "Surf_HWIDFJVLZN_1_Gain", "ControlSurfaceGroup_0"), 1.0 )
    vsp.SetParmVal(vsp.FindParm(cs_group_container_id, "DeflectionAngle", "ControlSurfaceGroup_0"), 0.0)
    aileron_deflection_angle = vsp.GetParmVal(cs_group_container_id, "DeflectionAngle", "ControlSurfaceGroup_0")
    print(f"Aileron Deflection Angle: {aileron_deflection_angle} degrees")

    #Flap Deflection

    vsp.SetParmVal(vsp.FindParm(cs_group_container_id, "Surf_ZXGFZXDAMI_0_Gain", "ControlSurfaceGroup_1"), 1.0 )
    vsp.SetParmVal(vsp.FindParm(cs_group_container_id, "Surf_ZXGFZXDAMI_1_Gain", "ControlSurfaceGroup_1"), -1.0 )
    vsp.SetParmVal(vsp.FindParm(cs_group_container_id, "DeflectionAngle", "ControlSurfaceGroup_1"), 15.0)
    flap_deflection_angle = vsp.GetParmVal(cs_group_container_id, "DeflectionAngle", "ControlSurfaceGroup_1")
    print(f"Flap Deflection Angle: {flap_deflection_angle} degrees")

    #Rudder Deflection

    vsp.SetParmVal(vsp.FindParm(cs_group_container_id, "Surf_GERLJOVLFG_0_Gain", "ControlSurfaceGroup_2"), 1.0 )
    vsp.SetParmVal(vsp.FindParm(cs_group_container_id, "Surf_GERLJOVLFG_1_Gain", "ControlSurfaceGroup_2"), -1.0 )
    vsp.SetParmVal(vsp.FindParm(cs_group_container_id, "DeflectionAngle", "ControlSurfaceGroup_2"), 0.0)
    rudder_deflection_angle = vsp.GetParmVal(cs_group_container_id, "DeflectionAngle", "ControlSurfaceGroup_2")
    print(f"Flap Deflection Angle: {rudder_deflection_angle} degrees")
    allResults = vsp.ExecAnalysis(VLMAnalysis)

    #vsp.WriteResultsCSVFile(allResults, "Results.csv")
    #print ("Aero analysis has been saved as Results.csv")

    #Save the file 
    vsp.WriteVSPFile("jetzero_edited_airfoil.vsp3" , vsp.SET_ALL)

    vsp.Update()
        
    # clear all the data in the code
    vsp.ClearVSPModel()

    print ("Beginning Stability analysis")
    stdout = vsp.cvar.cstdout
    errorMgr = vsp.ErrorMgrSingleton.getInstance()

    vsp.VSPCheckSetup()
    vsp.VSPRenew()
    vsp.ClearVSPModel()
    vsp.DeleteAllResults()

    # read the file again 
    vsp.ReadVSPFile("jetzero_edited_airfoil.vsp3")

    # Export the geometry that used for anaylsis
    compGeom = "VSPAEROComputeGeometry"
    vsp.SetAnalysisInputDefaults(compGeom)
    compGeom_results = vsp.ExecAnalysis(compGeom)

    VLMAnalysis = "VSPAEROSweep"

    vsp.SetAnalysisInputDefaults(VLMAnalysis)

    # Set the analysis type to vlm
    analysis_method = [vsp.VORTEX_LATTICE]
    vsp.SetIntAnalysisInput(VLMAnalysis, "AnalysisMethod", analysis_method)
    vsp.SetDoubleAnalysisInput(VLMAnalysis, "Mach", (0.8,), 0)
    vsp.SetDoubleAnalysisInput(VLMAnalysis, "ReCref", (6e+07,), 0)
    vsp.SetDoubleAnalysisInput(VLMAnalysis, "Sref", (339.545,), 0)
    vsp.SetDoubleAnalysisInput(VLMAnalysis, "bref", (46.228,), 0)
    vsp.SetDoubleAnalysisInput(VLMAnalysis, "cref", (9.399,), 0)
    vsp.SetDoubleAnalysisInput(VLMAnalysis, "AlphaStart", (0,), 0)
    vsp.SetDoubleAnalysisInput(VLMAnalysis, "AlphaEnd", (1,), 0)
    vsp.SetIntAnalysisInput(VLMAnalysis, "AlphaNpts", (1,), 0)
    vsp.SetDoubleAnalysisInput(VLMAnalysis, "ClMax" , (1.5,), 0)
    vsp.SetDoubleAnalysisInput(VLMAnalysis, "Vinf", (98.0,), 0)
    vsp.SetIntAnalysisInput(VLMAnalysis, "NCPU", [4], 0)
    vsp.SetIntAnalysisInput(VLMAnalysis, "WakeNumIter", [5], 0)
    stability_type = [vsp.STABILITY_DEFAULT]
    vsp.SetIntAnalysisInput(VLMAnalysis, "UnsteadyType", stability_type)

    #aileron control group
    a_index = vsp.CreateVSPAEROControlSurfaceGroup()
    aileron_group = "Aileron_Group"
    vsp.SetVSPAEROControlGroupName( aileron_group, a_index)
    available_cs_vec = vsp.GetAvailableCSNameVec( a_index)
    print (available_cs_vec)
    a_ind_vec = [1, 2]
    vsp.AddSelectedToCSGroup(a_ind_vec, a_index)
    active_aileron_vec = vsp.GetActiveCSNameVec(a_index)
    aileron_added_set = {available_cs_vec[i-1] for i in a_ind_vec}
    active_aileron_set = set(active_aileron_vec)
    if active_aileron_set != aileron_added_set:
        print("Incorrect Aileron set")
    print("The aileron group contains:")
    print(active_aileron_set)

    #Flap control group
    f_index = vsp.CreateVSPAEROControlSurfaceGroup()
    flap_group = "Flap_Group"
    vsp.SetVSPAEROControlGroupName(flap_group, f_index)
    f_ind_vec = [3, 4]
    vsp.AddSelectedToCSGroup(f_ind_vec, f_index)
    active_flap_vec = vsp.GetActiveCSNameVec(f_index)
    flap_added_set = {available_cs_vec[i-1] for i in f_ind_vec}
    active_flap_set = set(active_flap_vec)
    if active_flap_set != flap_added_set:
        print("Incorrect")
    print("The flap group contains:")
    print(active_flap_set)

    #Rudder control group
    r_index = vsp.CreateVSPAEROControlSurfaceGroup()
    rudder_group = "Rudder_Group"
    vsp.SetVSPAEROControlGroupName(rudder_group, r_index)
    r_ind_vec = [5, 6]
    vsp.AddSelectedToCSGroup(r_ind_vec, r_index)
    active_rudder_vec = vsp.GetActiveCSNameVec(r_index)
    rudder_added_set = {available_cs_vec[i-1] for i in r_ind_vec}
    active_rudder_set = set(active_rudder_vec)
    if active_rudder_set != rudder_added_set:
        print("Incorrect")
    print("The rudder group contains:")
    print(active_rudder_set)

    vsp.Update()

    #Control surface deflections

    cs_group_container_id = vsp.FindContainer("VSPAEROSettings", 0)

    vsp.Update()

    #Aileron Deflection

    vsp.SetParmVal(vsp.FindParm(cs_group_container_id, "Surf_HWIDFJVLZN_0_Gain", "ControlSurfaceGroup_0"), 1.0 )
    vsp.SetParmVal(vsp.FindParm(cs_group_container_id, "Surf_HWIDFJVLZN_1_Gain", "ControlSurfaceGroup_0"), 1.0 )
    vsp.SetParmVal(vsp.FindParm(cs_group_container_id, "DeflectionAngle", "ControlSurfaceGroup_0"), 0.0)
    aileron_deflection_angle = vsp.GetParmVal(cs_group_container_id, "DeflectionAngle", "ControlSurfaceGroup_0")
    print(f"Aileron Deflection Angle: {aileron_deflection_angle} degrees")

    #Flap Deflection

    vsp.SetParmVal(vsp.FindParm(cs_group_container_id, "Surf_ZXGFZXDAMI_0_Gain", "ControlSurfaceGroup_1"), 1.0 )
    vsp.SetParmVal(vsp.FindParm(cs_group_container_id, "Surf_ZXGFZXDAMI_1_Gain", "ControlSurfaceGroup_1"), -1.0 )
    vsp.SetParmVal(vsp.FindParm(cs_group_container_id, "DeflectionAngle", "ControlSurfaceGroup_1"), 15.0)
    flap_deflection_angle = vsp.GetParmVal(cs_group_container_id, "DeflectionAngle", "ControlSurfaceGroup_1")
    print(f"Flap Deflection Angle: {flap_deflection_angle} degrees")

    #Rudder Deflection

    vsp.SetParmVal(vsp.FindParm(cs_group_container_id, "Surf_GERLJOVLFG_0_Gain", "ControlSurfaceGroup_2"), 1.0 )
    vsp.SetParmVal(vsp.FindParm(cs_group_container_id, "Surf_GERLJOVLFG_1_Gain", "ControlSurfaceGroup_2"), -1.0 )
    vsp.SetParmVal(vsp.FindParm(cs_group_container_id, "DeflectionAngle", "ControlSurfaceGroup_2"), 0.0)
    rudder_deflection_angle = vsp.GetParmVal(cs_group_container_id, "DeflectionAngle", "ControlSurfaceGroup_2")
    print(f"Flap Deflection Angle: {rudder_deflection_angle} degrees")
    allResults = vsp.ExecAnalysis(VLMAnalysis)

    #vsp.WriteResultsCSVFile(allResults, "Results.csv")
    #print ("Aero analysis has been saved as Results.csv")

    #Save the file 
    vsp.WriteVSPFile("jetzero_edited_airfoil.vsp3" , vsp.SET_ALL)

    vsp.Update()

    vsp.ClearVSPModel()

aero_analysis()
