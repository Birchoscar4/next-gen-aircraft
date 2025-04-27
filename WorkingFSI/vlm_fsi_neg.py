def vlm_fsi_neg():  
    import openvsp as vsp
    import os
    import subprocess
    import platform

    vsp.ClearVSPModel()
    vsp.Update()
    vsp.ReadVSPFile("update_geom_neg.vsp3")
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
    print("\tExecuting Negative Defelction Analysis")
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
    ycuts = [8.5, 9.7782, 11.0564, 12.3346, 13.6128, 14.891, 16.169201, 17.447399, 18.725599]
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