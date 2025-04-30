def static_stability_vlm():
    import openvsp as vsp
    import numpy as np
    vsp.ClearVSPModel()
    vsp.Update()
    vsp.ReadVSPFile("update_geom_pos.vsp3")
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

    print("--> Generating Static Stability Analysis")
    #print("")

    # Set up Cp Slicer analysis
    analysis_name = "VSPAEROSweep"
    #print(analysis_name)

    # Set default inputs
    vsp.SetAnalysisInputDefaults(analysis_name)

    # Set to Vortex Lattice Method
    analysis_method = list(vsp.GetIntAnalysisInput(analysis_name, "AnalysisMethod"))
    analysis_method[0] = vsp.VORTEX_LATTICE
    vsp.SetIntAnalysisInput( analysis_name, "AnalysisMethod", analysis_method )
    stability_type = [vsp.STABILITY_DEFAULT]
    vsp.SetIntAnalysisInput(analysis_name, "UnsteadyType", stability_type)
    vsp.SetDoubleAnalysisInput(analysis_name, "MachStart", (0.3,), 0)
    vsp.SetDoubleAnalysisInput(analysis_name, "MachNpts", (1,), 0)
    vsp.SetDoubleAnalysisInput(analysis_name, "ReCref", (6e+07,), 0)
    vsp.SetDoubleAnalysisInput(analysis_name, "Sref", (296.529,), 0)
    vsp.SetDoubleAnalysisInput(analysis_name, "bref", (41.22,), 0)
    vsp.SetDoubleAnalysisInput(analysis_name, "cref", (4.975,), 0)
    vsp.SetIntAnalysisInput(analysis_name, "AlphaNpts", (1,), 0)
    vsp.SetDoubleAnalysisInput(analysis_name, "AlphaStart", (0,), 0)
    vsp.SetDoubleAnalysisInput(analysis_name, "AlphaEnd", (1,), 0)
    vsp.SetDoubleAnalysisInput(analysis_name, "Clmax" , (1.5,), 0)
    vsp.SetDoubleAnalysisInput(analysis_name, "Vinf", (98.0,), 0)
    vsp.SetIntAnalysisInput(analysis_name, "NCPU", [4], 0)
    vsp.SetIntAnalysisInput(analysis_name, "WakeNumIter", [3], 0)
    vsp.SetDoubleAnalysisInput(analysis_name, "Rho", [0.8], 0)

    # List inputs
    vsp.PrintAnalysisInputs(analysis_name)
    print("")

    # Execute
    print("\tExecuting...")
    rid = vsp.ExecAnalysis(analysis_name)
    print("COMPLETE")
    
    print("--> Generating QSTAB Analysis")
    #print("")

    # Set up Cp Slicer analysis
    analysis_name = "VSPAEROSweep"
    #print(analysis_name)

    # Set default inputs
    vsp.SetAnalysisInputDefaults(analysis_name)

    # Set to Vortex Lattice Method
    analysis_method = list(vsp.GetIntAnalysisInput(analysis_name, "AnalysisMethod"))
    analysis_method[0] = vsp.VORTEX_LATTICE
    vsp.SetIntAnalysisInput( analysis_name, "AnalysisMethod", analysis_method )
    stability_type = [vsp.STABILITY_Q_ANALYSIS]
    vsp.SetIntAnalysisInput(analysis_name, "UnsteadyType", stability_type)
    # List inputs
    #vsp.PrintAnalysisInputs(analysis_name)
    #print("")

    # Execute
    print("\tExecuting...")
    rid = vsp.ExecAnalysis(analysis_name)
    print("COMPLETE") 