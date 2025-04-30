import pandas as pd
import matplotlib.pyplot as plt
def vlm_fsi_pos():  
    import openvsp as vsp
    import os
    import subprocess
    import platform
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
    print("\tExecuting Positive Deflection Analysis")
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
    vsp.SetIntAnalysisInput(analysis_name, "AlphaNpts", (1,), 0)
    vsp.SetDoubleAnalysisInput(analysis_name, "AlphaStart", (1,), 0)
    vsp.SetDoubleAnalysisInput(analysis_name, "AlphaEnd", (10,), 0)
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
    ycuts = np.linspace(8.5, 18.72, 12)
    vsp.SetDoubleAnalysisInput(analysis_name, "YSlicePosVec", ycuts)
    #zcuts = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 3.9]
    #vsp.SetDoubleAnalysisInput(analysis_name, "ZSlicePosVec", zcuts)

    # List inputs
    #vsp.PrintAnalysisInputs(analysis_name)
    #print("")

    # Execute
    print("\tExecuting...")
    rid = vsp.ExecAnalysis(analysis_name)
    print("COMPLETE")
if __name__ == "__main__":
    vlm_fsi_pos()
"""
import pandas as pd

df = pd.read_csv('C:/Users/joepe/Documents/ProjectAero/OpenVSP-3.41.1-win64-Python3.11/OpenVSP-3.41.1-win64/python/update_geom_pos_DegenGeom.polar', delim_whitespace=True)

print(df.head())
print(df.shape)  # prints (rows, columns)

# To print all the data (careful if it's large!)
#print(df.to_string())
#df['CL'] = df['CL'] * (2 * df['AoA'])
df['CL'] = df['CL'] + 0.27483
print(df.head()) 
# Find the index where CL is minimum
closest_to_zero_index = (df['CL'].abs()).idxmin()

# Step 3: Get the 'cdtot' corresponding to that 'cl' value
cdtot_at_min_cl = df.loc[closest_to_zero_index, 'CDtot']

# Subtract that value from ALL CDtot entries
df['CDtot'] = (df['CDtot'] - cdtot_at_min_cl)+0.02018

print(df.head())  # just to check

df.loc[df['AoA'] < 0, 'CDtot'] += 0.016

df['AoA_Squared'] = df['AoA'] ** 2

df.loc[df['AoA'] > 0, 'CDtot'] = (df.loc[df['AoA'] > 0, 'AoA_Squared']*0.0018) + (df.loc[df['AoA'] > 0, 'AoA']*0.004) + df.loc[df['AoA'] > 0, 'CDtot'] + 0.02
df.loc[df['AoA'] == 0, 'CDtot'] = (df.loc[df['AoA'] == 0, 'AoA_Squared']*0.0018) + (df.loc[df['AoA'] == 0, 'AoA']*0.004) + df.loc[df['AoA'] == 0, 'CDtot'] + 0.02
print(df.head()) 

plt.figure(figsize=(8, 6))
plt.scatter(df['AoA'], df['CDtot'], color='blue', alpha=0.5)
plt.title("Plot of AoA vs CDtot")
plt.xlabel("AoA")
plt.ylabel("CDtot")
plt.grid(True)
plt.show()

plt.figure(figsize=(8, 6))
plt.scatter(df['AoA'], df['CL'], color='blue', alpha=0.5)
plt.title("Plot of AoA vs CL")
plt.xlabel("AoA")
plt.ylabel("CL")
plt.grid(True)
plt.show() """