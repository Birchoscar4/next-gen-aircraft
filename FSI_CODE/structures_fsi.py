import openvsp as vsp
import csv

# Load an existing OpenVSP model
vsp.ReadVSPFile("Prelim_2.vsp3")  # Replace with your OpenVSP model file

# Get the Wing Geom ID (assuming it's a wing)
wing_id = vsp.FindGeomsWithName("WING")  # Replace "WING" with your actual component name
def structures_fsi(wing_id, start_index = 2, end_index = 10):
    wing_id = wing_id[0]

    num_xsec_surfs = vsp.GetNumXSecSurfs(wing_id)
    print(f"Number of XSecSurfs: {num_xsec_surfs}")

    # Use the main XSecSurf (index 0)
    xsec_surf = vsp.GetXSecSurf(wing_id, 0)

    # Get the number of cross-sections in the XSecSurf
    num_xsec = vsp.GetNumXSec(xsec_surf)
    print(f"Number of XSec in XSecSurf 0: {num_xsec}")
    
    twist_values = []
    delta_y_values = []

    try:
        with open("E:\Engineering\Aerospace Engineering 2021-2025\Aerospace 2021-2025\\4th Year\MECH 5080 Team Project\CODES\Working Codes\def_twist.csv", mode="r") as file:  # Replace with your .csv file
            reader = csv.DictReader(file)
            print("CSV Columns:", reader.fieldnames)  # Debug: Print column names
            for row in reader:
                print("Row:", row)  # Debug: Print each row
                try:
                    delta_y = float(row["Deflection "]) # Column name in .csv
                    twist = float(row["Twist "]) # Column name in .csv
                    delta_y_values.append(delta_y)
                    twist_values.append(twist)
                except KeyError as e:
                    print(f"Error: Missing column in CSV - {e}")
                except ValueError as e:
                    print(f"Error: Invalid value in CSV - {e}")
    except FileNotFoundError:
        print("Error: CSV file not found. Ensure the file path is correct.")
    
    for i in range(num_xsec):
        # Get the XSec ID of the section to modify
        xsec_id = vsp.GetXSec(xsec_surf, i)

        if xsec_id:
                # Update twist for XSec 3 onwards (index 2)
            if i >= 2:  # XSec indices start from 0, so XSec 3 is index 2
                twist_index = i - 2  # Map XSec index to CSV index
                if twist_index < len(twist_values):
                    twist_parm_id = vsp.GetParm(xsec_id, "Twist", "XSec")
                    if twist_parm_id:
                        print(f"Updating Twist for XSec {i+1} to {twist_values[twist_index]}.")
                        vsp.SetParmVal(twist_parm_id, twist_values[twist_index])
                    else:
                        print(f"Error: Could not find Twist parameter for XSec {i+1}.")
                        
        # Apply changes
        vsp.Update()


    for xsec_index in range(start_index, end_index + 1):  # xsec_index is an integer
        
        xsec_id = vsp.GetXSec(xsec_surf, xsec_index)

        param_name = f"XSecCurve_{xsec_index}"  # Construct the parameter name as a string
        
        # Get the parm ID for the Deflection of the wing
        param_id = vsp.GetParm(wing_id, "DeltaY", param_name)
        
        # Define variable to get current value of parameter of interest
        deflection = vsp.GetParmVal(param_id)

        # Utilize parmID to update the value of the parameter of interest
        
        deflection_new = vsp.SetParmVal(param_id, delta_y_values[xsec_index])
        print(f"Updating Deflection for XSecCurve {xsec_index} to {delta_y_values[xsec_index]}.")

    return delta_y_values

structures_fsi(wing_id)
vsp.WriteVSPFile("structures_fsi_model.vsp3")
print("Twist and DeltaY updated successfully.")
print("Model saved correctly as structures_fsi_model.vsp3")

vsp.ClearVSPModel()
vsp.Update()
vsp.ReadVSPFile("structures_fsi_model.vsp3")
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
