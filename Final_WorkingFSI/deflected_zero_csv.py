import openvsp as vsp
import csv
def neg_geom_zero(wing_id, start_index = 2, end_index = 10):
    # Load an existing OpenVSP model
    #vsp.ReadVSPFile("Prelim_2.vsp3")  # Replace with your OpenVSP model file

    # Get the Wing Geom ID (assuming it's a wing)
    wing_id = vsp.FindGeomsWithName("WING")  # Replace "WING" with your actual component name
    wing_id = wing_id[0]

    num_xsec_surfs = vsp.GetNumXSecSurfs(wing_id)
    #print(f"Number of XSecSurfs: {num_xsec_surfs}")

    # Use the main XSecSurf (index 0)
    xsec_surf = vsp.GetXSecSurf(wing_id, 0)

    # Get the number of cross-sections in the XSecSurf
    num_xsec = vsp.GetNumXSec(xsec_surf)
    #print(f"Number of XSec in XSecSurf 0: {num_xsec}")
    
    twist_values = []
    delta_y_values = []

    # UPLOAD THE DIRECTORY OF THE NEGATIVE LOAD .CSV
    try:
        with open(r"C:\Users\joepe\Documents\ProjectAero\OpenVSP-3.41.1-win64-Python3.11\OpenVSP-3.41.1-win64\python\def_twist_neg_zero.csv", mode="r") as file:  # Replace with your .csv file
            reader = csv.DictReader(file)
            #print("CSV Columns:", reader.fieldnames)  # Debug: Print column names
            for row in reader:
                #print("Row:", row)  # Debug: Print each row
                try:
                    delta_y = float(row["Deflection"]) # Column name in .csv
                    twist = float(row["Twist"]) # Column name in .csv
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
        vsp.Update()
    return delta_y_values