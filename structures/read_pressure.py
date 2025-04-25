import numpy as np
import re
import pandas as pd
import matplotlib.pyplot as plt

file_path="pressure.slc"

# Function to compute centrol of pressure location at each slice
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
            x_CoP_percentage = (x_CoP - x_min) / chord_length * 100  # Convert to percentage
        else:
            x_CoP_percentage = np.nan  # Handle cases with zero sum

        results.append((y_slice, x_CoP_percentage))

    return results

# Function to compute sectional lift coefficient Cl at each slice
def compute_sectional_lift_coefficient(file_path):
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

    # Compute Cl for each slice
    results = []
    for y_slice, data in slices.items():
        if len(data) < 2:
            continue  # Need at least two points to compute lift

        # Sort data by x (in case they are unordered)
        data.sort(key=lambda point: point[0])
        x_vals, dCp_vals = np.array(data)[:, 0], np.array(data)[:, 1]

        # Compute delta x (spacing between adjacent points)
        dx_vals = np.diff(x_vals)

        # Compute chord length (max difference in x)
        chord_length = np.max(x_vals) - np.min(x_vals)
        if chord_length == 0:
            continue  # Avoid division by zero

        # Compute sectional lift coefficient Cl
        Cl = np.sum(-dCp_vals[:-1] * dx_vals) / chord_length

        results.append((y_slice, Cl))

    return results

# Function to compute sectional lift forces at each slice
def compute_total_lift_force(file_path, q_inf):
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

        lift_slices.append((y_slice, chord_length, Cl, L_prime))

    # Convert to numpy array for integration
    lift_slices = np.array(lift_slices)

    # Compute total lift force using trapezoidal integration
    if len(lift_slices) > 1:
        y_positions = lift_slices[:, 0]  # Spanwise positions
        L_prime_values = lift_slices[:, 3]  # Lift per unit span
        total_lift = np.trapz(L_prime_values, y_positions)  # Integrate L' over span
    else:
        total_lift = 0

    return lift_slices, total_lift

# Function to extract root chord, taper ratio, and sweep angle
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

    # Compute leading edge sweep angle
    delta_x_LE = x_LE_tip - x_LE_root
    delta_y = y_tip - y_root
    sweep_LE = np.degrees(np.arctan(delta_x_LE / delta_y)) if delta_y != 0 else np.nan

    return {
        "Root Chord (c_root)": c_root,
        "Tip Chord (c_tip)": c_tip,
        "Taper Ratio (λ)": taper_ratio,
        "Leading Edge Sweep Angle (Λ_LE)": sweep_LE
    }

wing_geometry = extract_wing_geometry(file_path)

cop_results_updated = compute_x_cop_percentage(file_path)
df_cop_updated = pd.DataFrame(cop_results_updated, columns=["Slice_Y", "X_CoP"])

cl_results = compute_sectional_lift_coefficient(file_path)
df_cl = pd.DataFrame(cl_results, columns=["Slice_Y", "C_L"])


q_inf_example = 500.0  # Dynamic pressure in N/m²
lift_slices_data, total_lift_force = compute_total_lift_force(file_path, q_inf_example)
df_lift_slices = pd.DataFrame(lift_slices_data, columns=["Slice_Y", "Chord Length", "C_L", "Lift Per Unit Span"])

# Display the results
df_cop_updated.plot(x="Slice_Y", y="X_CoP", kind="line", marker="o", title="CoP Plot")
df_cl.plot(x="Slice_Y", y="C_L", kind="line", marker="o", title="CL Plot")
df_lift_slices.plot(x="Slice_Y", y="Lift Per Unit Span", kind="line", marker="o", title="Lift Per Unit Span")
plt.show()
