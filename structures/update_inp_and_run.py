import math
import numpy as np
import os
import time
print("Time of start:", time.strftime("%Y-%m-%d %H:%M:%S"))

def equivalent_plate_properties(fw, tw, tb, ts, d, E, w, rho, v):
    A_stringer = (d - ts) * tw + 2 * tb * fw
    A_plate = ts * w

    term1 = ((d - ts - 2 * tb) ** 3) * tw / 12
    term2 = 2 * ((tb ** 3 * fw) / 12 + tb * fw * ((d - ts - tb) ** 2) / 4)
    Ix_stringer = term1 + term2

    Ix_plate = w * ts ** 3 / 12

    centroid = ((ts / 2) * A_plate + ((d - ts) / 2 + ts) * A_stringer) / (A_stringer + A_plate)

    dy_stringer = centroid - ((d - ts) / 2 + ts)
    dy_plate = centroid - ts / 2

    I = Ix_plate + A_plate * dy_plate ** 2 + Ix_stringer + A_stringer * dy_stringer ** 2

    te = math.sqrt(12 * I / (A_stringer + A_plate))
    E1e = E * (A_stringer + A_plate) / (w * te)
    rhoe = rho * (A_stringer + A_plate) / (w * te)
    E2e = E * ts / te
    G12 = E / (2 * (1 + v))
    G12e = G12 * ts / te
    return te, E1e, rhoe, E2e, G12e

def update_abaqus_sections_and_material(file_path, root_moment, new_thickness, new_radius,
                                        new_density, new_elastic_lines,
                                        original_thickness=0.052, original_radius=0.02,
                                        original_density=264.):
    root_moment_pos = root_moment * 3.58
    root_moment_neg = root_moment * -1.58

    with open(file_path, 'r') as file:
        lines = file.readlines()

    updated_lines = []
    current_step = None
    i = 0

    while i < len(lines):
        line = lines[i]
        stripped = line.strip()

        # Track step
        if stripped.startswith("** STEP:"):
            if "pos" in stripped:
                current_step = "pos"
            elif "neg" in stripped:
                current_step = "neg"

        # Root moment logic
        if stripped.startswith("** Name: root_moment") and "Type: Moment" in stripped:
            updated_lines.append(line)
            i += 1
            line = lines[i].strip()
            if line.startswith("*Cload"):
                updated_lines.append(lines[i])
                i += 1
                if current_step == "pos":
                    updated_lines.append(f"Set-33, 6, {root_moment_pos:.6e}\n")
                elif current_step == "neg":
                    updated_lines.append(f"Set-33, 6, {root_moment_neg:.6e}\n")
                else:
                    updated_lines.append(lines[i])
                i += 1
                continue

        # Section, material, and elastic logic copied from your working code:
        elif stripped.startswith("** Section: orthoskin"):
            updated_lines.append(line)
            i += 1
            updated_lines.append(lines[i])
            i += 1
            thickness_line = lines[i].strip()
            if thickness_line.startswith(str(original_thickness)):
                parts = thickness_line.split(',')
                parts[0] = str(new_thickness)
                updated_lines.append(', '.join(parts) + '\n')
            else:
                updated_lines.append(lines[i])
            i += 1
            continue

        elif stripped.startswith("** Section: sparcap"):
            updated_lines.append(line)
            i += 1
            updated_lines.append(lines[i])
            i += 1
            radius_line = lines[i].strip()
            if radius_line == str(original_radius):
                updated_lines.append(f"{new_radius}\n")
            else:
                updated_lines.append(lines[i])
            i += 1
            continue

        elif stripped == "*Material, name=orthoaluminium":
            updated_lines.append(line)
            i += 1
            while i < len(lines):
                line = lines[i].strip()
                if line == "*Density":
                    updated_lines.append(lines[i])
                    i += 1
                    updated_lines.append(f"{new_density},\n")
                elif line.startswith("*Elastic"):
                    updated_lines.append(lines[i])  # keep the *Elastic line
                    i += 1
                    updated_lines.extend([l + '\n' for l in new_elastic_lines])
                    i += 2  # skip both original elastic lines
                    break

                else:
                    updated_lines.append(lines[i])
                i += 1
            continue

        updated_lines.append(line)
        i += 1

    new_file_path = file_path.replace(".inp", "_updated.inp")
    with open(new_file_path, 'w') as file:
        file.writelines(updated_lines)

    print("Successfully updated root moments and material sections.")
    print(f"Updated file saved as: {new_file_path}")


def format_elastic_values(E1e, E2e, G12e, v):
    new_elastic_values = [
        f"{E1e:.2e}, {E2e:.2e}, {E2e:.2e}, {v:.1f}, {v:.1f}, {v:.1f}, {G12e:.2e}, {G12e:.2e}",
        f"{G12e:.2e},"
    ]
    return new_elastic_values

fw = 0.03 #flange width
tw = 0.006 #stiffener thickness
tb = tw
ts = 0.004 #skin thickness
d = 0.15 #stiffener depth
E = 68.947573e9
w = 0.1 # stiffener spacing
rho = 1550
v = 0.4
sparcap_radius = 0.06

root_moment = 1160941

for root_moment in np.linspace(580470,1741410,3):
    for sparcap_radius in np.linspace(0.035,0.105,3):
        for w in np.linspace(0.075,0.225,3):
            te, E1e, rhoe, E2e, G12e = equivalent_plate_properties(fw,tw,tb,ts,d,E,w,rho,v)
            print(te, E1e, rhoe, E2e, G12e)

            formatted = format_elastic_values(E1e, E2e, G12e, v)
            print(formatted)

            file_path = "130425_final.inp"
            new_thickness = te
            new_radius = sparcap_radius
            new_density = rhoe
            new_elastic_values = formatted
            update_abaqus_sections_and_material(file_path, root_moment, new_thickness, new_radius, new_density, new_elastic_values, original_thickness=0.052, original_radius=0.02, original_density=264.)

            
            os.system("echo y | abaqus job=130425_final input=130425_final_updated.inp") #or datacheck
            
            

            job_name = "130425_final"
            
            time.sleep(20)
            # Wait until the .lck file is gone
            while os.path.exists(f"{job_name}.lck"):
                print("Job running...")
                time.sleep(5)  # check every 5 seconds
            

            print("Job complete!")
            print("Current time:", time.strftime("%Y-%m-%d %H:%M:%S"))
   
            import re
            import csv

            node_id = 90

            # Run Abaqus script
            os.system(f"abaqus python read_odb.py")

            # Read the CSV results
            results = {}
            with open("results.csv", "r") as f:
                reader = csv.reader(f)
                next(reader)  # Skip header
                for row in reader:
                    if len(row) < 2:
                        continue  # Skip empty or malformed lines
                    results[row[0]] = float(row[1])

            def extract_summary_from_csv(csv_file):
                target_vm_node = 90
                disp_node_1 = 16140
                disp_node_2 = 27377

                data = {}

                with open(csv_file, "r") as f:
                    reader = csv.DictReader(f)
                    for row in reader:
                        node = int(row["node"])
                        data[node] = {
                            "vm_max": float(row["vm_max"]),
                            "u1": float(row["u1"]),
                            "u2": float(row["u2"]),
                            "u3": float(row["u3"])
                        }

                # Extract values
                max_vm_stress = data[target_vm_node]["vm_max"]
                max_disp = data[disp_node_1]["u2"] - data[disp_node_2]["u2"]

                return max_vm_stress, max_disp

            # Example usage
            max_vm_stress, max_disp = extract_summary_from_csv("results.csv")


            # Use the values
            print("Max displacement:", max_disp)
            print("Max von Mises stress:", max_vm_stress)
            
            

            def get_mass_properties_from_dat(datfile, ref_point=(0.0, 0.0, 0.0), csv_out=None):
                def parse_block(text):
                    mass_match = re.search(r"TOTAL MASS OF MODEL\s+([\d.+-Ee]+)", text)
                    cg_match = re.search(r"CENTER OF MASS OF THE MODEL\s+([\d.+-Ee]+)\s+([\d.+-Ee]+)\s+([\d.+-Ee]+)", text)
                    moi_com_match = re.search(r"MOMENTS OF INERTIA ABOUT THE CENTER OF MASS\s+I\(XX\)\s+I\(YY\)\s+I\(ZZ\)\s+([\d.+-Ee]+)\s+([\d.+-Ee]+)\s+([\d.+-Ee]+)", text)
                    prod_com_match = re.search(r"PRODUCTS OF INERTIA ABOUT THE CENTER OF MASS\s+I\(XY\)\s+I\(XZ\)\s+I\(YZ\)\s+([\d.+-Ee-]+)\s+([\d.+-Ee-]+)\s+([\d.+-Ee-]+)", text)

                    if not (mass_match and cg_match and moi_com_match and prod_com_match):
                        raise ValueError("Missing mass property fields in the .dat file.")

                    mass = float(mass_match.group(1))
                    cg = tuple(float(c) for c in cg_match.groups())
                    moi = tuple(float(i) for i in moi_com_match.groups())
                    prod = tuple(float(p) for p in prod_com_match.groups())
                    return mass, cg, moi, prod

                def parallel_axis(mass, cg, ref, moi_cm, prod_cm):
                    dx, dy, dz = cg[0] - ref[0], cg[1] - ref[1], cg[2] - ref[2]
                    Ixx = moi_cm[0]
                    Iyy = moi_cm[1] + mass * (dx**2)
                    Izz = moi_cm[2] + mass * (dx**2)
                    Ixy = prod_cm[0] - mass * dx * dy
                    Ixz = prod_cm[1] - mass * dx * dz
                    Iyz = prod_cm[2] - mass * dy * dz
                    return (Ixx, Iyy, Izz), (Ixy, Ixz, Iyz)

                with open(datfile, 'r') as f:
                    text = f.read()

                mass, cg, moi_cm, prod_cm = parse_block(text)
                moi_ref, prod_ref = parallel_axis(mass, cg, ref_point, moi_cm, prod_cm)
                cg_ref = tuple(c - r for c, r in zip(cg, ref_point))
                cg_ref = (cg_ref[2],0.0, cg_ref[1])
                moi_ref = (moi_ref[2],moi_ref[0], moi_ref[1])


                result = {
                    'mass': mass,
                    'center_of_mass': cg_ref,
                    'reference_point': ref_point,
                    'moment_of_inertia': moi_ref,
                    'product_of_inertia': prod_ref
                }

                def write_mass_properties_csv(csv_out, index, spacing, radius, cg_ref, moi_ref, mass):
                    file_exists = index != 0  # If index is 0, assume new file

                    with open(csv_out, 'a', newline='') as f:
                        writer = csv.writer(f)
                        
                        if not file_exists:
                            writer.writerow(['root_moment','spacing', 'radius', 'x_cg', 'y_cg', 'z_cg', 'Ixx', 'Iyy', 'Izz', 'mass','vm_max','u1','u2','u3'])

                        row = [
                            root_moment,
                            spacing,
                            radius,
                            *cg_ref,
                            *moi_ref,
                            mass,
                            max_vm_stress,
                            max_disp
                        ]
                        writer.writerow(row)
                index =1
                spacing = w
                radius = new_radius
                write_mass_properties_csv(csv_out, index, spacing, radius, cg_ref, moi_ref, mass)

                return result

            data = get_mass_properties_from_dat(
                datfile='130425_final.dat',
                ref_point=(0.0, 0.0, 12.4),
                csv_out='mass_properties.csv'
            )



"""print(data['mass'])
print(data['moment_of_inertia'])
print(data['center_of_mass'])"""