def format_elastic_values(E1e, E2e, G12e, v):
    new_elastic_values = [
        f"{E1e:.2e}, {E2e:.2e}, {E2e:.2e}, {v:.1f}, {v:.1f}, {v:.1f}, {G12e:.2e}, {G12e:.2e}",
        f"{G12e:.2e},"
    ]
    return new_elastic_values

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



formatted = format_elastic_values(9, 99, 999, 0.9)
new_elastic_lines = formatted
update_abaqus_sections_and_material("120425_final.inp", 1, 9, 99,
                                        999, new_elastic_lines,
                                        original_thickness=0.052, original_radius=0.02,
                                        original_density=264.)
    