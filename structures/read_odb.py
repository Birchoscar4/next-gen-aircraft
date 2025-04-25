def read_odb(filename, node_labels):
    from odbAccess import openOdb
    import csv

    # === Open ODB ===
    odb = openOdb(filename)
    instance = odb.rootAssembly.instances.values()[0]

    # === Map elements to nodes ===
    node_to_elements = {node: set() for node in node_labels}
    for elem in instance.elements:
        for node in elem.connectivity:
            if node in node_labels:
                node_to_elements[node].add(elem.label)

    # === Initialize data storage ===
    results = {node: {"vm_max": 0.0, "u": [0.0, 0.0, 0.0]} for node in node_labels}

    # === Loop through steps and frames ===
    for step in odb.steps.values():
        for frame in step.frames:
            stress_field = frame.fieldOutputs['S']
            disp_field = frame.fieldOutputs['U']

            # --- Displacement updates ---
            for val in disp_field.values:
                if val.instance == instance and val.nodeLabel in node_labels:
                    disp = val.data
                    node = val.nodeLabel
                    for i in range(3):
                        if abs(disp[i]) > abs(results[node]["u"][i]):
                            results[node]["u"][i] = disp[i]

            # --- Stress updates ---
            for val in stress_field.values:
                if val.instance == instance:
                    for node in node_labels:
                        if val.elementLabel in node_to_elements[node]:
                            vm = val.mises
                            if abs(vm) > abs(results[node]["vm_max"]):
                                results[node]["vm_max"] = vm

    odb.close()

    # === Write results to CSV ===
    with open("results.csv", "w") as f:
        writer = csv.writer(f)
        writer.writerow(["node", "vm_max", "u1", "u2", "u3"])
        for node, data in results.items():
            row = [node, data["vm_max"]] + data["u"]
            writer.writerow(row)



read_odb("130425_final.odb",node_labels=[16140, 90, 27377])