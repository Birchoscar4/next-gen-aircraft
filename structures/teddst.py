
import os
import csv
# Run Abaqus script
os.system("abaqus python read_odb.py")

# Read the CSV results
results = {}
with open("results.csv", "r") as f:
    reader = csv.reader(f)
    next(reader)  # Skip header
    for row in reader:
        if len(row) < 2:
            continue  # Skip empty or malformed lines
        results[row[0]] = float(row[1])

# Extract values
max_disp = [results["U1"], results["U2"], results["U3"]]
max_vm_stress = results["VM_max"]

# Use the values
print("Max displacement:", max_disp)
print("Max von Mises stress:", max_vm_stress)