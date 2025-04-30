import numpy as np
import pandas as pd
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score

# Load dataset with headers
df = pd.read_csv("mass_properties - Copy.csv")

# Drop all-zero columns
df = df.loc[:, (df != 0).any(axis=0)]

# Define which inputs to use for each output
output_dependency_map = {
    'x_cg': ['w', 'r'],
    'z_cg': ['w', 'r'],
    'Ixx': ['w', 'r'],
    'Iyy': ['w', 'r'],
    'Izz': ['w', 'r'],
    'mass': ['w', 'r'],
    'max_vm': ['w', 'r', 'root_moment'],
    'max_disp': ['w'],
}

# Fit and generate equations
def fit_models_per_output(df, dependency_map, degree=3):
    models = {}
    equations = {}
    
    for output, inputs in dependency_map.items():
        X = df[inputs].values
        y = df[output].values

        poly = PolynomialFeatures(degree=degree)
        X_poly = poly.fit_transform(X)

        model = LinearRegression()
        model.fit(X_poly, y)

        y_pred = model.predict(X_poly)
        r2 = r2_score(y, y_pred)

        feature_names = poly.get_feature_names_out(inputs)
        coef = model.coef_
        intercept = model.intercept_

        equation = f"{output} = {intercept:.6f}"
        for c, name in zip(coef, feature_names):
            if c != 0:
                equation += f" + ({c:.6f})*{name}"

        models[output] = {
            "model": model,
            "poly": poly,
            "r2": r2,
            "inputs": inputs
        }
        equations[output] = equation

    return models, equations

# Run model fitting
models, equations = fit_models_per_output(df, output_dependency_map)

# Print the equations
for name, eq in equations.items():
    print(f"\n{name} equation (R² = {models[name]['r2']:.4f}):\n{eq}")
