from sklearn.linear_model import LinearRegression
import numpy as np
import scipy as sp
import sys
import re
import matplotlib.pyplot as plt
path = r"C:\\temp\\python_self_crated\\packages"
sys.path.insert(0, path)

from allib import general as gl

def load_xyz(path):
    X = []
    Y = []
    Z = []

    with open(path, 'r') as f:
        for line in f:
            if line[0] == '#':
                continue
            else:
                xyz = line.split('\t')
                x = float(xyz[0])
                y = float(xyz[1])
                z = float(xyz[2])

                X.append(x)
                Y.append(y)
                Z.append(z)

    return X, Y, Z

def trim_edge(X, **kwagrs):
    bins = kwagrs.get('bins', 20)

    X = np.array(X)

    _, _, indizes = sp.stats.binned_statistic(X,X, statistic='count', bins=bins)

    bin_begin = indizes == 1
    bin_end = indizes == bins

    min_outer = min(X[bin_begin])
    min_inner = max(X[bin_begin])

    max_outer = max(X[bin_end])
    max_inner = min(X[bin_end])

    return (min_outer, min_inner), (max_outer, max_inner)


# %% INPUT

INPUT = gl.read_input_txt('Input.txt')

def fit_plane(X,Y,Z):


    # Prepare the input data for sklearn
    # Reshape X and Y to be 2D arrays and stack them horizontally
    XY = np.c_[X, Y]

    # Initialize and fit the linear regression model
    model = LinearRegression()
    model.fit(XY, Z)

    # Extract coefficients
    a = model.coef_[0]
    b = model.coef_[1]
    c = model.intercept_

    return a, b, c

def calc_plane_angle(a,b):

    c = 1
    mag_n = np.sqrt(a**2 + b**2 + c**2)

    cos_theta_x = a / mag_n
    cos_theta_y = b / mag_n

    theta_x_rad = np.arccos(cos_theta_x)
    theta_y_rad = np.arccos(cos_theta_y)

    theta_x_deg = np.degrees(theta_x_rad)
    theta_y_deg = np.degrees(theta_y_rad)

    phi_x_deg = 90 - theta_x_deg
    phi_y_deg = 90 - theta_y_deg

    return phi_x_deg, phi_y_deg

def orient_points(path, plot=False, marg=0.05):

    X, Y, Z = load_xyz(path)
    a, b, c = fit_plane(X, Y, Z)
    phi_x, phy_y = calc_plane_angle(a, b)
    marg_abs_x = max(X) * marg
    marg_abs_y = max(Y) * marg
    x_zone = (min(X) + marg_abs_x, max(X) - marg_abs_x)
    y_zone = (min(Y) + marg_abs_y, max(Y) - marg_abs_y)

    if plot:
        # Visualization of the points and the fitting plane
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        # Scatter plot of the points
        ax.scatter(X, Y, Z, color='b', label='Data Points', s=0.1, alpha=0.1)
        # Create a grid to plot the plane
        x_range = np.linspace(min(X), max(X), 10)
        y_range = np.linspace(min(Y), max(Y), 10)
        x_grid, y_grid = np.meshgrid(x_range, y_range)
        z_grid = a * x_grid + b * y_grid + c
        # Plot the plane
        ax.plot_surface(x_grid, y_grid, z_grid, color='r', label='Fitted Plane')
        # Set labels
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        ax.set_title('fitted points with \n' + f'angle from x: {phi_x:.6f}°, angle from y {phy_y:.6f}°' + '\n' + f'z - offset at (x=0, y=0): {c:.6f}m')
        # Show the plot
        # plt.show()

    out = {'phi_x': phi_x, 'phi_y': phi_x, 'z_offset': c, 'x_range': x_zone, 'y_range': y_zone}

    if plot:
        return out, fig
    else:
        return out
    
def replace_variables_ansysScript(code_text: str, replacements: dict) -> str:
    # Loop through the replacements dictionary
    for name, new_value in replacements.items():
        # Regex pattern to match 'name = value' with optional spaces around '='
        pattern = rf'(\b{name}\b)\s*=\s*.*'
        # Replace with 'name = new_value'
        replacement = rf'\1 = {new_value}'
        # Apply the substitution
        code_text = re.sub(pattern, replacement, code_text)

    return code_text