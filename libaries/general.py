import datetime
import os
import random
import sqlite3
import subprocess
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
import scipy as sc
import sklearn as skl
from sklearn.linear_model import LinearRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import PolynomialFeatures
from matplotlib.backends.backend_pdf import PdfPages
import decimal


def model_regression(x: pd.core.series.Series, y: pd.core.series.Series,
                     **kwargs) -> skl.linear_model._base.LinearRegression:
    """creates regression model of the given x- and y-Series of the certain degree with "include_bias = False"

    Parameters
    ----------
    x : pd.core.series.Series
        x-data.
    y : pd.core.series.Series
        y-data.

    optional:

    degree: int
        default: 4

    Returns
    -------
    model : skl.linear_model._base.LinearRegression
        regresssion-modell.

    """

    from sklearn.preprocessing import FunctionTransformer

    method = kwargs.get('method', 'poly')

    degree = kwargs.get('degree', 3)

    weights = kwargs.get('weights', float(0))

    reg_free_n = kwargs.get('reg_free_n', False)

    weights_regulation = kwargs.get('weights_regulation', float(1))

    x = np.array(x).reshape((-1, 1))
    y = np.array(y)

    functions = {
        'poly': (lambda f: f),
        'sqrt': (lambda f: np.sqrt(f)),
        'log': (lambda f: np.where(f > 0, np.log(f + 1), f)),
        'sin': (lambda f: np.sin(f))
    }

    model = make_pipeline(
        PolynomialFeatures(degree, include_bias=True),
        FunctionTransformer(functions[method]),
        LinearRegression(fit_intercept=reg_free_n)
    )

    # change weight
    if weights_regulation is not None:
        mean_weight = weights.mean()
        weights = weights + (1 - weights_regulation) * (mean_weight - weights)

    model.fit(x, y, linearregression__sample_weight=weights)

    linear_reg = model.named_steps['linearregression']
    coefficients = linear_reg.coef_

    return model, coefficients


def k_aus_omega(omega: np.ndarray, d: float) -> np.ndarray:
    """Itteratives lösen der Dispersionsgleichung nach der Wellenzahl
    
    berechnet aus omega als Einzelwert oder aus np.1darray 

    Parameters
    ----------
    omega: np.1darray oder float - 
        Kreisfrequenzen der zu bestimmenden Wellenzahlen.    
    d: float
        Wellenzahl
        
    Returns
    -------
    k:  korrespondierende Wellenzahl nach lösen der Dispersionsgleichung 
        -.  
        
    Notes
    -----    
    itteratives Lösen der Dispersionsgleichung der linearen Wellentheorie
    k = omega**2 / (np.tanh(k * d) * g)
   
        """
    g = 9.8070

    # Omega = 0 führt zu Fehlern in der Dispersionsgleichung und wird deshalb auf sehr kleinen Wert gesetzt
    if isinstance(omega, np.ndarray):
        omega[omega == 0] = 0.001

    # Startwert aus Tiefwasserbedingung
    k = omega ** 2 / g
    k_alt = 0
    i = 0

    # Itteration über Konvergenzgleichung bis zur Konvergenz oder zu vielen Iterationen
    while np.max(np.abs(k_alt - k)) > 0.001 and i < 1000:
        k_alt = k
        k = omega ** 2 / (np.tanh(k * d) * g)
        i += 1

    return k


def predict_regression(model: skl.linear_model._base.LinearRegression, x_pred: np.ndarray) -> np.ndarray:
    """creates the regression curve of the given regression-model "model" of the points x_pred and outputs the y-values
    y_pred as np.ndarray

    Parameters
    ----------
    model : skl.linear_model._base.LinearRegression
        regression model.

    x_pred : np.ndarray
        x-points to predict by regression-model.

    Returns
    -------
    y_pred : np.ndarray
         corresponding predicted y-points
    """

    x_pred = np.array(x_pred)

    x_pred = x_pred.reshape((-1, 1))

    y_pred = model.predict(x_pred)

    return y_pred


def filter_dataframe(dataframe, column_names: str | list, a_min: float | list,
                     a_max: float | list) -> pd.core.frame.DataFrame | list:
    """filters "dataframe" on the basis of one or more colums specified in column_names in the boundarys given in the corresponding entrys of "a_min" and "a_max"

    Parameters
    ----------
    dataframe : pd.core.frame.DataFrame
        dataframe (or series) to filter.

    column_names: str or list of stings
        names/keys of columns in dataframe to filter by.

    a_min: str or list
        lower boundaray of corresponduing values in column

    a_max: str or list
        upper boundaray of corresponduing values in column


    Returns
    -------
    df_out : pd.core.frame.DataFrame
         filtered dataframe

    """

    if type(dataframe) is pd.Series:
        dataframe = pd.DataFrame(dataframe)

    if type(column_names) is str:
        column_names = [column_names]
        a_min = [a_min]
        a_max = [a_max]

    df_out = dataframe

    for i, names in enumerate(column_names):
        if ~np.isnan(a_max[i]):
            if a_min[i] > a_max[i]:
                df_out = df_out.loc[(dataframe[names] >= a_min[i])
                                    | (dataframe[names] < a_max[i])]

            else:
                df_out = df_out.loc[(dataframe[names] >= a_min[i])
                                    & (dataframe[names] < a_max[i])]

        else:
            df_out = df_out.loc[(dataframe[names] >= a_min[i])]
    return df_out


def grid_pointcloud_in_x(x, y, grid, **kwargs):
    """sorts x,y points in bins in x-axis. Bin-edges given by "grid". Saves them into pandas Dataframe and calculates  Count, Mean, Standard Deviation, and the square weigthed mean. The last one is calculated by sqaring every value, adding them up, dividing them by the number of points and taking the squareroot 

    Parameters
    ----------
    x : pd.core.series.Series
        x-coordinate of points
    y : pd.core.series.Series
        x-coordinate of points
    grid : np.ndarray | list
        grid .
    
    *kwargs
    method: str, choose from: {'mean', 'weighted mean', 'median'}

    Returns
    -------
    averaged: pd.Series with averaged result over x-grid
    std: pd.Series with standatd deviation over x-grid
    count: pd.Series with count over x-grid
    """

    x_data = np.array(x)
    y_data = np.array(y)

    if len(x) == 0:
        x_data = float('nan')
        y_data = float('nan')

    method = kwargs.get('method', False)

    def mean_weight(i):
        # i[i == 0] = float('nan')
        mean_w = np.sqrt(np.sum(i ** 2) / np.size(i[~np.isnan(i)]))
        return mean_w

    if method == 'mean':
        averaged, x_edges, _ = sc.stats.binned_statistic(
            x_data, y_data, statistic='mean', bins=grid)

    elif method == 'weighted mean':
        averaged, x_edges, _ = sc.stats.binned_statistic(
            x_data, y_data, statistic=mean_weight, bins=grid)

    elif method == 'median':
        averaged, x_edges, _ = sc.stats.binned_statistic(x_data, y_data, statistic='median', bins=grid)

    else:
        print("choose method from: {'mean', 'weighted mean', 'median'}")
        return

    std, _, bin_ident = sc.stats.binned_statistic(
        x_data, y_data, statistic='std', bins=grid)

    count, _, _ = sc.stats.binned_statistic(
        x_data, y_data, statistic='count', bins=grid)

    x_grid = (grid[:-1] + grid[1:]) / 2

    averaged = pd.Series(averaged, index=x_grid)

    std = pd.Series(std, index=x_grid)

    count = pd.Series(count, index=x_grid)

    return averaged, std, count, bin_ident


def JONSWAP(f: list, T_p: float, H_s: float, gamma_mode='torset') -> list:
    """calculates the JONSWAP Spectrum
        
    Parameters
    ----------
    f : list
        input Frequency (x-Axis)
    T_p: float
        WavePeriod 
    H_s: float
        signifcant WaveHeight 
    
    Returns
    -------
    S : list
        coresponding spectal values"""

    f_p = 1 / T_p

    krit = np.sqrt(T_p) / H_s

    if gamma_mode == 'default':
        if krit <= 3.6:
            gamma = 5
        if 3.6 < krit <= 5:
            gamma = np.exp(5.75 - 1.15 * krit)
        if krit > 5:
            gamma = 1

    if gamma_mode == 'torset':
        gamma = 35 * ((2 * np.pi * H_s) / (9.81 * T_p ** 2)) ** (6 / 7)

    sigma_under = 0.07

    sigma_over = 0.09

    S_under = 0.3125 * H_s ** 2 * T_p * (f / f_p) ** -5 * np.exp(-1.25 * (f / f_p) ** -4) * (
            1 - 0.287 * np.log(gamma)) * gamma ** np.exp(-0.5 * ((f - f_p) / (sigma_under * f_p)) ** 2)

    S_over = 0.3125 * H_s ** 2 * T_p * (f / f_p) ** -5 * np.exp(-1.25 * (f / f_p) ** -4) * (
            1 - 0.287 * np.log(gamma)) * gamma ** np.exp(-0.5 * ((f - f_p) / (sigma_over * f_p)) ** 2)

    S = np.concatenate((S_under[f < f_p], S_over[f > f_p]))

    return S


def c_scatterplot(x: pd.core.series.Series, y: pd.core.series.Series):
    """calculates c-Vektor, conatining densitiy of point-distribution of the corresponding x,y pair.
    Can be used as input for scatter(x, y, c = c, cmap = ...)
        
    Parameters
    ----------
    x : pd.core.series.Series | np.array | list
        x - cooronates of point
    y : pd.core.series.Series | np.array | list
        y - cooronates of point
    
    Returns
    -------
    c : list
        coresponding c values"""

    x = np.array(x)
    y = np.array(y)

    # Punkthäufigkeit berechnen
    z, xedges, yedges = np.histogram2d(x, y, bins=200)

    # Farben für jeden Punkt festlegen
    colors = []
    for i in range(len(x)):
        bin_x = np.searchsorted(xedges, x[i]) - 1
        bin_y = np.searchsorted(yedges, y[i]) - 1
        colors.append(z[bin_x, bin_y])

    return colors


def linear_extrapolation_log_x(x, y, x_end):
    """
    Perform linear extrapolation on a log scale (x-axis) with a linear y-axis.

    Parameters:
    x (array-like): Original x data points.
    y (array-like): Original y data points.
    x_end (float): The end value for the x-direction to which extrapolation should occur.

    Returns:
    combined_x (numpy.ndarray): Combined original and extrapolated x data points.
    combined_y (numpy.ndarray): Combined original and extrapolated y data points.
    extrapolated (list of bool): Boolean list indicating which points are extrapolated.
    """
    # Logarithmic transformation for x-axis
    log_x = np.log(x)

    # Perform linear fit on the last segment of the data
    slope, intercept = np.polyfit(log_x[-2:], y[-2:], 1)

    # Generate new x values for extrapolation
    new_log_x = np.logspace(np.log10(x[-1]), np.log10(x_end), num=2)
    new_y = slope * np.log(new_log_x) + intercept

    # Combine original and extrapolated data
    combined_x = np.concatenate((x, new_log_x[1:]))  # Avoid duplicating the last original x value
    combined_y = np.concatenate((y, new_y[1:]))  # Avoid duplicating the last original y value

    # Generate boolean list for extrapolation
    extrapolated = [False] * len(x) + [True] * (len(new_log_x) - 1)

    return combined_x, combined_y, extrapolated


def calculate_histogram(x, **kwargs):
    x_max = kwargs.get('x_max', max(x))
    bin_size_fix = kwargs.get('bin_size_fix', None)
    x_min = kwargs.get('x_min', min(x))

    # get significant Digits
    x_str = [str(value) for value in x]
    lenths = [len(dig.split('.')[1]) for dig in x_str]
    sig_dig = max(lenths)
    if sig_dig > 3: sig_dig = 3
    # get minmal distance
    x_values = np.sort(x)
    # x_max = np.max(x_values)
    x_unique = np.unique(np.round(x_values, sig_dig))
    bin_size = round(min(x_unique[1:] - x_unique[:-1]), sig_dig)
    N_bins = int(np.ceil(x_max / bin_size))
    x_max = N_bins * bin_size
    if x_min < 0:
        N_bins_down = int(np.ceil(-x_min / bin_size))
        x_min = -N_bins_down * bin_size
        N_bins = N_bins + N_bins_down
    if x_min > 0:
        N_bins_down = int(np.floor(-x_min / bin_size))
        x_min = N_bins_down * bin_size
        N_bins = N_bins - N_bins_down
    if N_bins > 10000:
        bin_size = np.round(x_max / 1000, 3)
        N_bins = int(x_max / bin_size)

    # choose, if to supplied binsize is to small:

    if bin_size_fix is not None:
        bin_size = bin_size_fix
        N_bins = int(x_max / bin_size)

    edges = np.linspace(x_min - bin_size / 2, x_max + bin_size / 2, N_bins + 2)
    edges = np.round(edges, 4)
    count = [len(x_values[(x_values > edges[i]) & (x_values <= edges[i + 1])]) for i in range(len(edges) - 1)]
    center = (edges[1:] + edges[:-1]) / 2

    return bin_size, center, count


def fit_weibull_distribution(data, x_grid, **kwargs):
    """
    Fits a Weibull distribution to the given data

    Args:
        data (array-like): The data to fit the Weibull distribution to.
        x_grid (1D array-like): the grid on which to write the weibull distribution
    Optional:
        floc: loc parameter (0: is fixed at (0,0))

    Returns:
        weibull: the matching weibull distribution on x_grid
        params: A dictionary containing the parameters of the fitted Weibull distribution.
                {"shape": shape, "loc": loc, "scale": scale, "mean": mean, "std": std}
    """

    floc = kwargs.get("floc", None)
    # Fit the Weibull distribution to the data
    if floc is not None:
        shape, loc, scale = sc.stats.weibull_min.fit(data, floc=floc)
    else:
        shape, loc, scale = sc.stats.weibull_min.fit(data)

    # Create a Weibull distribution object
    weibull_dist = sc.stats.weibull_min(shape, loc=loc, scale=scale)

    # Calculate the mean and standard deviation of the fitted distribution
    mean = weibull_dist.mean()
    std = weibull_dist.std()

    weibull = weibull_dist.pdf(x_grid)
    params = {"shape": shape, "loc": loc, "scale": scale, "mean": mean, "std": std}

    return weibull, params


def write_dict(data):
    """
    Recursively writes all data of nested dictionaries to a string,
    separated with newlines and tabulators.

    :param data: The dictionary to process
    :return: A string representation of the dictionary
    """

    def process_dict(data, indent_level=0):
        result = ""
        indent = "\t" * indent_level

        for key, value in data.items():
            if isinstance(value, dict):
                result += f"{indent}{key}:\n"
                result += process_dict(value, indent_level + 1)
            else:
                result += f"{indent}{key}: {value}\n"

        return result

    return process_dict(data)


def interpolate_increasing_decreasing(new_x_values, x_values, y_values, kind='linear'):
    """
    Interpolates values on a given function defined by lists of x and y values.

    Parameters:
    - x_values: List or array of x values (must be in a monotonically increasing or decreasing order)
    - y_values: List or array of y values corresponding to x_values
    - new_x_values: List or array of new x values where interpolation is required
    - kind: Type of interpolation (default is 'linear'). Options include 'linear', 'nearest', 'zero', 'slinear', 'quadratic', 'cubic'.

    Returns:
    - Interpolated y values corresponding to new_x_values, with 'nan' for out-of-bounds values
    """

    # Ensure the x and y values are numpy arrays
    x_values = np.array(x_values)
    y_values = np.array(y_values)

    # Check if the function is decreasing and reverse if necessary
    if x_values[0] > x_values[-1]:
        x_values = x_values[::-1]
        y_values = y_values[::-1]

    # Create the interpolation function with 'nan' for out-of-bounds values
    interpolation_function = sc.interpolate.interp1d(x_values, y_values, kind=kind, bounds_error=False, fill_value=np.nan)

    # Interpolate the new x values
    new_y_values = interpolation_function(new_x_values)

    return new_y_values


def read_input_txt(file_path, encoding='utf-8'):
    data_dicts = []
    dict_names = []
    cd = {}
    dict_num = 0
    DICT_OUT = {}

    with open(file_path, 'r', encoding=encoding) as file:
        for line in file:

            # auskommentieren
            line = line.split('#', 1)[0]
            # Check for quotation marks to start a new dictionary

            if '*' in line:

                # Extract the name between quotation marks

                section_name = line[line.index('*') + 1:line.rindex('*')]
                dict_names.append(section_name)
                # Save the current dictionary and start a new one
                if dict_num != 0:
                    data_dicts.append(cd)

                cd = {}
                dict_num = dict_num + 1
            else:
                # Split the line into key and value using one or more spaces
                parts = line.split('>')

                for i, part in enumerate(parts):
                    parts[i] = part.strip()

                # Ensure the line is not empty
                if np.size(parts) > 1:
                    key = parts[0]
                    value = parts[1]
                    try:
                        value = eval(value)
                    except (SyntaxError, NameError, TypeError) as e:
                        print(
                            f"    WARNING: {e}, Error in Input File Value '{key}' in Section '{section_name}', entry set to 'None'")
                        value = None
                    # Add the key-value pair to the current dictionary
                    cd[key] = value

    # Append the last dictionary after reading the file
    data_dicts.append(cd)

    i = 0
    for name in dict_names:
        DICT_OUT[name] = data_dicts[i]
        i = i + 1

    return DICT_OUT


def string_to_latex(string):
    """converts sting to latex string by converting it into a raw string and putting $ at the begining and end"""

    string = repr(string)
    string = string.replace(r'\n', r'\\')
    string = "$" + string + "$"
    return string


import pandas as pd


def xlsx2dict(file_path):
    """
    Reads an Excel file and returns a dictionary containing each sheet as a DataFrame.

    Each sheet in the Excel file will be stored as a DataFrame in the dictionary,
    with the sheet name as the key. The first row of each sheet is used as column headers,
    and the first column is set as the index for each DataFrame.

    Parameters:
    ----------
    file_path : str
        The path to the Excel (.xlsx) file.

    Returns:
    -------
    dict
        A dictionary where each key is a sheet name, and each value is a DataFrame
        representing the sheet's data, with headers and indices correctly set.

    """

    # Read all sheets into a dictionary with pandas
    xls = pd.ExcelFile(file_path)
    sheets_dict = {}

    for sheet_name in xls.sheet_names:
        # Read each sheet, using the first row as headers and the first column as the index
        df = pd.read_excel(xls, sheet_name=sheet_name, header=0, index_col=0)
        sheets_dict[sheet_name] = df

    return sheets_dict


def export_df_from_sql(db_file, table_name, column_names=None, timeframe=None, indizes=None):
    """
    Load specified columns or all columns from a table in an SQLite database and return as a pandas DataFrame,
    if tables index is specified, it is set as the index of the dataframe

    Parameters:
    - db_file (str): Path to the SQLite database file.
    - table_name (str): Name of the table to query.
    - column_names (list of str, optional): List of column names to retrieve. If None, retrieves all columns.
    - timeframe: list, optional: [time_start, time_end], must be pd.DateTime object, also indizes of dataframe must be read as datetime objects, is ignored otherwise
    - indizes: list, optional: list of indizes, (dependent on index format of df) is applied AFTER timeframe! 
    
    Returns:
    - pandas.DataFrame: DataFrame containing the query results, including the row index as a column.

    Raises:
    - ValueError: If the table or columns do not exist.
    """
    try:
        conn = sqlite3.connect(db_file)
        cursor = conn.cursor()

        # Check if table exists
        cursor.execute(f"SELECT name FROM sqlite_master WHERE type='table' AND name='{table_name}' ORDER BY ROWID;")
        if not cursor.fetchone():
            raise ValueError(f"Table '{table_name}' does not exist in the database.")

        # Check if specified columns exist in the table
        cursor.execute(f"PRAGMA table_info({table_name});")
        columns_info = cursor.fetchall()
        table_columns = [info[1] for info in columns_info]

        if column_names:
            for col in column_names:
                if col not in table_columns:
                    raise ValueError(f"Column '{col}' does not exist in table '{table_name}'.")
            # Enclose column names with special characters in double quotes
            cols = ', '.join([f'"{col}"' for col in column_names])

        else:
            # Enclose column names with special characters in double quotes
            cols = ', '.join([f'"{col}"' for col in table_columns[1:]])

        # Add ROWID to the columns to retrieve
        query = f"SELECT {cols} FROM {table_name};"
        df = pd.read_sql_query(query, conn)

        # set index of sql database to index of dataframe
        cursor.execute(f"PRAGMA index_list('{table_name}');")
        indexes = cursor.fetchall()
        cursor.execute(f"PRAGMA index_info('{indexes[0][1]}');")

        index_info = cursor.fetchall()

        if len(index_info) > 0:

            index_col_name = index_info[0][2]

            cursor.execute(f'SELECT "{index_col_name}" FROM "{table_name}" ORDER BY ROWID;')
            index_col = cursor.fetchall()

            index_col = [info[0] for info in index_col]

            try:
                index_col = pd.to_datetime(index_col)
            except:
                print("could not convert index to datetime object")
                index_col = index_col
            df.index = index_col

            if timeframe is not None:
                df = df.loc[timeframe[0]:timeframe[1]]

        if indizes is not None:
            df = df.loc[indizes]

        return df

    except sqlite3.Error as e:
        print(f"An error occurred: {e}")
        raise
    finally:
        if conn:
            conn.close()


def export_colnames_from_db(database_path):
    # Connect to the SQLite database
    conn = sqlite3.connect(database_path)
    cursor = conn.cursor()

    # Dictionary to store table names and their respective column names
    table_columns = {}

    # Get the list of all tables in the database
    cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
    tables = cursor.fetchall()

    # Iterate over each table and get the column names
    for table in tables:
        table_name = table[0]

        # Use double quotes to escape the table name
        cursor.execute(f'PRAGMA table_info("{table_name}");')
        columns = cursor.fetchall()
        column_names = [column[1] for column in columns]  # Column names are in the 2nd index

        table_columns[table_name] = column_names

    # Close the connection
    conn.close()

    return table_columns


def auto_ticks(start, end, num_ticks=10, fix_end=False, edges=False):
    """generates an array of tick values similar to what would be automatically determined by Matplotlib for plotting purposes, given a specified start and end range. The function provides an optional feature to ensure that the last tick always aligns precisely with the end value.

    Parameters:

    start (float or int): The starting value of the range for which ticks should be generated.
    end (float or int): The ending value of the range for which ticks should be generated.
    num_ticks (int, optional): The maximum number of ticks to generate. The default is 10.
    fixed_end (bool, optional): If True, the last tick will be exactly equal to the end value. If False, the last tick is determined by Matplotlib's automatic algorithm. The default is False.
    edges (bool, optional): If True, outputs the edges with the tick as the midpoint. The default is False.
    
    Returns:

    ticks (numpy.ndarray): An array of tick values within the specified range."""

    # Create a MaxNLocator instance with a specified maximum number of ticks
    locator = ticker.MaxNLocator(nbins=num_ticks)

    # Generate the ticks using the locator
    ticks = locator.tick_values(start, end)

    if fix_end:
        # # Ensure the last tick is exactly the end value
        ticks = np.append(ticks[:-1], end)

    if edges:
        # Calculate the bin edges based on the ticks (treating ticks as midpoints)
        bin_width = np.diff(ticks) / 2
        edges = np.concatenate(([ticks[0] - bin_width[0]], ticks[:-1] + bin_width, [ticks[-1] + bin_width[-1]]))
        return ticks, edges

    else:
        return ticks


def range_stepfix(step, zone):
    """creates a range from zone[0] to at least zone[1]  and keeps the step width and does one step over the upper lim if necessary"""

    vm_length = zone[1] - zone[0]
    N_inner = int(np.ceil(vm_length / step) + 1)
    edges = np.linspace(
        zone[0], zone[0] + (N_inner - 1) * step, N_inner)
    return edges


def read_lua_values(file_path, keys):
    """
    Extracts specified key-value pairs from a Lua file.

    Parameters:
    ----------
    file_path : str
        The path to the Lua file from which values need to be extracted.

    keys : list of str
        A list of keys (as strings) for which the corresponding values should be extracted from the Lua file.

    Returns:
    -------
    dict
        A dictionary where the keys are the specified keys from the input list, and the values are the corresponding values
        found in the Lua file. The values are converted to their appropriate types (int, float, bool, or str) based on their
        format in the Lua file.
    """

    # Dictionary to store the values
    values_dict = {}

    # Regular expression to match key-value pairs in the Lua file
    pattern = re.compile(r'(\w+)\s*=\s*(.+)')

    with open(file_path, 'r') as file:
        for line in file:
            # Remove comments and strip any leading/trailing whitespace
            line = line.split('--')[0].strip()
            if not line:
                continue

            # Match the key-value pair
            match = pattern.match(line)
            if match:
                key, value = match.groups()

                # Remove trailing comma if present
                value = value.rstrip(',')

                # Check if the key is in the desired keys
                if key in keys:
                    # Attempt to convert to a number or leave as a string
                    try:
                        # Handle numbers and booleans
                        if value.lower() == "true":
                            value = True
                        elif value.lower() == "false":
                            value = False
                        elif "." in value:
                            value = float(value)
                        else:
                            value = int(value)
                    except ValueError:
                        # Keep as string if conversion fails
                        value = value.strip('"').strip("'")

                    values_dict[key] = value

    return values_dict


def write_lua_variables(file_path, variables):
    """writes lua varables found at 'filepath', specified by dict variables = {'var name': value, ...} """

    with open(file_path, 'r') as file:
        lines = file.readlines()
    for variable_name, new_value in variables.items():
        variable_pattern = re.compile(rf'^(\s*{variable_name}\s*=\s*).*?(\s*,\s*--.*)?$')

        for i, line in enumerate(lines):
            match = variable_pattern.match(line)
            if match:
                indentation = match.group(1)
                rest_of_line = match.group(2) if match.group(2) else ','
                lines[i] = f"{indentation}{new_value}{rest_of_line}\n"
                break

    with open(file_path, 'w') as file:
        file.writelines(lines)


def format_array(arr, formatter):
    """
    Formats the elements of a numpy array according to the provided formatter.

    Parameters:
    arr (numpy.ndarray): The input numpy array to be formatted.
    formatter (str): The format specifier, e.g., ".2f" for 2 decimal places.

    Returns:
    numpy.ndarray: A numpy array with formatted strings.
    """
    # Use vectorized string formatting for the entire array
    formatted_array = np.vectorize(lambda x: f"{x:{formatter}}")(arr)

    return formatted_array


def angle_deviation(angles_base, angles_comp):
    """"calculates the angle_deviation from angles_base and angles_comp in deg, always taking the smallest angle between the two values

    Input:
    angles_base: numpy array, list or series:  angle to base comparison on
    angles_comp: numpy array, list or series:  angle to compare to base

    Retruns:
    diff: numpy array: angle difference
    diff_abs: numpy array: angle difference abs

    Notes:
        diff =  angles_comp - angles_base, chossing negative or positive, which is smaller in abs
    """
    angles_base = np.array(angles_base)
    angles_comp = np.array(angles_comp)

    PITAU = 360 + 180  # for readablility
    DIFF_ANG = (angles_comp - angles_base + PITAU) % 360 - 180
    DIFF_ANG_ABS = abs(DIFF_ANG)

    return DIFF_ANG, DIFF_ANG_ABS


def angle_midpoints(angle1, angle2):
    """calculates the midpoint between pairs of angles provided in two input sequences.
    The angles are assumed to range between 0 and 360 degrees, and the midpoint is always calculated in the clockwise direction from the first angle to the second."""

    if isinstance(angle1, (list, np.ndarray, pd.Series)):
        angle1 = list(angle1)
    else:
        raise TypeError("list2 Input must be a list, numpy array, or pandas Series.")

    if isinstance(angle2, (list, np.ndarray, pd.Series)):
        angle2 = list(angle2)
    else:
        raise TypeError("Input must be a list, numpy array, or pandas Series.")

    if len(angle1) != len(angle2):
        raise ValueError("The lists must have the same length.")

    midpoints = []

    for angle1, angle2 in zip(angle1, angle2):
        angle1 = angle1 % 360
        angle2 = angle2 % 360
        # Calculate the difference
        diff = (angle2 - angle1) % 360
        # Calculate the midpoint by moving half the difference
        midpoint = (angle1 + diff / 2) % 360
        midpoints.append(midpoint)

    return midpoints


import re


def alias(input_data, original, alias):
    """Replaces every instance in the input string that matches the longest possible substring from the `original` dict values with the corresponding values in the `alias` dict, ensuring each part is replaced only once."""

    def replace_in_string(input_string):
        # Sort keys by the length of their corresponding values in the original dictionary in descending order
        sorted_keys = sorted(original, key=lambda k: len(original[k]), reverse=True)

        # Unique placeholder generator to avoid collisions
        placeholder_format = "__REPLACEMENT_PLACEHOLDER_{}__"
        placeholder_map = {}
        replacement_count = 0

        # First, replace using placeholders to avoid recursive changes
        for key in sorted_keys:
            if key in alias:
                # Generate a unique placeholder
                placeholder = placeholder_format.format(replacement_count)
                placeholder_map[placeholder] = alias[key]
                replacement_count += 1

                # Replace original value with the placeholder
                pattern = re.escape(original[key])
                input_string = re.sub(pattern, placeholder, input_string)

        # Finally, replace placeholders with their actual values
        for placeholder, replacement_value in placeholder_map.items():
            if replacement_value is not None:
                input_string = input_string.replace(str(placeholder), str(replacement_value))

        return input_string

    # Handle the input_data being a list or a string
    if isinstance(input_data, list):
        return [replace_in_string(s) for s in input_data]
    elif isinstance(input_data, str):
        return replace_in_string(input_data)
    else:
        raise ValueError("Input must be either a string or a list of strings.")


def save_figs_as_png(FIG, filename,dpi=600, lualatex_mode=False):


    i = 1
    for fig in FIG:
        fig.savefig(filename + f"_page_{i}.png", dpi=dpi)
        i = i + 1
        plt.close(fig)
    return

# def save_figs_as_png(figures, path_in, dpi=600):
#     """
#     Saves a list of Matplotlib figures as PNGs using a PDF file path.
#
#     Args:
#         figures (list): List of Matplotlib figure objects to save.
#         pdf_path (str): Path for the temporary PDF file (also used as a base name for PNGs).
#         dpi (int, optional): DPI for the PNG output. Defaults to 600.
#
#     Returns:
#         list: List of file paths for the saved PNG files.
#     """
#     pdf_path = path_in + '.pdf'
#
#     png_paths = []
#     pdf_dir, pdf_filename = os.path.split(pdf_path)
#     base_name = os.path.splitext(pdf_filename)[0]
#
#
#     # Save figures to a PDF file at the specified path
#     with PdfPages(pdf_path) as pdf:
#         for fig in figures:
#             pdf.savefig(fig, dpi=dpi)
#
#     # Open the PDF with fitz to convert each page to PNG
#     with fitz.open(pdf_path) as pdf_document:
#         for page_num in range(pdf_document.page_count):
#             page = pdf_document.load_page(page_num)
#             png_path = os.path.join(pdf_dir, f"{base_name}_page_{page_num + 1}.png")
#             page_pix = page.get_pixmap(dpi=dpi)
#             page_pix.save(png_path)
#             png_paths.append(png_path)
#
#     # Remove the PDF file after PNG conversion
#   #  os.remove(pdf_path)
#
#     return png_paths


def save_figs_as_pdf(FIG, filename, **kwargs):
    dpi = kwargs.get('dpi', 600)

    with PdfPages(filename + ".pdf") as pdf:
        for fig in FIG:
            pdf.savefig(fig, dpi=dpi)
            plt.close(fig)
    return


import numpy as np


# def make_monotone(x, increasing=True):
#     """Modifies an input sequence (either a Python list or a NumPy array) to make it monotone, ensuring that the sequence is either non-decreasing (monotone increasing) or non-increasing (monotone decreasing) based on an optional parameter.
#
#     Parameters:
#     x (list or numpy.ndarray):
#         The input sequence to be modified. This can be either a Python list or a NumPy array containing numeric elements.
#
#     increasing (bool, optional):
#         A boolean flag that determines the desired monotonicity of the sequence:
#         True (default): Modifies the sequence to be monotone increasing, where each element is at least as large as the one before it.
#         False: Modifies the sequence to be monotone decreasing, where each element is at most as large as the one before it.
#
#     Returns:
#     numpy.ndarray:
#         A NumPy array containing the modified sequence, either monotone increasing or decreasing based on the increasing parameter. The modification is done in place, meaning the original sequence is altered."""
#
#     # Convert to a NumPy array if it isn't one already
#     arr = np.array(x, copy=False)
#
#     # Iterate over the array starting from the second element
#     for i in range(1, len(arr)):
#         if np.isnan(arr[i]):  # Skip NaN values
#             continue
#
#         if increasing:
#             # For monotone increasing, replace current element if it's smaller than the previous one
#             if i > 0 and not np.isnan(arr[i - 1]) and arr[i] < arr[i - 1]:
#                 arr[i] = arr[i - 1]
#         else:
#             # For monotone decreasing, replace current element if it's larger than the previous one
#             if i > 0 and not np.isnan(arr[i - 1]) and arr[i] > arr[i - 1]:
#                 arr[i] = arr[i - 1]
#
#     return arr


def make_monotone(x, increasing=True):
    """modifies an input sequence (either a Python list or a NumPy array) to make it monotone, ensuring that the sequence is either non-decreasing (monotone increasing) or non-increasing (monotone decreasing) based on an optional parameter.

    Parameters:
    x (list or numpy.ndarray):
    The input sequence to be modified. This can be either a Python list or a NumPy array containing numeric elements.

    increasing (bool, optional):
    A boolean flag that determines the desired monotonicity of the sequence:

    True (default): Modifies the sequence to be monotone increasing, where each element is at least as large as the one before it.
    False: Modifies the sequence to be monotone decreasing, where each element is at most as large as the one before it.
    Returns:
    numpy.ndarray:
    A NumPy array containing the modified sequence, either monotone increasing or decreasing based on the increasing parameter. The modification is done in place, meaning the original sequence is altered."""

    # Convert to a NumPy array if it isn't one already
    arr = np.array(x, copy=False)

    # Iterate over the array starting from the second element
    for i in range(1, len(arr)):
        if increasing:
            # For monotone increasing, replace current element if it's smaller than the previous one
            if arr[i] < arr[i - 1]:
                arr[i] = arr[i - 1]
        else:
            # For monotone decreasing, replace current element if it's larger than the previous one
            if arr[i] > arr[i - 1]:
                arr[i] = arr[i - 1]

    return arr


def get_group_labels(series):
    return (series != series.shift()).cumsum() * series


def compare_values(value1, value2, operation="=="):
    """
    Compare two values based on the specified operation.

    Parameters:
        value1: The first value to compare. Can be of any type.
        value2: The second value to compare. Can be of any type.
        operation (str): The comparison operation to perform.
                         Supported operations are '==', '>', and '<'.
                         Defaults to '=='.

    Returns:
        - If either value is None, returns True.
        - If the types of value1 and value2 differ, returns False.
        - If the types match, returns the result of the comparison based on the specified operation:
            - '==': Returns True if value1 equals value2, otherwise False.
            - '>': Returns True if value1 is greater than value2, otherwise False.
            - '<': Returns True if value1 is less than value2, otherwise False.

    Raises:
        ValueError: If an unsupported operation is provided.
    """
    # If either value is None, return True
    if bool(value1 is None) or bool(value2 is None):
        return True

    if value1 is None and value2 is None:
        return False

    # If the types of the values are different, return False
    if type(value1) != type(value2):
        return False

    # Perform the comparison based on the selected operation
    if operation == "==":
        return value1 == value2
    elif operation == ">":
        return value1 > value2
    elif operation == "<":
        return value1 < value2
    else:
        raise ValueError("Invalid operation. Supported operations are '==', '>', and '<'.")


def round_to_significant_digit(arr, sig_digits):
    """
    Rounds each element in the numpy array to the specified number of significant digits and
    returns the results as strings.

    Parameters:
    arr (np.ndarray): The input numpy array of floats.
    sig_digits (int): The number of significant digits to round to.

    Returns:
    np.ndarray: A new numpy array with elements rounded to the specified significant digits, in string format.
    """

    def round_to_significant(x, sig_digits):
        if x == 0:
            return '0'
        else:
            # Calculate the number of decimal places needed
            from math import log10, floor
            precision = sig_digits - int(floor(log10(abs(x)))) - 1
            # Ensure precision is not negative
            precision = max(0, precision)
            # Format with scientific notation if necessary
            return f"{x:.{precision}e}" if precision < -3 or abs(x) >= 1e4 else f"{x:.{precision}f}"

    # Vectorize the rounding function to apply it to each element of the array
    vectorized_round = np.vectorize(lambda x: round_to_significant(x, sig_digits))

    return vectorized_round(arr)


def get_significant_digits(x):
    """
    Calculate the maximum number of significant digits after the decimal point
    in a list of floating-point numbers, including those in scientific notation.

    Args:
        x (list of float): A list of floating-point numbers, which may include
                            numbers in scientific notation.

    Returns:
        int: The maximum number of digits after the decimal point in the list.
             If no number has a decimal part, returns 0.

    Example:
        >>> get_significant_digits([-1e-06, -2e-06, -3e-06, 8e-05, -3e-05, 8e-06, 6e-06, 5e-05, 3e-05])
        6

    Notes:
        - The function uses the `decimal.Decimal` class to handle floating-point precision and scientific notation.
        - The input list can contain both small and large float values, and the function will correctly count significant digits for all valid inputs.
    """

    # Convert each number to a string
    x_str = [str(value) for value in x]

    # Use decimal module to handle floating point precision issues
    lenths = []
    for dig in x_str:
        # Try to parse scientific notation as a float, then convert to a string with fixed precision
        try:
            dec_value = decimal.Decimal(dig)
            # Split the number at the decimal point
            if '.' in str(dec_value):
                lenths.append(len(str(dec_value).split('.')[1]))
            else:
                lenths.append(0)
        except:
            lenths.append(0)

    # Return the maximum length of the decimal parts
    return max(lenths)


def merge_dataframes(df1, df2):
    """
    Merges two dataframes using an 'outer' join and fills missing values in common columns.

    The function performs an outer merge, retaining all indices from both dataframes. If columns
    have the same name in both dataframes, the function fills in any missing values ('None' or
    'NaN') by taking non-null values from the other dataframe where available.

    Args:
        df1 (pd.DataFrame): The first dataframe.
        df2 (pd.DataFrame): The second dataframe.

    Returns:
        pd.DataFrame: A new dataframe that results from the outer merge of df1 and df2, with
                      missing values in common columns filled in. """

    df1_cols = df1.columns
    df2_cols = df2.columns

    cols = list(df1_cols) + list(df2_cols)

    cols = list(set(cols))

    merged_df = pd.DataFrame(index=df1.index.union(df2.index), columns=cols)

    for col in cols:
        if col in df1.columns and col in df2.columns:
            merged_df[col] = df1[col].combine_first(df2[col])
        elif col in df1.columns:
            merged_df[col] = df1[col]
        elif col in df2.columns:
            merged_df[col] = df2[col]

    return merged_df


def find_string_with_substrings(strings, substrings):
    strings_out = []
    for string in strings:
        if any(substring in string for substring in substrings):
            strings_out.append(string)
    return strings_out


def filter_df_cols_by_keywords(df, keywords):
    """
    Filters columns of a dataframe based on a list of keywords.

    The function selects columns from the dataframe where the column names contain any of the keywords
    specified in the list. It returns a dataframe with only those columns.

    Args:
        df (pd.DataFrame): The dataframe to filter.
        keywords (list of str): A list of keywords to search for in column names.

    Returns:
        pd.DataFrame: A new dataframe containing only the columns with names that match any of the keywords.
    """

    # Create a boolean mask for columns containing any of the keywords
    mask = df.columns.to_series().str.contains('|'.join(keywords), case=False, na=False)

    # Filter the dataframe columns
    filtered_df = df.loc[:, mask]

    return filtered_df


def save_df_list_to_excel(path, dfs, sheet_names=None):
    """
    Save multiple DataFrames to an Excel file, with each DataFrame written to a separate sheet.

    Parameters:
    -----------
    path : str
        The file path where the Excel file will be saved (e.g., 'output.xlsx').

    dfs : list of pd.DataFrame
        A list of pandas DataFrames to save in the Excel file.

    sheet_names : list of str, optional
        A list of sheet names corresponding to each DataFrame.
        If not provided, default names 'Sheet1', 'Sheet2', etc., will be used.

    Raises:
    -------
    ValueError
        If the number of provided sheet names does not match the number of DataFrames.

    Returns:
    --------
    None
        The function saves the DataFrames to an Excel file and does not return anything.
    """

    path = path + '.xlsx'

    # If sheet names are not provided, generate default names like 'Sheet1', 'Sheet2', etc.
    if sheet_names is None:
        sheet_names = [f'Sheet{i + 1}' for i in range(len(dfs))]

    # Ensure that the number of DataFrames matches the number of sheet names
    if len(dfs) != len(sheet_names):
        raise ValueError("The number of DataFrames must match the number of sheet names")

    # Create a Pandas Excel writer object
    with pd.ExcelWriter(path, engine='xlsxwriter') as writer:
        # Iterate through DataFrames and corresponding sheet names
        for df, sheet in zip(dfs, sheet_names):
            # Write each DataFrame to a separate sheet
            df.to_excel(writer, sheet_name=sheet, index=True)


# Extreme Values
def gumbel_P(x, beta, mu):
    """"evaluates P = exp(-exp(-(x - mu) / beta)) for every x and constant beta and mu"""
    return np.exp(-np.exp(-(x - mu) / beta))


def gumbel_P_inv(P, beta, mu):
    """"evaluates x = - beta * np.log(-1 * np.log(P)) + mu for every P and constant beta and mu"""
    return - beta * np.log(-1 * np.log(P)) + mu


def gumbel_p(x, beta, mu):
    """"evaluates p = 1 / beta * np.exp(-1 / beta * (x - mu)) * np.exp(-1 * np.exp(-1 / beta * (x - mu))) for every x and constant beta and mu"""
    return 1 / beta * np.exp(-1 / beta * (x - mu)) * np.exp(-1 * np.exp(-1 / beta * (x - mu)))


def Intervall_conf_1(N, itter, beta, mu, perc_up, perc_down):
    P = np.empty((N, itter))
    X = np.empty((N, itter))

    for i in range(itter):
        P_samp = [random.uniform(0, 1) for i in range(N)]
        P_samp.sort()

        x_samp = gumbel_P_inv(P_samp, beta, mu)
        P[:, i] = P_samp
        X[:, i] = x_samp

    perc_up = np.percentile(X, perc_up, axis=1)
    perc_middle = np.percentile(X, 50, axis=1)
    perc_down = np.percentile(X, perc_down, axis=1)
    std_dev = np.std(X, axis=1)
    return perc_down, perc_middle, perc_up, std_dev


def Intervall_conf_2(N, itter, T_R_range, beta, mu, perc_up, perc_down):
    P = np.empty((N, itter))
    X = np.empty((N, itter))
    Q = np.empty((itter, len(T_R_range)))
    for i in range(itter):
        P_samp = [random.uniform(0, 1) for i in range(N)]
        P_samp.sort()

        x_samp = gumbel_P_inv(P_samp, beta, mu)
        P[:, i] = P_samp
        X[:, i] = x_samp

    Omega = [Gumbel_coeff(row) for row in X.T]

    F = T_R_to_F(T_R_range, 1)

    i = 0
    for Omega_curr in Omega:
        Q[i, :] = [gumbel_P_inv(F_temp, Omega_curr[0], Omega_curr[1]) for F_temp in F]
        i = i + 1

    perc_up = np.percentile(Q, perc_up, axis=0)
    perc_middle = np.percentile(Q, 50, axis=0)
    perc_down = np.percentile(Q, perc_down, axis=0)
    std_dev = np.std(Q, axis=0)
    return perc_down, perc_middle, perc_up, std_dev


def Gumbel_coeff(x):
    """"takes array, list or series of maximal values and calculates the gumbel coefficients beta and mu

    theorie:
    P_gumb = exp(-exp(-(x - mu) / beta))

    beta = beta - sigma * (sqrt(6)/pi)
    mu = mean(x) - gamma * beta
    gamma =  0.57721566490153286060651209008240243104215933593992
    """

    x = np.array(x)
    std = np.std(x)
    beta = std * np.sqrt(6) / np.pi
    gamma = 0.57721566490153286060651209008240243104215933593992
    mu = np.mean(x) - gamma * beta

    return beta, mu


def Hs_to_T_R(Hs, beta, mu, freq_samp):
    F = gumbel_P(Hs, beta, mu)
    T_R = [1 / ((1 - F_curr) * freq_samp) if F_curr != 1 else float("nan") for F_curr in F]

    return T_R


def T_R_to_F(T_R, freq_samp):
    F = [(1 - 1 / (T_R_temp * freq_samp)) for T_R_temp in T_R]
    return F


def x_max_sampeling(x, time_window_offeset, n_samp=None):
    """takes a series with values and datetime objects as indizes. Returns the n_samp biggest values in a one year period, ofsetted by time_window_offset

    arguments:
    x: Series with values and datetime objects as indizes
    time_window_offeset: float [0..1]: offset of the time window (relative)

    optional:
    n_samp: int, number if maximal values per year

    retruns:
    x_max: Series with maximal values and datetime objects as indizes
    years_newyear: fist datetime object of every year
    years_start: first datetime object of every evaluated timewindow (offset by time_window_offeset)
    """

    # todo: include higher/lower sampling rate than one time per year

    years = np.unique(x.index.year)
    window = pd.Timedelta(time_window_offeset * 365, "d")

    years_newyear = [datetime.date(year, 1, 1) for year in years]

    years_start = [new_year + window for new_year in years_newyear]

    x_max = pd.Series()

    for i in range(len(years_start) - 1):
        x_section = x.loc[years_start[i]:years_start[i + 1]]

        x_max_temp = max(x_section)
        indx_x_max = np.array(x_section).argmax()
        time_x_max = x_section.index[indx_x_max]

        x_max.at[time_x_max] = x_max_temp

    return x_max, years_newyear, years_start


def gumbel_conf_intervall(x_max, beta=None, mu=None, mode='percentile', algorithm='1', N_itter=1000, freq_samp=1, perc_up=95, perc_down=5, T_max=100):
    N_xmax = len(x_max)
    if beta is None or mu is None:
        beta, mu = Gumbel_coeff(x_max)

    if algorithm == '1':

        F_real = [(n + 0.5) / N_xmax for n in range(N_xmax)]
        T_R_grid = [1 / ((1 - F_curr) * freq_samp) if F_curr != 1 else float("nan") for F_curr in F_real]

        if mode == 'percentile':
            band_down, middle, band_up, std_dev = Intervall_conf_1(N_xmax, N_itter, beta, mu, perc_up, perc_down)

        if mode == 'std':
            _, middle, _, std_dev = Intervall_conf_1(N_xmax, N_itter, beta, mu, perc_up, perc_down)
            band_up = [middle[i] + std_dev_curr for i, std_dev_curr in enumerate(std_dev)]
            band_down = [middle[i] - std_dev_curr for i, std_dev_curr in enumerate(std_dev)]

    if algorithm == '2':

        T_R_grid = np.logspace(0, np.log10(T_max), 200)

        if mode == 'percentile':
            band_down, middle, band_up, std_dev = Intervall_conf_2(N_xmax, N_itter, T_R_grid, beta, mu, perc_up, perc_down)

        if mode == 'std':
            _, middle, _, std_dev = Intervall_conf_2(N_xmax, N_itter, T_R_grid, beta, mu, perc_up, perc_down)
            band_up = [middle[i] + std_dev_curr for i, std_dev_curr in enumerate(std_dev)]
            band_down = [middle[i] - std_dev_curr for i, std_dev_curr in enumerate(std_dev)]

    return band_down, middle, band_up, T_R_grid


# Function to check if df1 is in df2, and if so, output the row number of the matching row

def add_unique_row(df1, df2, exclude_columns=None):
    """
    Checks if a row from df1 exists in df2 (excluding specified columns).
    If the row exists, returns the updated df2 and the indices of the matching row(s).
    If the row does not exist, appends the row from df1 to df2 and returns the
    updated df2 with an empty list of matching indices.

    Parameters:
    df1 (pd.DataFrame): A DataFrame with a single row that will be checked
                        against df2.
    df2 (pd.DataFrame): A DataFrame that may contain one or more rows, which
                        will be compared to the row in df1.
    exclude_columns (list, optional): List of columns to exclude from the
                                      uniqueness check. Default is None.

    Returns:
    tuple:
        pd.DataFrame: The updated DataFrame (df2), either with or without the
                      new row from df1.
        list: A list of indices where the row from df1 matches any row in df2.
              If no match is found, the list is empty.
    """
    if exclude_columns is None:
        exclude_columns = []

    # Drop excluded columns from both df1 and df2 for comparison
    df1_comp = df1.drop(columns=exclude_columns, errors='ignore')
    df2_comp = df2.drop(columns=exclude_columns, errors='ignore')

    # Check for rows in df2 that match the row in df1
    matching_rows = df2_comp[df2_comp.eq(df1_comp.values[0]).all(axis=1)]

    if not matching_rows.empty:
        # If a match is found, get the row indices
        matching_indices = matching_rows.index.tolist()
    else:
        # If no match is found, append df1 to df2
        df2 = pd.concat([df2, df1], ignore_index=True)
        matching_indices = []

    return df2, matching_indices


# Validation
def calc_JBOOST(path_exe, proj_name, Hs, Tp, gamma):
    # exporting JBOOST input Files of Database
    write_JBOOST_wave(Hs, Tp, gamma, path_exe + 'wave.lua')

    # Run JBOOST
    subprocess.check_call(['JBOOST.exe', proj_name],
                          cwd=path_exe, shell=True)

    # import points of Database
    temp, debug = import_JBOOST(path_exe + 'Results_JBOOST_Text/JBOOST.out')

    out = {}
    for ck, cd in temp.items():
        out[ck] = pd.DataFrame(data=cd.values, index=Hs.index, columns=cd.columns)

    # write all node in one dataframe
    DEL_Data_full = pd.DataFrame()
    for ck, cd in temp.items():
        cd = cd.rename(columns={'IDLING': 'IDLING ' + ck, 'PRODUCTION': 'PRODUCTION ' + ck})

        DEL_Data_full = (DEL_Data_full.merge(cd, left_index=True, right_index=True, how='outer'))

    DEL_Data_full.index = Hs.index

    return DEL_Data_full


def import_JBOOST(path):
    lookup = 'Result Hindcast'

    ID = []
    NODE = []
    IDLING = []
    PRODUCTION = []

    with open(path) as myFile:
        for num, line in enumerate(myFile, 1):
            if lookup in line:
                line_num = num
                break

        for num_2, line_2 in enumerate(myFile, line_num):
            if (num_2 > line_num + 1) & (line_2 != '\n'):
                out = line_2.strip().split("\t")
                ID.append(int(out[0]))
                NODE.append(int(out[1]))
                IDLING.append(float(out[2]))
                PRODUCTION.append(float(out[3]))

    df = pd.DataFrame()
    ID_series = pd.Series(ID)
    Node_series = pd.Series(NODE)
    df['IDLING'] = IDLING
    df['PRODUCTION'] = PRODUCTION

    Nodes = Node_series.unique()

    JBOOST_OUT = {}
    for node in Nodes:
        temp = df[Node_series == node]
        temp.index = ID_series[Node_series == node]
        JBOOST_OUT['node ' + str(node)] = temp

    return JBOOST_OUT, IDLING


def write_JBOOST_wave(Hs, Tp, gamma, path):
    if isinstance(gamma, float):
        gamma = [gamma for i in range(Hs.shape[0])]

    i = 0
    with open(path, "w") as text_file:
        for Hs_curr, Tp_curr, gamma_curr in zip(Hs, Tp, gamma):
            string = str()
            string += 'os_Hindcast{'

            string += f'id = {i}, '
            string += f'Hs = {Hs_curr:.2E}, '
            string += f'Tp = {Tp_curr:.2E}, '
            string += f'gamma = {gamma_curr:.2E}'

            string += '}\n'

            text_file.write(string)

            i = i + 1
    return


def write_DEL_base(path_DataBase, DEL_data, Meta_Data):
    # write meta data in dataframe
    df_Meta = pd.DataFrame(columns=list(Meta_Data.keys()))

    for col in df_Meta.columns:
        df_Meta.loc[0, col] = Meta_Data[col]

    # connect to db
    conn = sqlite3.connect(path_DataBase)
    cursor = conn.cursor()

    query = "SELECT name FROM sqlite_master WHERE type='table';"

    # Execute the query
    cursor.execute(query)

    # Fetch all table names
    table_names = cursor.fetchall()
    table_names = [name[0] for name in table_names]

    # If the table doesn't exist, create it
    if 'DEL_Meta' not in table_names:
        # Use the DataFrame's to_sql() method to create the table
        DEL_config_name = 'DEL_config_1'
        df_Meta = df_Meta.set_index(pd.Index([DEL_config_name]))
        df_Meta.to_sql('DEL_Meta', conn, if_exists='fail', index=True)

    else:
        df_meta_sql = pd.read_sql(f"SELECT * FROM 'DEL_Meta'", conn, index_col='index')

        duplicate_name = check_meta_in_valid_db(path_DataBase, Meta_Data)
        if len(duplicate_name) > 0:
            DEL_config_name = duplicate_name[0]
        else:
            DEL_config_name = f'DEL_config_{len(df_meta_sql) + 1}'
            df_Meta.set_index(pd.Index([DEL_config_name]))
            df_Meta = df_Meta.set_index(pd.Index([DEL_config_name]))
            df_meta_combined = pd.concat((df_meta_sql, df_Meta), axis=0)

            df_meta_combined.to_sql('DEL_Meta', conn, if_exists='replace', index=True)

    DEL_data.to_sql(DEL_config_name, conn, if_exists='replace', index=True)

    conn.close()

    return DEL_config_name


def check_meta_in_valid_db(db_path, Meta, exclude_columns=None):
    # write meta data in dataframe
    df_Meta = pd.DataFrame(columns=list(Meta.keys()))

    for col in df_Meta.columns:
        df_Meta.loc[0, col] = Meta[col]

    # connect to db
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # check if DEL_Meta tabel exists
    table_name = 'DEL_Meta'
    cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name=?", (table_name,))

    meta_exists = cursor.fetchone()

    if meta_exists:
        # check, if calculation is already been run in the database
        df_Meta_sql = pd.read_sql_query(f"SELECT * FROM {table_name}", conn, index_col='index')
        _, idx_in_meta = add_unique_row(df_Meta, df_Meta_sql, exclude_columns=exclude_columns)

    else:
        idx_in_meta = []

    return idx_in_meta


import pandas as pd


def median_sample_rate(index) -> pd.Timedelta:
    """
    Calculate the median sample rate of a series of datetime objects.

    The sample rate is determined by the median of the time differences
    between consecutive datetime values.

    Parameters:
    index (Union[pd.DatetimeIndex, List[pd.Timestamp]]): A pandas DatetimeIndex or a list of datetime objects.

    Returns:
    pd.Timedelta: The median sample rate as a Timedelta object.
    """
    # Convert list of datetime objects to DatetimeIndex if necessary
    if isinstance(index, list):
        index = pd.DatetimeIndex(index)

    # Calculate the time differences between consecutive datetime values
    time_diffs = pd.Series(index).diff().dropna()

    # Calculate the median sample rate and convert it to a Timedelta
    median_diff = time_diffs.median()

    return median_diff


def fill_nans_constant(x_vector, mask=None):
    """
    Fills NaN values in x_vector using the last non-NaN value,
    using a boolean mask to determine where to fill.
    If no mask is provided, it assumes all values are to be filled.

    Parameters:
    x_vector (np.ndarray): The input vector containing NaNs.
    mask (np.ndarray, optional): A boolean mask of the same length as x_vector.
                                 Defaults to an array of True values.

    Returns:
    np.ndarray: The x_vector with NaNs filled according to the mask.

    Notes:
    - If the x_vector starts with NaNs and the mask is all True, those NaNs will be filled with the fist not nan value unfilled.
    - If the mask is not provided, all NaNs will be filled with the last valid value encountered.
    """
    if mask is None:
        mask = np.ones_like(x_vector, dtype=bool)  # Create a mask of all True

    filled_vector = np.copy(x_vector)
    last_valid = None
    leading_nan = True
    i_leading_nan = 0

    for i in range(len(filled_vector)):
        if mask[i]:  # Only fill if the mask is True
            if np.isnan(filled_vector[i]):
                if leading_nan:
                    i_leading_nan += 1

                if last_valid is not None:
                    filled_vector[i] = last_valid
            else:
                last_valid = filled_vector[i]  # Update last valid value
                if leading_nan:
                    filled_vector[i - i_leading_nan:i] = last_valid
                    leading_nan = False


        else:
            # If mask is False, do not fill and reset last_valid
            last_valid = None if np.isnan(filled_vector[i]) else filled_vector[i]

    return filled_vector


def fill_nan_with_linspace(vector):
    # Convert the input to a NumPy array if it's not already one
    vector = np.asarray(vector)

    # Find the indices of non-NaN values
    non_nan_indices = np.where(~np.isnan(vector))[0]

    # Check if there are at least two non-NaN values
    if len(non_nan_indices) < 2:
        raise ValueError("At least two non-NaN values are required to perform linear interpolation.")

    # Get the first and last non-NaN values
    start_value = vector[non_nan_indices[0]]
    end_value = vector[non_nan_indices[-1]]

    # Create an array of the same length as the original
    filled_vector = np.copy(vector)

    # Generate linearly spaced values
    linspace_values = np.linspace(start_value, end_value, num=len(vector))

    # Fill in the NaN values
    filled_vector[np.isnan(vector)] = linspace_values[np.isnan(vector)]

    return filled_vector


def xlsx2csv(excel_file, output_dir, exclude_sheets=None):
    """
    Save each sheet of an Excel file as a CSV file in the specified output directory.

    Parameters:
    excel_file (str): Path to the input Excel file.
    output_dir (str): Directory where CSV files will be saved. Created if it doesn't exist.
    exclude_sheets (list of int, optional): List of sheet indices (1-based) to exclude.
                                            For example, [1, 3] will exclude the first and third sheets.

    Returns:
    None
    """
    # Load the Excel file
    excel_data = pd.ExcelFile(excel_file)

    # Ensure output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Normalize exclude_sheets to 0-based indexing
    if exclude_sheets:
        exclude_sheets = [i - 1 for i in exclude_sheets]

    # Loop through each sheet in the Excel file
    for sheet_index, sheet_name in enumerate(excel_data.sheet_names):
        if exclude_sheets and sheet_index in exclude_sheets:
            print(f"Skipping sheet {sheet_name} (index {sheet_index + 1})")
            continue

        # Read the sheet into a DataFrame
        df = pd.read_excel(excel_file, sheet_name=sheet_name)

        # Define the output CSV file path
        csv_file = os.path.join(output_dir, f"{sheet_name}.csv")

        # Save the DataFrame to a CSV file
        df.to_csv(csv_file, index=False)
        print(f"Saved {sheet_name} to {csv_file}")


def e57_2_txt(file_path, save_path =None):
    """
    Converts E57 point cloud files to TXT format by extracting and saving both point color
    and point coordinate data.

    Parameters:
    -----------
    file_path : str
        The path to the directory containing E57 files or a specific E57 file.
    save_path : str, optional
        The directory where the TXT files will be saved. If not provided, files are saved
        in the same directory as the input E57 files.

    Returns:
    --------
    None

    Notes:
    ------
    - Only processes files with the '.e57' extension found in the specified directory.
    - For each E57 file:
        - The color data of the points is scaled to a range of 0-255, rounded,
          and saved in a file named `<basename>_color.txt`.
        - The point coordinate data is saved in a file named `<basename>_points.txt`
          with four decimal places.
    - Both TXT files are saved in the specified or default output directory.
    - Prints progress messages for each file being processed.
    """