import numpy as np
import scipy as sc
import pandas as pd
import sys
import sqlite3
from sqlite3 import Error
import os
import shutil
import chardet


from libaries import general as gl


# %%classes
class Segment:
    def __init__(self, num, angles, angle_name, result, indizes=None, colnames=None):
        self.num = num
        self.angles = angles
        self.result = result
        self.indizes = indizes
        if angles is not None:
            self.angles_mod = [angles[0] % 360, angles[1] % 360]
        else:
            self.angles_mod = None

        self.colnames = colnames
        if colnames is None:
            self.colnames = dict()

        if indizes is not None:
            self.N = len(indizes)
        else:
            self.N = None
        self.angle_name = angle_name

    def FilterData(self, Data):
        Data = Data.loc[(Data[self.angle_name] >= self.angles_mod[0])
                        | (Data[self.angle_name] < self.angles_mod[1])]
        return Data


class Calculation:

    def __init__(self, result=None, basedata=None, filt=None):
        self.result = result

        self.basedata = basedata
        if basedata is None:
            self.basedata = dict()

        self.filt = filt

    def initilize_from_db(self, db_path, table_name, colnames, **kwargs):
        """initilizes dataframe from sql database at "db_path" by loading the colims specified by colnames in the table "table_name". If colnames is None, all columns are loaded
         Returns dataframe

         optional:
         timeframe: list of two datetime objects specifying the start and end time of the data to be evaluated
         """

        timeframe = kwargs.get('timeframe', None)
        indizes = kwargs.get('indizes', None)

        df = gl.export_df_from_sql(db_path, table_name, column_names=colnames)
        df.index = pd.to_datetime(df.index)

        if timeframe is not None:
            df = df.loc[timeframe[0]:timeframe[1]]

        self.basedata["dbname"] = db_path
        self.basedata["tablename"] = table_name
        self.basedata["colnames_ini"] = colnames
        self.basedata["dbdate"] = None
        self.basedata["db_timeframe"] = [df.index[0], df.index[-1]]
        self.basedata["N_rows"] = len(df)
        self.basedata["sample_rate"] = gl.median_sample_rate(df.index)

        if indizes is None:
            self.basedata["indizes"] = df.index
        else:
            self.basedata["indizes"] = indizes

        return df

    def initilize_filter(self, colnames, mode='range', ranges=None):
        """filters data specified in "basedata" propterty, has to be initilized by using "initilize_from_db". creats "filt" dictionary, that contains:
         indizes_in: datetime index object, with all indizies still in the filtered list
         indizes_out: datetime index object, with all indizies excluded by filtering
         colnames: colnames
         ranges: ranges

        INPUT:
        colnames: list of stings, columnames to filter by in Dataframe
        ranges: list of lists, list containing len(colnames) elements of length 2, in fomrmat [a_min, a_max] specifing the range of the filtered data
                if a_max is None, maximal value of the correspnding values of column spcified by colname in with same index used

        retrun:
        indizes: datetime index object, with all indizies still in the filtered list"""

        indizes_full = self.basedata["indizes"]

        if colnames is None:
            df_filt = gl.export_df_from_sql(self.basedata["dbname"], self.basedata["tablename"])
            colnames = list(df_filt.columns)
        else:
            df_filt = gl.export_df_from_sql(self.basedata["dbname"], self.basedata["tablename"], column_names=colnames)

        self.filt = dict()

        if mode == 'range':
            for i, range_curr in enumerate(ranges):
                if range_curr[1] is None:
                    range_curr[1] = np.max(df_filt[colnames[i]])

            a_min = [curr[0] for curr in ranges]
            a_max = [curr[1] for curr in ranges]

            df_filt = gl.filter_dataframe(df_filt, colnames, a_min, a_max)
            self.filt["ranges"] = ranges
            self.filt["colnames"] = colnames

        if mode == 'nans':
            df_filt = df_filt.dropna(how='any')

        index_filt = df_filt.index
        indizes_in = index_filt.intersection(indizes_full)
        indizes_out = indizes_full.difference(index_filt)

        self.filt["indizes_in"] = indizes_in
        self.filt["indizes_out"] = indizes_out
        self.filt["mode"] = mode

        return indizes_in
    
    def create_segment_title(self, mode='verbose', latex=True):
        """"if the data stored in result is list of "Segment",  it exports a list of title in Latex format
            if mode='general', Segment.basedata has to be initilized
        """""

        titles = []
        for segment_curr in self.result:
            underscore = str()
            header = str()
            N_exp = segment_curr.N

            if segment_curr.angle_name is not None:
                if latex:
                    header += r"\small " f"'{segment_curr.angle_name}': {segment_curr.angles[0]}° to {segment_curr.angles[1]}°"
                else:
                    header += f"'{segment_curr.angle_name}':{segment_curr.angles[0]}° to {segment_curr.angles[1]}°"

            else:
                header += "omnidirectional"

            if mode == "verbose":
                N_ges = self.basedata["N_rows"]
                timeframe = self.basedata["db_timeframe"]
                sample_rate = self.basedata["sample_rate"]
                if latex:
                    underscore = (r"\scriptsize " + f"samples: {N_exp:.2e} ({round(N_exp / N_ges * 100, 1)}\%), " +
                                  f"{timeframe[0].round('1d').date()} to {timeframe[1].round('1d').date()}, " +
                                  f"d_t: {sample_rate.total_seconds()} s")
                else:
                    underscore = (f"samples: {N_exp:.2e} ({round(N_exp / N_ges * 100, 1)}%), " +
                                  f"timeframe: {timeframe[0].round('1d').date()} to {timeframe[1].round('1d').date()}, " +
                                  f"time step: {sample_rate.total_seconds()} s")
            elif mode == "standard":
                N_ges = self.basedata["N_rows"]

                if latex:
                    underscore = r"\scriptsize " + f"samples: {N_exp:.2e} ({round(N_exp / N_ges * 100, 1)}%)"
                else:
                    underscore = f"samples: {N_exp:.2e} ({round(N_exp / N_ges * 100, 1)}%)"

            elif mode == "sparse":
                underscore = None

            else:
                underscore = None

            if underscore is not None:
                titles.append(header + "\n" + underscore)
            else:
                titles.append(header)
        return titles

    def load_from_db(self, column_names=None, applie_filt=True, colnames_ini=False, **kwargs):
        """wrapper for export_df_from_sql in general lib, takes db_name and tablename from Calculation information, applies filter if it is there"""
        if colnames_ini:
            column_names = self.basedata["colnames_ini"]

        df = gl.export_df_from_sql(self.basedata["dbname"], self.basedata["tablename"], column_names=column_names)

        if (self.filt is None) or not applie_filt:
            df = df[df.index.isin(self.basedata["indizes"])]
        else:
            df = df[df.index.isin(self.filt["indizes_in"])]

        return df


class DataCol:
    def __init__(self, name_data=str(), name_plot=str(), db_name=None, table_raw=None, symbol=None):
        self.name_data = name_data
        self.name_plot = name_plot
        self.db_name = db_name

        if table_raw is None:
            self.table_raw = []

        self.symbol = symbol


# %% general functions

def percentiles(df, percent: list):
    for n_perc, perc in enumerate(percent):

        df[f'{perc}th percentile'] = np.percentile(df['mean'], perc)
    return


def condensation(x, y, grid,
                 reg_model='poly',
                 deg_reg=3,
                 cut_reg=0,
                 reg_weighting=0,
                 zone_reg=None,
                 zone_line=None,
                 perc=None,
                 bin_min=0,
                 average_correction=1.0,
                 avrg_method='mean',
                 make_monotone=False):

    # Set default values if None is passed
    zone_reg = zone_reg if zone_reg is not None else [None, None]
    zone_line = zone_line if zone_line is not None else [None, None]
    perc = perc if perc is not None else []

    if zone_reg[1] is None:
        zone_reg[1] = max(x)

    if zone_reg[0] is None:
        zone_reg[0] = min(x)

    if zone_line[1] is None:
        zone_line[1] = max(x)

    if zone_line[0] is None:
        zone_line[0] = min(x)

    # find x-zone
    n_bin = len(grid) - 1
    x_zone = [min(grid), max(grid)]

    # avearge points (only middle)
    averaged, std, count, bin_ident = gl.grid_pointcloud_in_x(
        x, y, grid, method=avrg_method)

    averaged.name = 'averaged'

    x_bins = averaged.index

    # percentiles
    percentiles = []
    for perc_curr in perc:
        perc_data_curr = pd.Series(index=x_bins, name=f"{perc_curr} percentile")
        for curr_ident in np.unique(bin_ident):
            perc_data_curr.iloc[curr_ident-1] = np.percentile(y.iloc[np.where(bin_ident == curr_ident)], perc_curr)

        percentiles.append(perc_data_curr)

    # average correction
    averaged= averaged * average_correction

    for perc_curr in percentiles:
        perc_curr = perc_curr * average_correction

    #zone plot

    plot_line = (x_bins > zone_line[0]) & (
            x_bins < zone_line[1])

    # make monoton
    if make_monotone:
        averaged.loc[plot_line & ~np.isnan(averaged.values)] = gl.make_monotone(averaged.loc[plot_line & ~np.isnan(averaged.values)])

        for perc_curr in percentiles:
            perc_curr.loc[plot_line & ~np.isnan(perc_curr.values)] = gl.make_monotone(perc_curr.loc[plot_line & ~np.isnan(perc_curr.values)])

    nanMask = count > bin_min

    # cut_reg bereich
    use_regression = pd.Series(index=x_bins, name='use_regression', dtype=bool)
    use_regression[:] = False

    N_upper = round(cut_reg / 100 * len(x))
    x_points_sorted = np.sort(x)

    if N_upper == 0:
        use_regression[:] = True

    elif N_upper != len(x):
        x_lim_upper = x_points_sorted[N_upper]

        # find bin in which the regression cut is and set use_regression to 1 wenn regressionosbereich
        _, edges, vs_bin = sc.stats.binned_statistic(
            x_lim_upper, x_lim_upper, statistic='count', bins=n_bin, range=x_zone)

        use_regression.iloc[vs_bin[0]:] = True

    else:
        use_regression.iloc[:] = False
    # regression
    # regressionsbereich

    reg_zone = (x_bins > zone_reg[0]) & (x_bins < zone_reg[1])


    combined = []
    regression = []
    combined_plot = []
    Coeffs = {}
    lines = [averaged.copy()] + [perc.copy() for perc in percentiles]
    for line_curr in lines:

        if any(use_regression):

            x_reg_zone = line_curr.index[reg_zone & nanMask]
            y_reg_zone = line_curr.values[reg_zone & nanMask]
            counts_curr = count.values[reg_zone & nanMask]

            reg_model, coeffs = gl.model_regression(x_reg_zone, y_reg_zone, degree=deg_reg, weights=counts_curr, weights_regulation=reg_weighting, reg_model=reg_model)

            line_curr_regression = gl.predict_regression(reg_model, x_bins)

            line_curr_regression = pd.Series(data=line_curr_regression, index=x_bins, name=f"{line_curr.name} regression")

        else:
            line_curr_regression = pd.Series(index=x_bins, name=f"{line_curr.name} regression")
            line_curr_regression[:] = float('nan')
            coeffs = None

        Coeffs[line_curr.name] = coeffs

        # combine averaged with regression
        line_curr_combined = line_curr.copy()
        line_curr_combined.loc[use_regression] = line_curr_regression[use_regression]
        line_curr_combined.name = f"{line_curr.name} result"

        # clip to plot zone
        line_curr_combined_plot = line_curr_combined.copy()
        line_curr_combined_plot.name = f"{line_curr.name} plot"
        line_curr_combined_plot[~plot_line] = float("nan")

        combined.append(line_curr_combined)
        regression.append(line_curr_regression)
        combined_plot.append(line_curr_combined_plot)

    OUT = pd.DataFrame()
    OUT["x"] = x_bins
    OUT["mean"] = averaged.values
    OUT["std"] = std.values
    OUT["count"] = count.values

    for perc_curr in percentiles:
        OUT[perc_curr.name] = perc_curr.values

    OUT["mean regression"] = regression[0].values
    OUT["mean result"] = combined[0].values
    OUT["mean result plot"] = combined_plot[0].values

    if len(combined) > 1:
        for perc_combined, perc_combined_plot in zip(combined[1:], combined_plot[1:]):
            OUT[perc_combined.name] = perc_combined.values
            OUT[perc_combined_plot.name] = perc_combined_plot.values
        for perc_regression in regression[1:]:
            OUT[perc_regression.name] = perc_regression.values

    OUT['isData'] = nanMask.values.astype(int)
    OUT['use_regression'] = use_regression.values.astype(int)
    OUT['bool_reg_zone'] = reg_zone.astype(int)
    OUT['bool_plot_zone'] = plot_line.astype(int)
    OUT.loc[OUT['count'] <= bin_min, 'isData'] = 0

    return OUT, Coeffs


def quantiles(perc_low, middle, perc_up, quant_low, quant_up):
    """" expects 2 percentiles and HSTP Data """
    quantile = np.empty(len(middle))
    quantile[:] = float('nan')

    try:
        T_up = 1 / quant_up
        T_low = 1 / quant_low

        bool_LOW = (perc_up < T_low)
        bool_MIDDLE = (middle >= T_low) & (middle <= T_up)

        bool_UP = (perc_low > T_up)

        bool_CONNECT_LOW = ~bool_LOW & ~bool_MIDDLE & ~bool_UP & (middle < (T_low + T_up) / 2)
        bool_CONNECT_HIGHT = ~bool_LOW & ~bool_MIDDLE & ~bool_UP & (middle > (T_low + T_up) / 2)

        quantile[bool_CONNECT_LOW] = T_low
        quantile[bool_CONNECT_HIGHT] = T_up

        quantile[bool_LOW] = perc_up[bool_LOW]

        quantile[bool_MIDDLE] = middle[bool_MIDDLE]
        quantile[bool_UP] = perc_low[bool_UP]

    except:
        print(
            f"    quantile not possible for segment, check if graph is monotone in 'zone_line' or if percentiles cross in the frequency band. quantile is set to mean")

    return quantile


def csv_to_sqlDB(path_csvs, db_name, resample_rate, data_kind='MetOcean', encoding='auto', nans=None, skiprows=None, delimiter=';', dayfirst=False, datetime_mode='single_col', low_memory=True, drop_rows=None):
    """accepts Dict of csv paths (Metocean Format!) and stores the individual data as well a resampled version in
    a combined table to a sql database. If the database already exists, it will be overwritten

    Parameters:
        Paths: dict with the keys: {path_wave, path_atmo, path_ocean} and the corresponding paths to the csv files, can be None
        db_name: sting to name the database, with path if necessary
        resample_rate: resample rate in the fomrat: float {y,d,s,m} for year, day, second, month. Example: "1d"

    Return:
        sql Database at the desired path
    """

    def join_dataframes_by_index(dataframes):
        # Perform the join operation on the non-empty dataframes
        res = dataframes[0]
        for df in dataframes[1:]:
            res = res.join(df, how='inner')
        return res

    def read_MetaData_Metocean(file):
        global name
        Names = []
        Values = []
        with open(file) as input_file:
            for _ in range(14):
                line = input_file.readline()
                line = line.replace('\n', '')
                line = line.replace('"', '')
                splited = line.split('\t')
                name = splited[0]
                value = splited[-1][2:]
                Names.append(name)
                Values.append(value)

        return Names, Values

    def check_encoding(path, sample_size=10000):
        with open(path, 'rb') as file:
            # Read a sample of the file
            raw_data = file.read(sample_size)
            # Detect the encoding
            result = chardet.detect(raw_data)
            encoding = result
            return encoding

    data_resampled = []
    NAMES = []
    VALUES = []
    csv_files = [os.path.join(path_csvs, f) for f in os.listdir(path_csvs) if f.endswith(('.csv', '.CSV'))]

    db_path = path_csvs + "/" + db_name
    db_exists = os.path.exists(db_path)

    if db_exists:
        overwrite = input(f"    Database {db_path} already exists, overwrite? (y/n)")
        if overwrite == 'y':
            os.remove(db_path)
        else:
            return None

    conn = sqlite3.connect(db_path)

    for file in csv_files:
        
        if encoding == 'auto':
            encoding = check_encoding(file)
    
        file_name = os.path.basename(file)
        print(f"adding {file_name} to database")

        if data_kind == 'MetOcean':
            # values
            skiprows = 15
            df_NAN = pd.read_csv(file, skiprows=skiprows, na_filter=False, encoding=encoding)

            # metadata
            meta = read_MetaData_Metocean(file)

        elif data_kind == 'APGMer':
            # values
            df_NAN = pd.read_csv(file, sep=';', low_memory=False, na_filter=False, encoding=encoding)
            df_NAN = df_NAN.drop(0)
            df_datetime = df_NAN[df_NAN.columns[[0, 1, 2, 3, 4]]]

            df_NAN.drop(df_NAN.columns[[0, 1, 2, 3, 4]], axis=1, inplace=True)

            # metadata
            meta = [[], []]
        
        elif data_kind is None:
            df_NAN = pd.read_csv(file, skiprows=skiprows, encoding=encoding, delimiter=delimiter, na_values=nans, low_memory=low_memory)

            if datetime_mode == 'single_col':
                df_datetime = df_NAN.loc[:, df_NAN.columns[0]]
                df_NAN = df_NAN.drop(df_NAN.index[drop_rows])
                df_datetime = df_datetime.drop(df_datetime.index[drop_rows])
                df_NAN.drop(df_NAN.columns[[0]], axis=1, inplace=True)
                df_datetime = pd.to_datetime(df_datetime, dayfirst=dayfirst)
                
            elif datetime_mode == 'multi_col':
                df_datetime = df_NAN[df_NAN.columns[[0, 1, 2, 3, 4]]]

                df_NAN.drop(df_NAN.columns[[0, 1, 2, 3, 4]], axis=1, inplace=True)
                df_NAN = df_NAN.drop(drop_rows)
                df_datetime = df_datetime.drop(drop_rows)
                df_datetime = pd.to_datetime(df_datetime, dayfirst=dayfirst)

            else:
                print("please choose from 'single_col' or 'multi_col'")

            df_NAN = df_NAN.astype(float)
            meta = [[], []]
        else:
            print('   please choose "APGMer","MetOcean" or None as DataBase mode')
            return None

        df_NAN.index = df_datetime
        df = df_NAN.dropna(how='all')
        meta[0].append("NANs")
        meta[1].append(len(df_NAN) - len(df))
        meta[0].insert(0, "DataSet")
        meta[1].insert(0, file)
        NAMES.append(meta[0])
        VALUES.append(meta[1])

        df.to_sql('Hindcast_raw_' + file_name, conn)

        df_resample = df.resample(resample_rate).mean()

        data_resampled.append(df_resample)

    df_ges = join_dataframes_by_index(data_resampled)

    df_ges.to_sql('Hindcast_combined', conn)

    # Metadata
    df_Meta = pd.DataFrame(columns=NAMES[0][1:], index=[values[0] for values in VALUES])

    for values in VALUES:
        df_Meta.loc[values[0], :] = values[1:]

    df_Meta.to_sql('MetaData', conn)

    conn.close()

    return db_path


def extract_data(db_file, table_name, column_names=None):
    """
    Load specified columns or all columns from a table in an SQLite database and return as a pandas DataFrame,
    if tables index is specified, it is set as the index of the dataframe

    Parameters:
    - db_file (str): Path to the SQLite database file.
    - table_name (str): Name of the table to query.
    - column_names (list of str, optional): List of column names to retrieve. If None, retrieves all columns.

    Returns:
    - pandas.DataFrame: DataFrame containing the query results, including the row index as a column.

    Raises:
    - ValueError: If the table or columns do not exist.
    """
    try:
        conn = sqlite3.connect(db_file)
        cursor = conn.cursor()

        # Check if table exists
        cursor.execute(f"SELECT name FROM sqlite_master WHERE type='table' AND name='{table_name}';")
        if not cursor.fetchone():
            raise ValueError(f"Table '{table_name}' does not exist in the database.")

        if column_names:
            # Check if specified columns exist in the table
            cursor.execute(f"PRAGMA table_info({table_name});")
            columns_info = cursor.fetchall()
            table_columns = [info[1] for info in columns_info]
            for col in column_names:
                if col not in table_columns:
                    raise ValueError(f"Column '{col}' does not exist in table '{table_name}'.")
            # Enclose column names with special characters in double quotes
            cols = ', '.join([f'"{col}"' for col in column_names])
        else:
            # If no columns are specified, select all columns
            cols = '*'

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

            cursor.execute(f'SELECT "{index_col_name}" FROM "{table_name}";')
            index_col = cursor.fetchall()

            index_col = [info[0] for info in index_col]

            df.index = index_col

        return df

    except Error as e:
        print(f"An error occurred: {e}")
        raise
    finally:
        if conn:
            conn.close()


def RWI(Hs, Tp, f_0, gamma_mode='default'):
    """calculates the RWI, using the JONSWAP sepectrum for Hs/Tp lists, arrays or Series"""
    RWI_list = []
    for Hs_curr, Tp_curr in zip(Hs, Tp):
        if Tp_curr != 0:
            S = gl.JONSWAP(f_0, Tp_curr, Hs_curr, gamma_mode=gamma_mode)[0]
        else:
            S = float('nan')
        RWI_list.append(np.sqrt(S))
    return RWI_list


def WaveBreak_steep(Hs, Tp, d, steep_crit):
    """"calculates, if the Wave (or seastate) is likley to breaking, eavaluated by wave steepness steep_crit

    Parameters:
        Hs: List, np array or Series, Wave Heigt
        Hs: List, np array or Series, Wave Period
        steep_crit: float, steepness criteria for wave steepness

    Returns:
        break_steep_bool_list: list, list of bools, True if Wave Breaks
        lamda_list: list, list of wavelengths, calculated with dispersion equation
    """

    break_steep_bool_list = []
    lamda_list = []
    steepness_list = []
    for HS, TP in zip(Hs, Tp):

        if TP != 0:
            omega = 1 / TP * 2 * np.pi

            k = gl.k_aus_omega(omega, d)
            lamda = 2 * np.pi / k
            steepness = HS / lamda
            break_steep_bool = steepness > steep_crit

        else:
            lamda = float('nan')
            break_steep_bool = float('nan')
            steepness = float('nan')
        break_steep_bool_list.append(break_steep_bool)
        steepness_list.append(steepness)
        lamda_list.append(lamda)

    return break_steep_bool_list, lamda_list, steepness_list


def angles(mode, N, start, **kwargs):
    width = kwargs.get("width", None)

    ang = []
    ang_mod = []

    if mode == 'full':
        width = 360 / N

    for i, _ in enumerate(range(N)):
        ang.append([start, start + width])
        ang_mod.append([start % 360, (start + width) % 360])
        start = start + width

    return ang, ang_mod


def ExtremeValues(x, N_itter=1000, intervall_mode='percentile', intervall_algorithm='1', T_Return_single=None, time_window_offset=0.5, perc_up=95, perc_down=5, freq_samp=1):
    if T_Return_single is None:
        T_Return_single = [1, 10, 40, 50]

    # sample maxiamal values
    x_max, years_newyear, years_start = gl.x_max_sampeling(x, time_window_offset)

    # sort values by size
    x_sorted_series = x_max.sort_values(ascending=True)
    x_sortet = np.array(x_sorted_series)
    N = len(x_sortet)

    # calculate gumbel coeffs of real data
    beta, mu = gl.Gumbel_coeff(x_max)

    # grid theorie (gumbel)
    x_grid = np.linspace(-(max(x_max) - min(x_max)) * 0.2 + min(x_max),
                         (max(x_max) - min(x_max)) * 0.2 + max(x_max), 100)

    F_theroie = gl.gumbel_P(x_grid, beta, mu)
    f_theroie = gl.gumbel_p(x_grid, beta, mu)

    # real
    F_real = [(n + 0.5) / N for n in range(N)]

    # theoretische erwartungswerte nach gumbelverteilung
    x_theorie = gl.gumbel_P_inv(F_real, beta, mu)

    # confidence intervall über x_theorie
    band_down_xth, _, band_up_xth, _ = gl.Intervall_conf_1(len(x_sortet), N_itter, beta, mu, perc_up, perc_down)

    T_R_x_real = [1 / ((1 - F_curr) * freq_samp) if F_curr != 1 else float("nan") for F_curr in F_real]
    # confidence intervall über TP
    band_down, middle, band_up, T_R_grid = gl.gumbel_conf_intervall(x_max, beta=beta, mu=mu, algorithm=intervall_algorithm, mode=intervall_mode, N_itter=N_itter,
                                                                    T_max=T_Return_single[-1])

    # single Returnperiods
    single_middle = [np.interp(temp_T, T_R_grid, middle) for temp_T in T_Return_single]
    single_down = [np.interp(temp_T, T_R_grid, band_down) for temp_T in T_Return_single]
    single_up = [np.interp(temp_T, T_R_grid, band_up) for temp_T in T_Return_single]

    points = pd.DataFrame({
        "x_max": x_sortet,
        "x_theorie": x_theorie,
        "F_real": F_real,
        "band_down": band_down_xth,
        "band_up": band_up_xth,
        "T_R_x_max": T_R_x_real})

    points.index = x_sorted_series.index

    grid = pd.DataFrame({
        "F_theorie": F_theroie,
        "f_theorie": f_theroie,
        "x_grid": x_grid})

    T_return = pd.DataFrame({"band_down": band_down,
                             "middle": middle,
                             "band_up": band_up,
                             "T_R_grid": T_R_grid})

    T_return_single = pd.DataFrame({
        "T_Return": T_Return_single,
        "middle": single_middle,
        "up": single_up,
        "down": single_down})

    out = {"points": points,
           "grid": grid,
           "years": {"newyear": years_newyear,
                     "start_window": years_start},
           "T_return_single": T_return_single,
           "T_return": T_return,
           "meta": {"N_itter": N_itter,
                    "intervall_mode": intervall_mode,
                    "intervall_algorithm": intervall_algorithm,
                    "perc_up": perc_up,
                    "perc_down": perc_down,
                    "freq_samp": freq_samp}}

    return out


# def write_DEL_database(db_path, DEL):

# alt
def Weibull_fit(Data_Sec, Input):
    global INFO_LOG

    def Weibull_fit(x):
        bin_size, center, count = gl.calculate_histogram(x)

        Weibull, params = gl.fit_weibull_distribution(x, center)

        Weibull = pd.Series(Weibull, index=center)

        return {'curve': Weibull, 'params': params, 'hist_count': count, 'hist_bin_size': bin_size}

    WEIBULL = {}

    INFO_LOG += "Weibull Parameters  \n"

    for key, dict_curr in Data_Sec.items():
        x = dict_curr[Input["col_name_values"]]
        WEIBULL[key] = Weibull_fit(x)

        INFO_LOG += f"{key}: " + "\n"
        INFO_LOG += gl.write_dict(WEIBULL[key]["params"])

    return WEIBULL


def DEl_Condensed(v_m, H_s, T_p, gamma, count, proj_path, exe_path):
    Table = pd.DataFrame(index=v_m)
    Added = pd.DataFrame()

    if gamma is None:
        gamma = 3.3

    isnan = np.isnan(H_s)

    proj_name = os.path.basename(proj_path)

    DEL_raw = gl.calc_JBOOST(exe_path, proj_name, H_s[~isnan], T_p[~isnan], gamma)
    Meta_data = gl.read_lua_values(proj_path, ["design_life", "N_ref", "SN_slope"])

    skal_time = (Meta_data["design_life"] * 365.25 * 24)

    for col in DEL_raw.columns:

        vm_vise = []
        weighted = []

        for DEl_curr, count_curr in zip(DEL_raw[col], count[~isnan]):
            # calculation DEL from count
            DEL_weigted = (DEl_curr ** Meta_data["SN_slope"] * Meta_data["N_ref"]) / skal_time * count_curr

            # including N_ref and SN-slope
            DEL_normalised = (DEL_weigted / Meta_data["N_ref"]) ** (1 / Meta_data["SN_slope"])

            vm_vise.append(DEL_normalised)
            weighted.append(DEL_weigted)

        Table[col] = 0
        Table[col] = Table[col].astype(float)
        Table.loc[isnan.values == False, col] = vm_vise
        Added[col] = pd.Series((np.nansum(weighted) / Meta_data["N_ref"]) ** (
                1 / Meta_data["SN_slope"]))

    Table["count"] = count
    Table.index = v_m

    return Table, Added


def DEL_points(DEL, v_m, v_m_edges, design_life, N_ref, SN_slope):
    Table = pd.DataFrame()
    Added = pd.DataFrame()

    count, _, bin_id = sc.stats.binned_statistic(v_m, v_m, statistic='count', bins=v_m_edges)
    skal_time = (design_life * 365.25 * 24)

    for col in DEL.columns:

        DEL_curr = DEL[col]

        vm_vise = []

        for n_v_bin in range(1, len(v_m_edges)):
            DEL_curr_vm_sec = DEL_curr.loc[bin_id == n_v_bin]

            # calculation DEL from count
            temp = (DEL_curr_vm_sec ** SN_slope * N_ref) / skal_time

            # including N_ref and SN-slope
            DEL_normalised = (np.nansum(temp) / N_ref) ** (1 / SN_slope)

            vm_vise.append(DEL_normalised)

        added = (np.nansum((DEL_curr.values ** SN_slope * N_ref) / skal_time) / N_ref) ** (
                1 / SN_slope)

        Table[col] = vm_vise
        Added[col] = pd.Series(added)

    Table["count"] = count

    v_m_edges = np.array(v_m_edges)
    v_m_mids = (v_m_edges[1:] + v_m_edges[:-1]) / 2

    Table.index = v_m_mids

    return Table, Added


def histogramm(x, x_max=None, x_min=None, auto_size=True, bin_size=None):
    bin_size_soll = bin_size

    if x_min is None:
        x_min = np.min(x)

    if x_max is None:
        x_max = np.max(x)

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

    if auto_size is False:
        if bin_size_soll is None:
            print("   if auto_size is False, bin_size has to be definded")
        else:
            bin_size = bin_size_soll
            N_bins = int(x_max / bin_size)

    edges = np.linspace(x_min - bin_size / 2, x_max + bin_size / 2, N_bins + 2)
    edges = np.round(edges, 4)
    count = [len(x_values[(x_values > edges[i]) & (x_values <= edges[i + 1])]) for i in range(len(edges) - 1)]
    center = (edges[1:] + edges[:-1]) / 2

    return bin_size, center, count


def weibull_fit(x):
    bin_size, center, count = gl.calculate_histogram(x)

    Weibull, params = gl.fit_weibull_distribution(x, center)

    weibull = pd.Series(Weibull, index=center)

    prob = np.array(count) / ((sum(count)) * bin_size)

    return bin_size, center, prob, weibull, params


def cross_correlation(VM_grid, HS_values, HS_grid, TP_values, fill_range=None):
    """"does a cross correlation between the VMHS HS_values(VM_grid) and HSTP TP_values(HS_grid)
    the resulting VMTP TP(VM) condensation the length of HSTP. "fill range" fills the not correlated values by setting TP constant in the specified range

    paramters:
    VM_grid, HS_values: numpy array, same length, can contain nans and dobble values, is expected to be (non strictyl) INCEASING in both vektors
    HS_grid, TP_values, fill_range: numpy array, same length, can contain nans and dobble values, dont has to be increasing!
    fill_range: Series, lenght of HS_grid/TP_values conatining the vm-infomation of the new VMTP as index and a bool for desired filling, no holes!


    return:
    Vm_res, TP_res: numpy arrays containing the new correlation
    """


    is_data_HSTP = ~np.isnan(HS_grid)

    TP_res = np.empty(len(TP_values))
    TP_res[:] = float('nan')

    Vm_res = np.empty(len(TP_values))
    Vm_res[:] = float('nan')

    # Handle NaNs in x vector
    HS_values = gl.fill_nans_constant(HS_values)

    # handle duplicate HS_values in a row
    # Find indices of duplicate x values
    _, unique_indices = np.unique(HS_values, return_index=True)
    duplicates = set(range(len(HS_values))) - set(unique_indices)

    # Remove duplicate x values and corresponding y values for interpolation
    HS_unique = np.delete(HS_values, list(duplicates))
    VM_unique = np.delete(VM_grid, list(duplicates))

    # Perform interpolation on unique values
    Vm_res = np.interp(HS_grid, HS_unique, VM_unique, left=float('nan'), right=float('nan'))

    TP_res[np.where(~np.isnan(Vm_res))[0]] = TP_values[np.where(~np.isnan(Vm_res))[0]]

    # fill not interpolated values to data limits
    if fill_range is not None:
        idx_data_left = np.where(fill_range.values)[0][0]
        idx_data_right = np.where(fill_range.values)[0][-1]

        idx_interp_left = np.where(~np.isnan(Vm_res))[0][0]
        idx_interp_right = np.where(~np.isnan(Vm_res))[0][-1]

        # fill TP_res with constant
        TP_res[idx_data_left:idx_interp_left] = TP_res[idx_interp_left]
        TP_res[idx_interp_right:idx_data_right] = TP_res[idx_interp_right]

        # fill Vm_res with linspaces
        Vm_grid_VMTP = fill_range.index
        Vm_data_left = Vm_grid_VMTP[idx_data_left]
        Vm_data_right = Vm_grid_VMTP[idx_data_right]

        Vm_interp_left = Vm_res[idx_interp_left]
        Vm_interp_right = Vm_res[idx_interp_right]

        Vm_res[idx_data_left:idx_interp_left] = np.linspace(Vm_data_left, Vm_interp_left, idx_interp_left - idx_data_left)
        Vm_res[idx_interp_right:idx_data_right] = np.linspace(Vm_interp_right, Vm_data_right, idx_data_right - idx_interp_right)

    return Vm_res, TP_res


def extreme_contures_blackbox(Hs, Tp, T_return):

    def DVN_steepness(df, h, t, periods, interval):
        import scipy.stats as stats
        ## steepness
        max_y = max(periods)
        X = max_y  # get max 500 year
        period = X * 365.2422 * 24 / interval
        shape, loc, scale = Weibull_method_of_moment(df.hs.values)  # shape, loc, scale
        rve_X = stats.weibull_min.isf(1 / period, shape, loc, scale)

        h1 = []
        t1 = []
        h2 = []
        t2 = []
        h3 = []
        t3 = []
        g = 9.80665
        j15 = 10000
        for j in range(len(t)):
            if t[j] <= 8:
                Sp = 1 / 15
                temp = Sp * g * t[j] ** 2 / (2 * np.pi)
                if temp <= rve_X:
                    h1.append(temp)
                    t1.append(t[j])

                j8 = j  # t=8
                h1_t8 = temp
                t8 = t[j]
            elif t[j] >= 15:
                Sp = 1 / 25
                temp = Sp * g * t[j] ** 2 / (2 * np.pi)
                if temp <= rve_X:
                    h3.append(temp)
                    t3.append(t[j])
                if j < j15:
                    j15 = j  # t=15
                    h3_t15 = temp
                    t15 = t[j]

        xp = [t8, t15]
        fp = [h1_t8, h3_t15]
        t2_ = t[j8 + 1:j15]
        h2_ = np.interp(t2_, xp, fp)
        for i in range(len(h2_)):
            if h2_[i] <= rve_X:
                h2.append(h2_[i])
                t2.append(t2_[i])

        h_steepness = np.asarray(h1 + h2 + h3)
        t_steepness = np.asarray(t1 + t2 + t3)

        return t_steepness, h_steepness
    def Weibull_method_of_moment(X):
        import scipy.stats as stats
        X = X + 0.0001;
        n = len(X);
        m1 = np.mean(X);
        cm1 = np.mean((X - np.mean(X)) ** 1);
        m2 = np.var(X);
        cm2 = np.mean((X - np.mean(X)) ** 2);
        m3 = stats.skew(X);
        cm3 = np.mean((X - np.mean(X)) ** 3);

        from scipy.special import gamma
        def m1fun(a, b, c):
            return a + b * gamma(1 + 1 / c)

        def cm2fun(b, c):
            return b ** 2 * (gamma(1 + 2 / c) - gamma(1 + 1 / c) ** 2)

        def cm3fun(b, c):
            return b ** 3 * (gamma(1 + 3 / c) - 3 * gamma(1 + 1 / c) * gamma(1 + 2 / c) + 2 * gamma(1 + 1 / c) ** 3)

        def cfun(c):
            return abs(np.sqrt(cm3fun(1, c) ** 2 / cm2fun(1, c) ** 3) - np.sqrt(cm3 ** 2 / cm2 ** 3))

        from scipy import optimize
        cHat = optimize.fminbound(cfun, -2, 5)  # shape

        def bfun(b):
            return abs(cm2fun(b, cHat) - cm2)

        bHat = optimize.fminbound(bfun, -5, 30)  # scale

        def afun(a):
            return abs(m1fun(a, bHat, cHat) - m1)

        aHat = optimize.fminbound(afun, -5, 30)  # location

        return cHat, aHat, bHat  # shape, location, scale

    def joint_distribution_Hs_Tp(data, var_hs='hs', var_tp='tp', periods=None, adjustment=None):
        """
        This fuction will plot Hs-Tp joint distribution using LogNoWe model (the Lognormal + Weibull distribution)
        df : dataframe,
        var1 : Hs: significant wave height,
        var2 : Tp: Peak period
        file_out: Hs-Tp joint distribution, optional
        """
        if periods is None:
            periods = [1, 10, 100, 10000]

        if adjustment == 'NORSOK':
            periods_adj = np.array([x * 6 for x in periods])
        else:
            periods_adj = periods

        df = data
        pd.options.mode.chained_assignment = None  # default='warn'
        df.loc[:, 'hs'] = df[var_hs].values
        df.loc[:, 'tp'] = Tp_correction(df[var_tp].values)

        import scipy.stats as stats
        from scipy.optimize import curve_fit
        from scipy.signal import find_peaks

        # calculate lognormal and weibull parameters and plot the PDFs
        mu = np.mean(np.log(df.hs.values))  # mean of ln(Hs)
        std = np.std(np.log(df.hs.values))  # standard deviation of ln(Hs)
        alpha = mu
        sigma = std

        h = np.linspace(start=0.01, stop=30, num=1500)

        if 0 < mu < 5:
            pdf_Hs1 = 1 / (np.sqrt(2 * np.pi) * alpha * h) * np.exp(-(np.log(h) - sigma) ** 2 / (2 * alpha ** 2))
        else:
            param = stats.lognorm.fit(df.hs.values, )  # shape, loc, scale
            pdf_lognorm = stats.lognorm.pdf(h, param[0], loc=param[1], scale=param[2])
            pdf_Hs1 = pdf_lognorm

        param = Weibull_method_of_moment(df.hs.values)  # stats.weibull_min.fit(df.hs.values) # shape, loc, scale
        pdf_Hs2 = stats.weibull_min.pdf(h, param[0], loc=param[1], scale=param[2])

        # Find the index where two PDF cut, between P60 and P99
        for i in range(len(h)):
            if abs(h[i] - np.percentile(df.hs.values, 60)) < 0.1:
                i1 = i

            if abs(h[i] - np.percentile(df.hs.values, 99)) < 0.1:
                i2 = i

        epsilon = abs(pdf_Hs1[i1:i2] - pdf_Hs2[i1:i2])
        param = find_peaks(1 / epsilon)
        try:
            index = param[0][1]
        except:
            try:
                index = param[0][0]
            except:
                index = np.where(epsilon == epsilon.min())[0]
        index = index + i1

        # Merge two functions and do smoothing around the cut
        eta = h[index]
        pdf_Hs = h * 0
        for i in range(len(h)):
            if h[i] < eta:
                pdf_Hs[i] = pdf_Hs1[i]
            else:
                pdf_Hs[i] = pdf_Hs2[i]

        for i in range(len(h)):
            if eta - 0.5 < h[i] < eta + 0.5:
                pdf_Hs[i] = np.mean(pdf_Hs[i - 10:i + 10])

        #####################################################
        # calcualte a1, a2, a3, b1, b2, b3
        # firstly calcualte mean_hs, mean_lnTp, variance_lnTp
        Tp = df.tp.values
        Hs = df.hs.values
        maxHs = max(Hs)
        if maxHs < 2:
            intx = 0.05
        elif 2 <= maxHs < 3:
            intx = 0.1
        elif 3 <= maxHs < 4:
            intx = 0.2
        elif 4 <= maxHs < 10:
            intx = 0.5
        else:
            intx = 1.0

        mean_hs = []
        variance_lnTp = []
        mean_lnTp = []

        hs_bin = np.arange(0, maxHs + intx, intx)
        for i in range(len(hs_bin) - 1):
            idxs = np.where((hs_bin[i] <= Hs) & (Hs < hs_bin[i + 1]))
            if Hs[idxs].shape[0] > 15:
                mean_hs.append(np.mean(Hs[idxs]))
                mean_lnTp.append(np.mean(np.log(Tp[idxs])))
                variance_lnTp.append(np.var(np.log(Tp[idxs])))

        mean_hs = np.asarray(mean_hs)
        mean_lnTp = np.asarray(mean_lnTp)
        variance_lnTp = np.asarray(variance_lnTp)

        # calcualte a1, a2, a3
        parameters, covariance = curve_fit(Gauss3, mean_hs, mean_lnTp)
        a1 = parameters[0]
        a2 = parameters[1]
        a3 = 0.36

        # calcualte b1, b2, b3
        start = 1
        x = mean_hs[start:]
        y = variance_lnTp[start:]
        parameters, covariance = curve_fit(Gauss4, x, y)
        b1 = 0.005
        b2 = parameters[0]
        b3 = parameters[1]

        # calculate pdf Hs, Tp
        t = np.linspace(start=0.01, stop=40, num=2000)

        f_Hs_Tp = np.zeros((len(h), len(t)))
        pdf_Hs_Tp = f_Hs_Tp * 0

        for i in range(len(h)):
            mu = a1 + a2 * h[i] ** a3
            std2 = b1 + b2 * np.exp(-b3 * h[i])
            std = np.sqrt(std2)

            f_Hs_Tp[i, :] = 1 / (np.sqrt(2 * np.pi) * std * t) * np.exp(-(np.log(t) - mu) ** 2 / (2 * std2))
            pdf_Hs_Tp[i, :] = pdf_Hs[i] * f_Hs_Tp[i, :]

        interval = ((df.index[-1] - df.index[0]).days + 1) * 24 / df.shape[0]  # in hours

        t3 = []
        h3 = []
        X = []
        hs_tpl_tph = pd.DataFrame()

        # Assuming Hs_Tp_curve() returns four values, otherwise adjust accordingly
        for i in range(len(periods)):
            t3_val, h3_val, X_val, hs_tpl_tph_val = Hs_Tp_curve(df.hs.values, pdf_Hs, pdf_Hs_Tp, f_Hs_Tp, h, t, interval, X=periods_adj[i])
            t3.append(t3_val)
            h3.append(h3_val)
            X.append(X_val)
            hs_tpl_tph_val.columns = [f'{col}_{periods[i]}' for col in hs_tpl_tph_val.columns]
            hs_tpl_tph = pd.concat([hs_tpl_tph, hs_tpl_tph_val], axis=1)

        # if save_rve:
        #    hs_tpl_tph[3].to_csv(str(param[2])+'_year.csv', index=False)

        return a1, a2, a3, b1, b2, b3, pdf_Hs, h, t3, h3, X, hs_tpl_tph

    def Tp_correction(Tp):
        """
        This function will correct the Tp from ocean model which are vertical straight lines in Hs-Tp distribution
        """
        new_Tp = 1 + np.log(Tp / 3.244) / 0.09525
        index = np.where(Tp >= 3.2)  # indexes of Tp
        r = np.random.uniform(low=-0.5, high=0.5, size=len(Tp[index]))
        Tp[index] = np.round(3.244 * np.exp(0.09525 * (new_Tp[index] - 1 - r)), 1)
        return Tp

    def Hs_Tp_curve(data, pdf_Hs, pdf_Hs_Tp, f_Hs_Tp, h, t, interval, X=100):
        import scipy.stats as stats
        from scipy.signal import find_peaks

        # RVE of X years
        shape, loc, scale = Weibull_method_of_moment(data)  # shape, loc, scale

        if X == 1:
            period = 1.5873 * 365.2422 * 24 / interval
        else:
            period = X * 365.2422 * 24 / interval
        rve_X = stats.weibull_min.isf(1 / period, shape, loc, scale)

        # Find index of Hs=value
        epsilon = abs(h - rve_X)
        param = find_peaks(1 / epsilon)  # to find the index of bottom
        index = param[0][0]  # the  index of Hs=value

        # Find peak of pdf at Hs=RVE of X year
        pdf_Hs_Tp_X = pdf_Hs_Tp[index, :]  # Find pdf at RVE of X year
        param = find_peaks(pdf_Hs_Tp_X)  # find the peak
        index = param[0][0]
        f_Hs_Tp_100 = pdf_Hs_Tp_X[index]

        h1 = []
        t1 = []
        t2 = []
        for i in range(len(h)):
            f3_ = f_Hs_Tp_100 / pdf_Hs[i]
            f3 = f_Hs_Tp[i, :]
            epsilon = abs(f3 - f3_)  # the difference
            para = find_peaks(1 / epsilon)  # to find the bottom
            index = para[0]
            if t[index].shape[0] == 2:
                h1.append(h[i])
                t1.append(t[index][0])
                t2.append(t[index][1])

        h1 = np.asarray(h1)
        t1 = np.asarray(t1)
        t2 = np.asarray(t2)
        t3 = np.concatenate((t1, t2[::-1]))  # to get correct circle order
        h3 = np.concatenate((h1, h1[::-1]))  # to get correct circle order
        t3 = np.concatenate((t3, t1[0:1]))  # connect the last to the first point
        h3 = np.concatenate((h3, h1[0:1]))  # connect the last to the first point

        df = pd.DataFrame()
        df['hs'] = h1
        df['t1'] = t1
        df['t2'] = t2

        return t3, h3, X, df

    def Gauss3(x, a1, a2):
        y = a1 + a2 * x ** 0.36
        return y

    def Gauss4(x, b2, b3):
        y = 0.005 + b2 * np.exp(-x * b3)
        return y

    df = pd.concat((Hs,Tp), axis=1)

    a1, a2, a3, b1, b2, b3, pdf_Hs, h, t3, h3, X, hs_tpl_tph = joint_distribution_Hs_Tp(df, var_hs=Hs.name, var_tp=Tp.name, periods=T_return)

    # calculate pdf Hs, Tp
    t = np.linspace(start=0.01, stop=40, num=2000)

    f_Hs_Tp = np.zeros((len(h), len(t)))
    pdf_Hs_Tp = f_Hs_Tp * 0

    for i in range(len(h)):
        mu = a1 + a2 * h[i] ** a3
        std2 = b1 + b2 * np.exp(-b3 * h[i])
        std = np.sqrt(std2)

        f_Hs_Tp[i, :] = 1 / (np.sqrt(2 * np.pi) * std * t) * np.exp(-(np.log(t) - mu) ** 2 / (2 * std2))
        pdf_Hs_Tp[i, :] = pdf_Hs[i] * f_Hs_Tp[i, :]

    interval = ((df.index[-1] - df.index[0]).days + 1) * 24 / df.shape[0]  # in hours
    t_steepness, h_steepness = DVN_steepness(df, h, t, T_return, interval)

    out = {}
    for i in range(len(X)):
        out[f"{X[i]} years"] = pd.DataFrame()
        out[f"{X[i]} years"]["x"] = h3[i]
        out[f"{X[i]} years"]["y"] = t3[i]

    return out


# %% macro functions
def calc_VMHS(Vm, Hs, angle, angle_grid,
              N_grid=100,
              weight_y=False,
              deg_reg=3,
              model_reg='poly',
              cut_reg=0,
              weighting_reg=0,
              zone_reg=None,
              zone_line=None,
              bin_min=0,
              average_correction=1.0,
              avrg_method='mean',
              make_monotone=False):
    """Returns a list of VMHS segment objects for all segments in angle_grid.
    If angle_grid is None, omnidirectional is returned.

    Arguments:
        Vm: Series, Wind-Speed data; index is stored in Segment.indize Object
            to link used data. Datetime format recommended; Series Name is saved
            in Segment.colnames['x'].
        Hs: Series, Wave-Height data; same index as Vm required; Series Name is
            saved in Segment.colnames['y'].
        angle: Series, Angle data; same index as Vm required; Series Name is
            saved in Segment.angle_name.
        angle_grid: List of List (,2) with angle pairs describing the edges
            of the segments; if None, omnidirectional is returned.

    Optional:
        N_grid: int, default: 100.
        weight_y: bool, default: False.
        deg_reg: int, default: 3.
        model_reg: str, default: 'poly'.
        cut_reg: int, default: 0.
        weighting_reg: int, default: 0.
        zone_reg: list, default: [None, None].
        zone_line: list, default: [None, None].
        bin_min: int, default: 0.
        perc_mean: int, default: 50.
        avrg_method: str, default: 'mean'.
        make_monotone: str, default: 'mean'.

    Returns:
        Data_Out: list of segment objects.
    """

    # Ensure default values for mutable parameters
    zone_reg = zone_reg if zone_reg is not None else [None, None]
    zone_line = zone_line if zone_line is not None else [None, None]


    Data_Out = []
    grid = np.linspace(0, max(Vm), N_grid + 1)

    if angle_grid is None:
        # omni
        df = pd.concat([Vm, Hs], axis=1)
        VMHS_DATA, Coeffs = condensation(df[Vm.name], df[Hs.name], grid,
                                 reg_model=model_reg,
                                 deg_reg=deg_reg,
                                 cut_reg=cut_reg,
                                 reg_weighting=weighting_reg,
                                 zone_reg=zone_reg.copy(),  # Pass the list directly
                                 zone_line=zone_line.copy(),  # Pass the list directly
                                 bin_min=bin_min,
                                 average_correction=average_correction,
                                 avrg_method=avrg_method,
                                 make_monotone=make_monotone)

        temp = Segment(0, angles=None,
                       result={'data': VMHS_DATA, 'coeffs': Coeffs},
                       colnames={'x': Vm.name, 'y': Hs.name},
                       indizes=list(df.index),
                       angle_name=None)

        Data_Out.append(temp)

    else:
        num = 1
        # Grid festlegen

        for angle_segment in angle_grid:
            # Make copies of the mutable parameters to avoid overwriting


            df = pd.concat([Vm, Hs, angle], axis=1)
            df_filt = gl.filter_dataframe(df, angle.name, angle_segment[0], angle_segment[1])

            VMHS_DATA, Coeffs = condensation(df_filt[Vm.name], df_filt[Hs.name], grid,
                                     deg_reg=deg_reg,
                                     reg_model=model_reg,
                                     cut_reg=cut_reg,
                                     reg_weighting=weighting_reg,
                                     zone_reg=zone_reg.copy(),
                                     zone_line=zone_line.copy(),
                                     bin_min=bin_min,
                                     average_correction=average_correction,
                                     avrg_method=avrg_method,
                                     make_monotone=make_monotone)

            temp = Segment(num,
                           angles=[angle_segment[0],
                                   angle_segment[1]],
                           indizes=list(df_filt.index),
                           result={'data': VMHS_DATA, 'coeffs': Coeffs},
                           colnames={'x': Vm.name, 'y': Hs.name},
                           angle_name=angle.name)

            Data_Out.append(temp)
            num = num + 1

    return Data_Out


def calc_HSTP(Hs, Tp, angle, angle_grid, **kwargs):
    """retruns list of HSTP segment objects for all segments in angle_grid, if angle_grid is None, omnidirectional is returned

    Arguments:
        Hs: Series, Wave-Height data, index is stored in Segment.indize Object to link used data, Datetime format recomended, Series Name is saved in Segment.colnames['x']
        Tp: Series, Wave-Period data, same index as Vm required, Series Name is saved in Segment.colnames['y']
        angle: Series, Wave-Height data, same index as Vm required, Series Name is saved in Segment.angle_name
        angle_grid: List of List (,2) with angle pairs discribing the edges of the segments, if None, omnidirectional is returned

    optional:
        N_grid: int, default: 100
        weight_y: bool, default: False
        deg_reg: int, default: 3
        model_reg: str, int, default: 'poly'
        cut_reg: int, default: 0
        weighting_reg: int, default: 0
        zone_reg: list, default: [None,None]
        zone_line: list, default: [None,None]
        bin_min: int, default: 0

    return:
        Data_Out: list of segment objects
    """

    N_grid = kwargs.get('N_grid', 100)
    deg_reg = kwargs.get('deg_reg', 3)
    model_reg = kwargs.get('model_reg', 'poly')
    cut_reg = kwargs.get('cut_reg', 0)
    weighting_reg = kwargs.get('weighting_reg', 0)
    zone_reg = kwargs.get('zone_reg', [None, None])
    zone_line = kwargs.get('zone_line', [None, None])
    bin_min = kwargs.get('bin_min', 0)
    quantile = kwargs.get('quantile', False)
    quant_up = kwargs.get('quant_up', None)
    quant_low = kwargs.get('quant_low', None)
    perc = kwargs.get('percentiles', [33, 66])
    avrg_method = kwargs.get('avrg_method', 'mean')
    make_monotone = kwargs.get('make_monotone', False)

    Data_Out = []

    # Grid festlegen
    grid = np.linspace(0, max(Hs), N_grid + 1)

    if angle_grid is None:
        # omni
        Table_cond, Coeffs = condensation(Hs, Tp, grid,
                                  deg_reg=deg_reg,
                                  reg_model=model_reg,
                                  cut_reg=cut_reg,
                                  reg_weighting=weighting_reg,
                                  zone_reg=zone_reg.copy(),
                                  zone_line=zone_line.copy(),
                                  bin_min=bin_min,
                                  perc=perc,
                                  avrg_method=avrg_method,
                                  make_monotone=make_monotone
                                  )

        if quantile:

            if (len(perc) > 1) & (quant_up is not None) & (quant_low is not None):
                # todo: fehlernachricht
                key_low = f'{perc[0]} percentile plot'
                key_up = f'{perc[-1]} percentile plot'

                Table_cond["quantile"] = quantiles(Table_cond[key_low], Table_cond["mean result plot"], Table_cond[key_up], quant_low, quant_up)

        temp = Segment(0, angles=None,
                       result={'data': Table_cond, 'coeffs': Coeffs},
                       angle_name=None,
                       colnames={'x': Hs.name, 'y': Tp.name},
                       indizes=list(Hs.index))

        Data_Out.append(temp)

    else:
        df = pd.concat([Hs, Tp, angle], axis=1)

        for num, angle_segment in enumerate(angle_grid):

            df_filt = gl.filter_dataframe(df, angle.name, angle_segment[0], angle_segment[1])

            Table_cond, Coeffs = condensation(df_filt[Hs.name], df_filt[Tp.name], grid,
                                      deg_reg=deg_reg,
                                      reg_model=model_reg,
                                      cut_reg=cut_reg,
                                      reg_weighting=weighting_reg,
                                      zone_reg=zone_reg.copy(),
                                      zone_line=zone_line.copy(),
                                      bin_min=bin_min,
                                      perc=perc,
                                      avrg_method=avrg_method,
                                      make_monotone=make_monotone
                                      )

            if quantile:

                if (len(perc) > 1) & (quant_up is not None) & (quant_low is not None):
                    # todo: fehlernachricht
                    key_low = f'{perc[0]} percentile plot'
                    key_up = f'{perc[-1]} percentile plot'

                    Table_cond["quantile"] = quantiles(Table_cond[key_low], Table_cond["mean result plot"], Table_cond[key_up], quant_low, quant_up)

            temp = Segment(num=num,
                           angles=[angle_segment[0], angle_segment[1]],
                           result={'data': Table_cond, 'coeffs': Coeffs},
                           indizes=list(df_filt.index),
                           angle_name=angle.name,
                           colnames={'x': Hs.name, 'y': Tp.name})

            Data_Out.append(temp)

    return Data_Out


def calc_VMTP(vmhs, hstp, vm_points=None, fill_range=False):
    """takes output from VMHS_calc und HSTP_calc (list of angle Sgements) and cross-correlates the results to get the VMTP correlation
    - Takes angle information from vmhs, assumes vmhs and hstp in fitting order of angle-segments!!
    - needs Segment.result object to be a pd.Dataframe with "mean result plot" and "x" in vmhs result dataframes and "x" and "mean result plot" or "quantile" in hstp dataframes ("quantile" overwrites "mean result plot")!)
    - indizes (from basedata) needs to be intialized for counts!
    """

    import matplotlib.pyplot as plt

    VMTP = []
    num = 0
    for vmhs_curr, hstp_curr in zip(vmhs, hstp):
        vmtp_curr_data = pd.DataFrame()

        vmhs_curr_data = vmhs_curr.result["data"]
        hstp_curr_data = hstp_curr.result["data"]

        if "quantile" in hstp_curr_data.keys():
            key_tp = "quantile"
        else:
            key_tp = "mean result plot"

        # VMHS
        HS_values = vmhs_curr_data['mean result plot']
        VM_grid = vmhs_curr_data['x']

        # HSTP
        TP_values = hstp_curr_data[key_tp]
        HS_grid = hstp_curr_data['x']

        # vm_grid aus vm_points
        if vm_points is not None and fill_range:
            vm_points_curr = vm_points[vmhs_curr.indizes]
            grid = np.linspace(0, max(vm_points_curr), len(HS_grid) + 1)
            count, vm_edges, _ = sc.stats.binned_statistic(vm_points_curr, vm_points_curr, statistic='count', bins=grid)
            vm_grid = (vm_edges[:-1] + vm_edges[1:]) / 2
            is_data = count != 0

            first_data = np.where(is_data)[0][0]
            last_data = np.where(is_data)[0][-1]
            is_data[first_data:last_data] = True

            vm_grid = pd.Series(is_data, index=vm_grid)

        else:
            vm_grid = None

        Vm_res, TP_res = cross_correlation(VM_grid.values, HS_values,HS_grid.values, TP_values, fill_range=vm_grid)

        vmtp_curr_data['x'] = Vm_res
        vmtp_curr_data['mean result plot'] = TP_res
        vmtp_curr_data['iscondensation'] = False
        vmtp_curr_data['iscondensation'] = ~np.isnan(Vm_res)

        vmtp_curr = Segment(num, result={'data': vmtp_curr_data}, colnames={'x': vmhs_curr.colnames['x'], 'y': hstp_curr.colnames['y']}, angle_name=vmhs_curr.angle_name, angles=vmhs_curr.angles, indizes=vmhs_curr.indizes)

        VMTP.append(vmtp_curr)
        num = num + 1

    return VMTP


def calc_Roseplot(angle, magnitude, angle_segments):
    r_edges = gl.auto_ticks(0, max(magnitude), fix_end=False, edges=False)
    r_middle = (r_edges[1:] - r_edges[:-1]) / 2
    counts = pd.DataFrame(index=r_middle)
    for angle_segment in angle_segments:
        df = pd.concat([angle, magnitude], axis=1)
        df_filt = gl.filter_dataframe(df, angle.name, angle_segment[0], angle_segment[1])

        count, _, _ = sc.stats.binned_statistic(df_filt[magnitude.name], df_filt[magnitude.name], statistic='count', bins=r_edges)

        angle_midpoint = gl.angle_midpoints([angle_segment[0]], [angle_segment[1]])

        counts[str(angle_midpoint[0])] = pd.Series(count, index=counts.index)

    return counts, r_edges


def calc_RWI(Hs, Tp, angle, angle_grid, f_0):
    """retruns list of RWI segment objects for all segments in angle_grid, if angle_grid is None, omnidirectional is returned

    Arguments:
        Hs: Series, Wave-Height data, index is stored in Segment.indize Object to link used data, Datetime format recomended, Series Name is saved in Segment.colnames['x']
        Tp: Series, Wave-Period data, same index as Vm required, Series Name is saved in Segment.colnames['y']
        angle: Series, Wave-Height data, same index as Vm required, Series Name is saved in Segment.angle_name
        angle_grid: List of List (,2) with angle pairs discribing the edges of the segments, if None, omnidirectional is returned
        f_0: first resonance freq

    return:
        Data_Out: list of segment objects
        RWI:info:  {'RWI_max': RWI_max, 'HS_max': HS_max, 'TP_max': TP_max}
    """
    Data_Out = []
    HS_max = 0
    TP_max = 0
    RWI_max = 0

    if angle_grid is None:

        RWI_list = RWI(Hs, Tp, f_0)

        RWI_df = pd.DataFrame(RWI_list, index=Hs.index)

        temp = Segment(0, angles=None, result=RWI_df, colnames={'x': Hs.name, 'y': Tp.name}, angle_name=None, indizes=list(Hs.index))
        Data_Out.append(temp)

        RWI_max = max(RWI_list)
        HS_max = Hs[RWI_list == RWI_max]
        TP_max = Tp[RWI_list == RWI_max]

    else:

        df = pd.concat([Hs, Tp, angle], axis=1)

        for num, angle_segment in enumerate(angle_grid):

            df_filt = gl.filter_dataframe(df, angle.name, angle_segment[0], angle_segment[1])

            RWI_list = RWI(df_filt[Hs.name], df_filt[Tp.name], f_0)

            RWI_df = pd.DataFrame(RWI_list, index=df_filt[Hs.name].index)

            temp = Segment(num, angles=[angle_segment[0], angle_segment[1]], result=RWI_df, colnames={'x': Hs.name, 'y': Tp.name}, angle_name=angle.name,
                           indizes=list(df_filt[Hs.name].index))
            Data_Out.append(temp)

            if max(RWI_list) > RWI_max:
                RWI_max = max(RWI_list)
                HS_max = df_filt[Hs.name][RWI_list == RWI_max]
                TP_max = df_filt[Tp.name][RWI_list == RWI_max]

    return Data_Out, {'RWI_max': RWI_max, 'HS_max': HS_max, 'TP_max': TP_max}


def calc_tables(vmhs, vm_grid, vm_data):
    VMHS = vmhs
    #    VMTP = hstp

    VMHS_OUT = []
    VMTP_OUT = []

    vm_step = vm_grid[0] - vm_grid[1]
    vm_center = (vm_grid[1:] + vm_grid[:-1]) / 2

    Hsdata = pd.DataFrame()

    Tpdata = pd.DataFrame()
    isDATA_VMHS = pd.DataFrame()
    isDATA_VMTP = pd.DataFrame()
    Count = pd.DataFrame()

    Hsdata['v_m'] = [
        f"{vm_grid[i]} to < {vm_grid[i + 1]}" for i in range(len(vm_grid) - 1)]
    Tpdata['v_m'] = [
        f"{vm_grid[i]} to < {vm_grid[i + 1]}" for i in range(len(vm_grid) - 1)]
    isDATA_VMHS['v_m'] = [
        f"{vm_grid[i]} to < {vm_grid[i + 1]}" for i in range(len(vm_grid) - 1)]
    isDATA_VMTP['v_m'] = [
        f"{vm_grid[i]} to < {vm_grid[i + 1]}" for i in range(len(vm_grid) - 1)]
    Count['v_m'] = [
        f"{vm_grid[i]} to < {vm_grid[i + 1]}" for i in range(len(vm_grid) - 1)]

    # VMHS
    i = 0
    for vmhs_curr in VMHS:
        vmhs_result = vmhs_curr.result["data"]

        # VM and HS data from Regression curves
        HS_VMHS = vmhs_result["mean result plot"]
        VM_VMHS = vmhs_result["x"]

        # get VM grid and count infomation for Table
        index_vmhs_curr = vmhs_curr.indizes
        vm_curr = vm_data[vm_data.index.isin(pd.to_datetime(index_vmhs_curr))]

        count_table, x_edges, _ = sc.stats.binned_statistic(
            vm_curr, vm_curr, statistic='count', bins=vm_grid)

        # 'isdata' for table where count != 0
        Vm_table_center = (x_edges[:-1] + x_edges[1:]) / 2
        isdata = count_table != 0
        count_table = np.array(count_table)

        # set span between fist and last isdata==True to True to catch values in between
        spanned_isdata_table = np.empty(len(Hsdata['v_m']), dtype=bool)
        spanned_isdata_table[:] = False
        spanned_isdata_table[min(
            np.where(isdata)[0]):max(np.where(isdata)[-1] + 1)] = True

        HS = np.empty(len(Hsdata['v_m']))
        HS[:] = np.nan

        # interplate in span for the non nan values of the regression data
        HS[spanned_isdata_table] = sc.interpolate.interp1d(
            VM_VMHS[~np.isnan(HS_VMHS)], HS_VMHS[~np.isnan(HS_VMHS)], fill_value='extrapolate')(
            Vm_table_center[spanned_isdata_table])

        # vm_zone_VMHS = [Vm_table_center[spanned_isdata_table].iloc[0] -
        #                vm_step / 2, Vm_table_center[spanned_isdata_table].iloc[-1] + vm_step / 2]

        result = pd.DataFrame()

        result["vm"] = vm_center
        result["vm_edges"] = [[vm_grid[i], vm_grid[i + 1]] for i in range(len(vm_grid) - 1)]
        result["value"] = HS
        result["isdata"] = isdata
        result["count"] = count_table
        result["iscondensation"] = False
        result["iscondensation"] = spanned_isdata_table

        temp = Segment(i, angles=vmhs_curr.angles, result=result, colnames=vmhs_curr.colnames, indizes=index_vmhs_curr, angle_name=vmhs_curr.angle_name)
        VMHS_OUT.append(temp)
        i = i + 1

    # VMTP
    i = 0

    # for vmtp_curr in VMTP:
    #
    #     vmtp_result = vmtp_curr.result
    #
    #     # VM and TP data from cross correlation, VM has no leading nans and starts late, not gridded,
    #     # TP gridded and starts with leading nans
    #
    #     TP_VMTP = VMTP[ck]["grid"]["mean plot"]
    #     VM_VMTP = VMTP[ck]["grid"]["x"]
    #
    #     # get VM grid and count infomation for Table
    #     df = gl.grid_pointcloud_in_x(
    #         DATA_SEC[ck][COLNAMES["Vm"]], DATA_SEC[ck][COLNAMES["Hs"]], vm_grid)
    #
    #     # 'isdata' for table where count != 0
    #     count_table = df["count"]
    #     Vm_table_center = df["x"]
    #     isdata = count_table != 0
    #
    #     # 'isdata' for table where there is no data because of no cross correlation
    #     bool_corr = (Vm_table_center > min(VM_VMTP[~np.isnan(VM_VMTP)]) - vm_step / 2) & (
    #             Vm_table_center < (max(VM_VMTP[~np.isnan(VM_VMTP)]) + vm_step / 2))
    #     isdata_table = isdata & bool_corr
    #
    #     # set span between fist and last isdata==True to True to catch values in between
    #
    #     spanned_isdata_table = np.empty(len(Hsdata['v_m']), dtype=bool)
    #     spanned_isdata_table[:] = False
    #     spanned_isdata_table[min(
    #         np.where(isdata_table)[0]):max(np.where(isdata_table)[-1] + 1)] = True
    #
    #     TP = np.empty(len(Hsdata['v_m']))
    #     TP[:] = np.nan
    #     TP[spanned_isdata_table] = sc.interpolate.interp1d(
    #         VM_VMTP[~np.isnan(TP_VMTP)], TP_VMTP[~np.isnan(TP_VMTP)], fill_value='extrapolate')(
    #         Vm_table_center.iloc[spanned_isdata_table])
    #
    #     vm_zone_VMTP = [Vm_table_center[spanned_isdata_table].iloc[0] -
    #                     vm_step / 2, Vm_table_center[spanned_isdata_table].iloc[-1] + vm_step / 2]
    #
    #     VMZone_VMTP[ck] = vm_zone_VMTP
    #
    #     Tpdata[ck] = TP
    #     isDATA_VMTP[ck] = isdata_table
    #
    #     i = i + 1
    #
    # OUT = {"VMHS": {}}
    # OUT["VMHS"]["table_content"] = Hsdata
    # OUT["VMHS"]["VM_Zone"] = VMZone_VMHS
    # OUT["VMHS"]["isData"] = isDATA_VMHS
    #
    # OUT["VMTP"] = {}
    # OUT["VMTP"]["table_content"] = Tpdata
    # OUT["VMTP"]["VM_Zone"] = VMZone_VMTP
    # OUT["VMTP"]["isData"] = isDATA_VMTP
    #
    # OUT["count"] = Count
    # OUT["Vm_edges"] = Vm_table_edges

    return VMHS_OUT


def calc_WaveBreak_Steep(Hs, Tp, angle, angle_grid, steep_crit, d):
    """retruns list of WaveBreakSteep segment objects for all segments in angle_grid, if angle_grid is None, omnidirectional is returned

    Arguments:
        Hs: Series, Wave-Height data, index is stored in Segment.indize Object to link used data, Datetime format recomended, Series Name is saved in Segment.colnames['x']
        Tp: Series, Wave-Period data, same index as Vm required, Series Name is saved in Segment.colnames['y']
        angle: Series, Wave-Height data, same index as Vm required, Series Name is saved in Segment.angle_name
        angle_grid: List of List (,2) with angle pairs discribing the edges of the segments, if None, omnidirectional is returned
        steep_crit: speeness criteria
        d: water depth

    return:
        Data_Out: list of segment objects
    """

    Data_Out = []

    if angle_grid is None:

        break_steep_bool_list, lamda_list, steepness = WaveBreak_steep(Hs, Tp, d, steep_crit)
        break_steep_bool = pd.Series(break_steep_bool_list, name='bool_break')
        lamda = pd.Series(lamda_list, name='lamda')
        steepness = pd.Series(steepness, name='steepness')

        result = pd.concat([break_steep_bool, lamda, steepness, pd.Series(data=Hs.values, name=Hs.name)], axis=1, ignore_index=False)

        temp = Segment(0, angles=None, result=result, colnames={'x': Hs.name, 'y': Tp.name}, angle_name=None, indizes=list(Hs.index))
        Data_Out.append(temp)

    else:
        df = pd.concat([Hs, Tp, angle], axis=1)

        for num, angle_segment in enumerate(angle_grid):
            df_filt = gl.filter_dataframe(df, angle.name, angle_segment[0], angle_segment[1])

            break_steep_bool_list, lamda_list, steepness = WaveBreak_steep(df_filt[Hs.name], df_filt[Tp.name], d, steep_crit)
            break_steep_bool = pd.Series(break_steep_bool_list, name='bool_break')
            lamda = pd.Series(lamda_list, name='lamda')
            steepness = pd.Series(steepness, name='steepness')

            result = pd.concat([break_steep_bool, lamda, steepness], axis=1)

            temp = Segment(num, angles=[angle_segment[0], angle_segment[1]], result=result, colnames={'x': Hs.name, 'y': Tp.name}, angle_name=angle.name,
                           indizes=list(df_filt.index))
            Data_Out.append(temp)

    return Data_Out


def calc_angle_deviation_tables(angle_orig, angle_comp, v_m, angle_grid, **kwargs):
    """creates missalinment tables for every angle section of angle_orig stored in angle_grid. 
    The rows are definded by a histogram of the v_m vektor. 
    The columns are named after the missalinment of the compare angle angle_comp and is calcultate with: diff =  angle_com - angle_orig"""

    v_m_zone = kwargs.get('v_m_zone', None)
    v_m_step = kwargs.get('v_m_step', 1)
    N_angle_comp_sec = kwargs.get('N_angle_comp_sec', 12)
    Data_Out = []
    indexes = angle_orig.index

    # v_m zone generation
    if v_m_zone is None:
        v_m_zone = [0, max(v_m.values)]

    if v_m_zone[1] is None:
        v_m_zone[1] = max(v_m.values)

    v_m_range = gl.range_stepfix(v_m_step, v_m_zone)

    comp_ang, comp_ang_mod = angles('full', N_angle_comp_sec, -15)

    df = pd.DataFrame({"angle_orig": angle_orig.values, "angle_comp": angle_comp.values, "v_m": v_m}, index=indexes)

    for num, angle_segment in enumerate(angle_grid):

        df_filt = gl.filter_dataframe(df, "angle_orig", angle_segment[0], angle_segment[1])
        midpoint_orig = gl.angle_midpoints([angle_segment[0]], [angle_segment[1]])[0]
        AngleDev = pd.DataFrame()
        AngleDev.index = (v_m_range[1:] + v_m_range[:-1]) / 2
        missalignments = []
        midpoints_comp = []

        for comp_ang_mod_curr in comp_ang_mod:

            midpoint_comp = gl.angle_midpoints([comp_ang_mod_curr[0]], [comp_ang_mod_curr[1]])[0]

            missalignment = (midpoint_comp - midpoint_orig) % 360

            missalignments.append(missalignment)
            midpoints_comp.append(midpoint_comp)

            df_filt_filt = gl.filter_dataframe(df_filt, "angle_comp", comp_ang_mod_curr[0], comp_ang_mod_curr[1])

            if len(df_filt_filt) != 0:
                Count, _, _ = sc.stats.binned_statistic(df_filt_filt["v_m"].values, df_filt_filt["v_m"].values, statistic='count', bins=v_m_range)
            else:
                Count = np.zeros(len(v_m_range) - 1)

            Count = Count / len(df_filt)

            AngleDev[midpoint_comp] = Count

        # sorting

        _, count_cols = zip(*sorted(zip(missalignments, midpoints_comp)))
        count_cols = count_cols
        AngleDev = AngleDev.reindex(columns=count_cols)

        missalinment_rename = {midpoints_comp[i]: missalignments[i] for i in range(len(missalignments))}

        AngleDev = AngleDev.rename(columns=missalinment_rename)

        temp = Segment(num,
                       angles=[angle_segment[0], angle_segment[1]],
                       result=AngleDev,
                       colnames={'ang_orig': angle_orig.name, 'ang_comp': angle_comp.name, 'v_m': v_m.name},
                       angle_name=angle_orig.name,
                       indizes=list(df_filt.index))

        Data_Out.append(temp)

    return Data_Out


def calc_ExtemeValues(x, angles, angle_grid, T_return_single=None, conf_inter_mode=None, conf_inter_algorithm=None, N_itter=None, freq_samp=None, perc_up=None, perc_down=None,
                      time_window_offset=None):
    Data_Out = []

    if angle_grid is None:
        result = ExtremeValues(x,
                               intervall_mode=conf_inter_mode,
                               intervall_algorithm=conf_inter_algorithm,
                               T_Return_single=T_return_single,
                               time_window_offset=time_window_offset,
                               perc_up=perc_up,
                               perc_down=perc_down,
                               freq_samp=freq_samp,
                               N_itter=N_itter)

        temp = Segment(0, angles=None, result=result, colnames={'x': x.name, 'angle':  None}, angle_name=None, indizes=list(x.index))
        Data_Out.append(temp)

    else:
        df = pd.concat([x, angles], axis=1)
        for num, angle_segment in enumerate(angle_grid):
            df_filt = gl.filter_dataframe(df, angles.name, angle_segment[0], angle_segment[1])
            result = ExtremeValues(df_filt[x.name],
                                   intervall_mode=conf_inter_mode,
                                   intervall_algorithm=conf_inter_algorithm,
                                   T_Return_single=T_return_single,
                                   time_window_offset=time_window_offset,
                                   perc_up=perc_up,
                                   perc_down=perc_down,
                                   freq_samp=freq_samp,
                                   N_itter=N_itter)

            temp = Segment(0, angles=angle_segment, result=result, colnames={'x': x.name, 'angle': angles.name}, angle_name=angles.name, indizes=list(df_filt.index))
            Data_Out.append(temp)
    return Data_Out


def update_DEL_db(db_path, Hs, Tp, gamma, proj_path=None, input_path=None, exe_path=None):
    """function, that handels updating DEL values from a Database or calculates the missing information using JBOOST

        tries to calculate missing information using JBOOST, but also compares the Metadata provided in proj_path and input_path for no unnessesary calculations
        """

    if os.path.normpath(proj_path) == os.path.normpath(exe_path):
        print('   invalid JBOOST procject file location, do no place inside ./JBOOST')
        return None
    else:
        JBOOST_proj_path_new = shutil.copy(proj_path, exe_path)
    proj_name = os.path.basename(proj_path)

    if os.path.normpath(input_path) == os.path.normpath(exe_path):
        print('   invalid JBOOST procject input file location, do no place inside ./JBOOST')
        return None
    else:
        shutil.copy(input_path, exe_path)

    Var_lua = gl.read_lua_values(JBOOST_proj_path_new, ['seabed_level', 'design_life', 'N_ref', 'SN_slope', 'res_Nodes'])

    Meta_Curr = {"d": -Var_lua["seabed_level"],
                 "design_life": Var_lua["design_life"],
                 "N_ref": Var_lua["N_ref"],
                 "SN_slope": Var_lua["SN_slope"],
                 "Hs": Hs.name,
                 "Tp": Tp.name,
                 "gamma": gamma.name,
                 'input': input_path}

    print(f'   used variables from {proj_path}:')
    print(f"   {gl.write_dict(Var_lua)}")

    table_name = gl.check_meta_in_valid_db(db_path, Meta_Curr)

    if len(table_name) != 0:
        table_name = table_name[0]
        print('   data with same input parameters found, checking for nodes and timeframe')
        DEL_base_sql = gl.export_df_from_sql(db_path, table_name)

        node_coles = DEL_base_sql.columns
        columnnames_data = list(node_coles)

        # node difference
        nodes = []
        for col in node_coles:
            node = ''.join([s for s in col if s.isdigit()])
            node = float(node)
            nodes.append(node)

        nodes_db = np.unique(nodes)

        nodes_proj = Var_lua["res_Nodes"][1:-1].split(',')
        nodes_proj = [float(node) for node in nodes_proj]

        nodes_proj = set(nodes_proj)
        nodes_db = set(nodes_db)

        new_nodes = nodes_proj.difference(nodes_db)
        both_nodes = nodes_proj.intersection(nodes_db)
        if len(new_nodes) > 0:
            nodes_to_add = True
        else:
            nodes_to_add = False

        # find data, which covers all nans of the nodes provided in nodes_proj
        SQL_temp = gl.filter_df_cols_by_keywords(DEL_base_sql, [str(int(node)) for node in nodes_proj])
        SQL_temp = SQL_temp.dropna()

        # timeframe difference
        in_df = Hs.index

        in_db = SQL_temp.index
        df2_in_df1 = Hs.index.isin(in_db)
        in_df_not_in_db = Hs.index[~df2_in_df1]

        if len(in_df_not_in_db) > 0:
            time_to_add = True
        else:
            time_to_add = False

        if time_to_add and not nodes_to_add:
            print(f'   New time found in run (N = {len(in_df_not_in_db)}), but no new nodes. Calculating missing timepoints to add to db {table_name}, this might take a long time')
            # write database data, which will not be calculated in db frame

            # wirte calculated data, which was not in database in db frame
            DEL_temp = gl.calc_JBOOST(exe_path, proj_name, Hs.loc[in_df_not_in_db], Tp.loc[in_df_not_in_db], gamma.loc[in_df_not_in_db])

        if not time_to_add and nodes_to_add:
            print(f'   New nodes found in run ({new_nodes}), but no new timestep. Calculating the new nodes, this might take a long time')

            gl.write_lua_variables(JBOOST_proj_path_new, {'res_Nodes': new_nodes})

            # write calculated dataframe where it is calculated in db frame
            DEL_temp = gl.calc_JBOOST(exe_path, proj_name, Hs, Tp, gamma)

        if time_to_add and nodes_to_add:
            print(f'   New nodes found in run ({new_nodes}).')
            print(f'   New time found in run (N = {len(in_df_not_in_db)}). Calculating nodes and timepoints, this might take a long time')

            gl.write_lua_variables(JBOOST_proj_path_new, {'res_Nodes': new_nodes})
            # calculate all new nodes and timepoints in dataframe and write in db_frame
            DEL_temp_01 = gl.calc_JBOOST(exe_path, proj_name, Hs.loc[in_df], Tp.loc[in_df], gamma.loc[in_df])

            # calculate the new data for the other nodes
            gl.write_lua_variables(JBOOST_proj_path_new, {'res_Nodes': both_nodes})
            DEL_temp_02 = gl.calc_JBOOST(exe_path, proj_name, Hs.loc[in_df_not_in_db], Tp.loc[in_df_not_in_db], gamma.loc[in_df])

            DEL_temp = gl.merge_dataframes(DEL_temp_02, DEL_temp_01)

        if time_to_add or nodes_to_add:

            #  wait_statement = input("overwrite data in database with current calculation? \n write y for yes and n for no: ")
            wait_statement = 'y'
            if wait_statement == 'y':
                DEL_save = gl.merge_dataframes(DEL_temp, DEL_base_sql)

                write_db = True
            else:
                write_db = False
                DEL_save = None

        else:
            print('   database with same input parameters, same nodes and same or wider timeframe found, no calculation needed')
            write_db = False

    else:
        print(f"   crate new database table and calculate all {len(Hs)} datapoints for all nodes, this might take a long time")
        DEL_save = gl.calc_JBOOST('./JBOOST/', proj_name, Hs, Tp, gamma)

        write_db = True
        columnnames_data = list(DEL_save.columns)

    if write_db:
        table_name = gl.write_DEL_base(db_path, DEL_save, Meta_Curr)

    return table_name, columnnames_data


def calc_Validation(df_DEL, Vm, angle, vmhs_table, vmtp_table, proj_path, input_path, exe_path):
    """retruns list of Validation segment objects for all segments in angle_grid, if angle_grid is None, omnidirectional is returned

    Arguments:
        Vm: Series, Wind-Speed data, index is stored in Segment.indize Object to link used data, Datetime format recomended, Series Name is saved in Segment.colnames['x']
        angle: Series, Wave-Height data, same index as Vm required, Series Name is saved in Segment.angle_name

    optional:

    return:
        Data_Out: list of segment objects
    """
    # extract Nodes from DEL data
    node_coles = df_DEL.columns
    Data_Out = []
    # node difference
    nodes = []
    for col in node_coles:
        node = ''.join([s for s in col if s.isdigit()])
        node = float(node)
        nodes.append(node)

    nodes_db = set(nodes)

    if os.path.normpath(proj_path) == os.path.normpath(exe_path):
        print('   invalid JBOOST procject file location, do no place inside ./JBOOST')
        return None
    else:
        JBOOST_proj_path_new = shutil.copy(proj_path, exe_path)

    if os.path.normpath(input_path) == os.path.normpath(exe_path):
        print('   invalid JBOOST procject input file location, do no place inside ./JBOOST')
        return None
    else:
        shutil.copy(input_path, exe_path)

    gl.write_lua_variables(JBOOST_proj_path_new, {'res_Nodes': nodes_db})

    Meta_data = gl.read_lua_values(JBOOST_proj_path_new, ["design_life", "N_ref", "SN_slope"])

    num = 0
    for vmhs_curr, vmtp_curr in zip(vmhs_table, vmtp_table):

        table_condensed, added_condensed = DEl_Condensed(vmtp_curr.result["vm"],
                                                         vmhs_curr.result["value"],
                                                         vmtp_curr.result["value"],
                                                         None,
                                                         vmhs_curr.result["count"],
                                                         JBOOST_proj_path_new,
                                                         exe_path)

        angle_segment = vmhs_curr.angles

        if angle_segment is not None:
            id_filt = gl.filter_dataframe(angle, [angle.name], [angle_segment[0]], [angle_segment[1]]).index
            angle_name = angle.name
        else:
            id_filt = angle.index
            angle_name = None

        v_m_edges = [vm_curr[0] for vm_curr in vmhs_curr.result["vm_edges"].values]
        v_m_edges.append(vmhs_curr.result["vm_edges"].values[-1][1])

        table_points, added_points = DEL_points(df_DEL.loc[id_filt], Vm.loc[id_filt], v_m_edges, Meta_data["design_life"], Meta_data["N_ref"], Meta_data["SN_slope"])

        out = {"condensed": {"vm_vise": table_condensed,
                             "added": added_condensed},
               "hindcast": {"vm_vise": table_points,
                            "added": added_points},
               "meta": Meta_data}

        V_m_name = vmhs_curr.colnames["x"]
        H_s_name = vmhs_curr.colnames["y"]
        T_p_name = vmtp_curr.colnames["y"]

        temp = Segment(num, angles=angle_segment, indizes=list(id_filt), result=out, angle_name=angle_name,
                       colnames={'Hindcast': {'v_m': V_m_name, 'H_s': H_s_name, 'T_p': T_p_name},
                                 'DELs': list(table_points.columns)})

        Data_Out.append(temp)
        num = num + 1

    return Data_Out


def calc_histogram(x, angle, angle_grid):
    Data_Out = []

    if angle_grid is None:
        # omni
        bin_size, center, count = histogramm(x)
        histo_data = {"bin_size": bin_size, "center": center, "count": count}

        temp = Segment(0, angles=None,
                       result=histo_data,
                       angle_name=None,
                       colnames={'x': x.name},
                       indizes=list(x.index))

        Data_Out.append(temp)

    else:
        num = 1
        # Grid festlegen

        for angle_segment in angle_grid:
            df = pd.concat([x, angle], axis=1)
            df_filt = gl.filter_dataframe(df, angle.name, angle_segment[0], angle_segment[1])

            bin_size, center, count = histogramm(df_filt[x.name])

            histo_data = {"bin_size": bin_size, "center": center, "count": count}

            temp = Segment(num, angles=[angle_segment[0], angle_segment[1]], indizes=list(df_filt.index), result=histo_data, angle_name=angle.name,
                           colnames={'x': x.name})
            Data_Out.append(temp)
            num = num + 1

    return Data_Out


def calc_weibull(x, angle, angle_grid):
    Data_Out = []

    if angle_grid is None:
        # omni
        bin_size, center, prob, weibull, weibull_params = weibull_fit(x)
        histo_data = {"bin_size": bin_size,
                      "center": center,
                      "prob": prob,
                      "weibull": weibull,
                      "weibull_params": weibull_params
                      }

        temp = Segment(0, angles=None,
                       result=histo_data,
                       angle_name=None,
                       colnames={'x': x.name},
                       indizes=list(x.index))

        Data_Out.append(temp)

    else:
        num = 1
        # Grid festlegen

        for angle_segment in angle_grid:
            df = pd.concat([x, angle], axis=1)
            df_filt = gl.filter_dataframe(df, angle.name, angle_segment[0], angle_segment[1])

            bin_size, center, prob, weibull, weibull_params = weibull_fit(df_filt[x.name])
            histo_data = {"bin_size": bin_size,
                          "center": center,
                          "prob": prob,
                          "weibull": weibull,
                          "weibull_params": weibull_params
                          }

            temp = Segment(num, angles=[angle_segment[0], angle_segment[1]], indizes=list(df_filt.index), result=histo_data, angle_name=angle.name,
                           colnames={'x': x.name})
            Data_Out.append(temp)
            num = num + 1

    return Data_Out


def calc_extreme_contures(Hs, Tp, angle, angle_grid, T_return):
    Data_Out = []

    if angle_grid is None:
        # omni
        out = extreme_contures_blackbox(Hs, Tp, T_return)

        temp = Segment(0, angles=None,
                       result=out,
                       angle_name=None,
                       colnames={'x': Hs.name, 'y': Tp.name},
                       indizes=list(Hs.index))

        Data_Out.append(temp)

    else:
        num = 1
        # Grid festlegen

        for angle_segment in angle_grid:

            df = pd.concat([Hs, Tp, angle], axis=1)
            df_filt = gl.filter_dataframe(df, angle.name, angle_segment[0], angle_segment[1])

            out = extreme_contures_blackbox(df_filt[Hs.name], df_filt[Tp.name], T_return)

            temp = Segment(num, angles=[angle_segment[0], angle_segment[1]], indizes=list(df_filt.index), result=out, angle_name=angle.name,
                           colnames={'x': Hs.name, 'y': Tp.name})
            Data_Out.append(temp)
            num = num + 1

    return Data_Out

# %% DataBaseHandling


# load data f
