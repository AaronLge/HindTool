import argparse
import inspect
import os
import sqlite3
from datetime import datetime
import chardet
import pandas as pd
from libaries import general as gl

# %% startup

pd.options.mode.chained_assignment = None  # default='warn'

script_name = os.path.basename(__file__)

parser = argparse.ArgumentParser(description=script_name)
parser.add_argument('-i', metavar='path_in', required=False, type=str, default=None,
                    help='the filepath to the input file, if empty "Input.txt"')
parser.add_argument('-o', metavar='path_out', required=False, type=str,
                    help='the filepath to the output dir, if empty, taken from Input file')

args = parser.parse_args()

path_in = args.i

if path_in is None:
    path_in = 'C:/Users/aaron.lange/Desktop/Projekte/Hindcast_Tool/HindTool/db_example/db_meta.txt'
    path_in = 'C:/Users/aaron.lange/Desktop/Projekte/Hindcast_Tool/Data_Omexon/db_meta.txt'

filename = inspect.getframeinfo(inspect.currentframe()).filename
path_main = os.path.dirname(os.path.abspath(filename))

# %% reading input file
print(f"reading Inputfile ({path_in})...")

INPUT = gl.read_input_txt(path_in)

timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")


# %% functions
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


def csv_to_df(path_csvs, resample_rate, data_kind='MetOcean', encoding='auto', nans=None, skiprows=None, delimiter=';', dayfirst=False, datetime_mode='single_col',
              low_memory=True, drop_rows=None, drop_cols=None):
    """accepts Dict of csv paths (Metocean Format!) and stores the individual data as well a resampled version in
    a combined table to a sql database. If the database already exists, it will be overwritten

    Parameters:
        Paths: dict with the keys: {path_wave, path_atmo, path_ocean} and the corresponding paths to the csv files, can be None
        db_name: sting to name the database, with path if necessary
        resample_rate: resample rate in the fomrat: float {y,d,s,m} for year, day, second, month. Example: "1d"

    Return:
        sql Database at the desired path
    """
    data_resampled = []
    DF_Raw = {}
    csv_files = [os.path.join(path_csvs, f) for f in os.listdir(path_csvs) if f.endswith(('.csv', '.CSV'))]

    for file in csv_files:

        file_name = os.path.basename(file)
        print(f"reading  {file_name}")

        if encoding == 'auto':
            encoding_full = check_encoding(file)

            if encoding_full["confidence"] == 1:
                encoding = encoding_full["encoding"]

            else:
                print(f"encoding confidence less than 1 ({encoding_full['confidence']}), 'utf-8' used.")
                encoding = 'utf-8'

        if data_kind == 'MetOcean':
            nans = None
            skiprows = 15
            delimiter = ','
            dayfirst = False
            datetime_mode = 'single_col'
            drop_rows = None
            drop_cols = None

        if data_kind == 'APGMer_floating':
            nans = ['', ' ', '-']
            skiprows = 6
            delimiter = ','
            dayfirst = False
            datetime_mode = 'single_col'
            low_memory = False
            drop_rows = 0
            drop_cols = 0

        if data_kind == 'APGMer_tight':
            nans = ['', ' ', 'NaN']
            skiprows = None
            delimiter = ','
            dayfirst = False
            datetime_mode = 'multi_col'
            low_memory = False
            drop_rows = 0
            drop_cols = None

        df_NAN = pd.read_csv(file, skiprows=skiprows, encoding=encoding, delimiter=delimiter, na_values=nans, low_memory=low_memory, na_filter=True)

        if drop_rows is not None:
            df_NAN.drop(index=drop_rows, axis=0, inplace=True)

        if drop_cols is not None:
            df_NAN.drop(df_NAN.columns[drop_cols], axis=1, inplace=True)

        if datetime_mode == 'single_col':
            df_datetime = df_NAN.loc[:, df_NAN.columns[0]]
            df_NAN = df_NAN.drop(df_NAN.index[drop_rows])
            df_datetime = df_datetime.drop(df_datetime.index[drop_rows])
            df_NAN.drop(df_NAN.columns[[0]], axis=1, inplace=True)
            df_datetime = pd.to_datetime(df_datetime, dayfirst=dayfirst)

        elif datetime_mode == 'multi_col':
            df_datetime = df_NAN[df_NAN.columns[[0, 1, 2, 3, 4]]]

            df_NAN.drop(df_NAN.columns[[0, 1, 2, 3, 4]], axis=1, inplace=True)
            df_datetime = pd.to_datetime(df_datetime, dayfirst=dayfirst)

        else:
            print("please choose from 'single_col' or 'multi_col'")
            return

        df_NAN.index = df_datetime
        df = df_NAN.dropna(how='all')

        df = df.astype(float)
        df.index.name = 'index'
        DF_Raw[file_name] = df

        df_resample = df.resample(resample_rate).mean()

        data_resampled.append(df_resample)

    df_ges = join_dataframes_by_index(data_resampled)
    df_ges.index.name = 'index'
    DF_comb = df_ges

    return {'raw': DF_Raw, 'combined': DF_comb}

# %%Main
if INPUT["General"]["xlsx_concvert"]:
    gl.xlsx2csv(INPUT["General"]["xlsx_path"], INPUT["General"]["csv_path"], exclude_sheets=INPUT["General"]["exclude_sheets"])

DFs = csv_to_df(INPUT["General"]["csv_path"],
                INPUT["General"]["resample_rate"],
                data_kind=INPUT["General"]["type"],
                encoding=INPUT["General"]["encoding"],
                nans=INPUT["General"]["nans"],
                skiprows=INPUT["General"]["skiprows"],
                delimiter=INPUT["General"]["delimiter"],
                dayfirst=INPUT["General"]["dayfirst"],
                datetime_mode=INPUT["General"]["datetime_mode"],
                low_memory=INPUT["General"]["low_memory"],
                drop_rows=INPUT["General"]["drop_rows"],
                drop_cols=INPUT["General"]["drop_cols"])

if INPUT["General"]["type"] is None:
    INPUT["General"]["type"] = 'manual'

META = pd.DataFrame(columns=["type",
                             "Source",
                             "Date created",
                             "Point Name",
                             "Longitude",
                             "Latitude",
                             "Projection",
                             "Water Depth",
                             "Start Date",
                             "End Date",
                             "Time Step",
                             "Number of samples",
                             "Cell size"],
                    index=list(DFs["raw"].keys()) + ["Combined"])

for name, df in DFs["raw"].items():

    for paramerter in META.columns:
        if INPUT[name].get(paramerter, {}):
            META.loc[name, paramerter] = INPUT[name][paramerter]

    META.loc[name, "type"] = INPUT["General"]["type"]

    if INPUT[name]["Time Step"] == 'auto':
        META.loc[name, "Time Step"] = (df.index[1] - df.index[0]).total_seconds()

    if INPUT[name]["Number of samples"] == 'auto':
        META.loc[name, "Number of samples"] = len(df)

    if INPUT[name]["Start Date"] == 'auto':
        META.loc[name, "Start Date"] = str(DFs["combined"].index[0])

    if INPUT[name]["End Date"] == 'auto':
        META.loc[name, "End Date"] = str(DFs["combined"].index[-1])

META.loc["Combined", "Date created"] = str(datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
META.loc["Combined", "Start Date"] = str(DFs["combined"].index[0])
META.loc["Combined", "End Date"] = str(DFs["combined"].index[-1])
META.loc["Combined", "Time Step"] = (DFs["combined"].index[1] - DFs["combined"].index[0]).total_seconds()
META.loc["Combined", "Number of samples"] = len(DFs["combined"])

db_path = INPUT["General"]["db_path"]
db_exists = os.path.exists(db_path)

if db_exists:
    overwrite = input(f"Database {db_path} already exists, overwrite? (y/n)")
    if overwrite == 'y':
        os.remove(db_path)

conn = sqlite3.connect(db_path)

META.to_sql('Hind_MetaData', conn)

for name, df in DFs["raw"].items():
    print(f"adding {name} to database")

    df.to_sql("Hind_raw_" + name, conn)

print(f"adding Combined data to database")

DFs["combined"].to_sql("Hind_combined", conn)

conn.close()

print(f"finished, saved at {INPUT['General']['db_path']}")
