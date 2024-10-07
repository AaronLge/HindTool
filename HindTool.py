import argparse
import datetime
import inspect
import os
import re
import shutil
import subprocess
from warnings import simplefilter
import numpy as np
import pandas as pd
import scipy as sc
import warnings
import re
import sys
import matplotlib
# matplotlib.use('TkAgg',force=True)
from matplotlib.colors import LinearSegmentedColormap

path = r"C:\\temp\\python_self_crated\\packages"
sys.path.insert(0, path)

from allib import general as gl
from allib import hindtoolcalc as hc_calc
from allib import hindtoolplot as hc_plt


# %% FUNCTIONS - General

def modify_lua_variables(file_path, variables):
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


def Series_to_txt(series, path):
    with open(path, "w") as text_file:
        for time, row in series.items():
            string = f"{time}" + '\t' + f"{row}" + '\n'
            text_file.write(string)
    return


def Data_out_csv(Data_Out, Name, path_csv):
    global path_out

    i = 0

    for dict_dict_name, cd_dict in Data_Out.items():
        dict_dict_name = dict_dict_name.replace(':', '=')
        dict_dict_name = dict_dict_name.replace(' ', '')

        i = i + 1
        i_print = str(i).zfill(2)

        path_csv_full = path_csv + Name + '_' + i_print + '_' + dict_dict_name + '.csv'

        cd_dict.to_csv(path_csv_full)


def Data_out_xls(Data_Out, Name, path_csv):
    global path_out

    i = 0

    with pd.ExcelWriter(path_csv + Name) as writer:
        for dict_dict_name, cd_dict in Data_Out.items():
            dict_dict_name = dict_dict_name.replace(':', '=')
            dict_dict_name = dict_dict_name.replace(' ', '')

            cd_dict.to_excel(writer, sheet_name=dict_dict_name, index=False)

            i = i + 1


# %% Startup

simplefilter(action="ignore", category=pd.errors.PerformanceWarning)

warnings.filterwarnings("ignore", category=RuntimeWarning)

pd.options.mode.chained_assignment = None  # default='warn'

script_name = os.path.basename(__file__)

parser = argparse.ArgumentParser(description=script_name)
parser.add_argument('-i', metavar='path_in', required=False, type=str, default='Input.txt',
                    help='the filepath to the input file, if empty "Input.txt"')
parser.add_argument('-o', metavar='path_out', required=False, type=str,
                    help='the filepath to the output dir, if empty, taken from Input file')

args = parser.parse_args()

path_in = args.i

filename = inspect.getframeinfo(inspect.currentframe()).filename
path_main = os.path.dirname(os.path.abspath(filename))

DATA_OUT = {}

Colors = {
    'JBO_grey': '#8e8778',
    'JBO_green': '#008f85',
    'Black': '#000000',
    'Pink': '#FF00FF',
    'White': '#FFFFFF',
    'Aqua': '#00FFFF',
}

INFO_LOG = str()

print("\n***Starting " + f"{script_name}" +
      ", this might take a few minutes***\n")

# %% UserInput
print(f"reading Inputfile ({path_in})...")

INPUT = gl.read_input_txt(path_in)
timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

# #ReportInput
# if INPUT.get("Report", {}):
#     try:
#         INPUT_REPORT = gl.read_input_txt(INPUT["Report"]["ReportInput"])
#
#     except:
#         print(f"Report input file {INPUT['Report']['ReportInput']} not fond")

if args.o is None:
    if INPUT['DataOut']['dir_name'] is None:
        path_out = INPUT['DataOut']['path_out'] + 'HindCast_' + timestamp + '/'
    else:
        path_out = INPUT['SELECT']['path_out'] + \
                   INPUT['SELECT']['dir_name'] + '/'
else:
    path_out = args.o + '/'

INFO_LOG += f"Path_out = {path_out}" + "\n \n"

print(f"Path_out = {path_out}")

if not os.path.exists(path_out):
    os.makedirs(path_out)

shutil.copy(path_in, path_out + 'Input.txt')

# load from lua?
if INPUT['Structure']['d_from_proj']:
    JBOOST_proj_Path = INPUT['DataBase']['JBOOST_proj']
    Var_lua = gl.read_lua_values(JBOOST_proj_Path, ['seabed_level'])

    INPUT["Structure"]["d"] = -Var_lua["seabed_level"]

    print(f'loaded d = {INPUT["Structure"]["d"]} from {JBOOST_proj_Path}')

if INPUT['Filter']['timeframe']:
    timeframe = (pd.to_datetime(INPUT['Filter']['datetime_start']), pd.to_datetime(INPUT['Filter']['datetime_end']))
else:
    timeframe = None
# %% DataRead

if INPUT['DataBase']["CreateDataBase"]:
    db_path = hc_calc.csv_to_sqlDB(INPUT["DataBase"]["path_csvs"],
                                   INPUT["DataBase"]["name_DataBase"],
                                   INPUT["DataBase"]["resample_rate"],
                                   data_kind=INPUT["DataBase"]["mode"],
                                   encoding=INPUT["DataBase"]["encoding"],
                                   nans=INPUT["DataBase"]["nans"],
                                   skiprows=INPUT["DataBase"]["skiprows"],
                                   delimiter=INPUT["DataBase"]["delimiter"],
                                   dayfirst=INPUT["DataBase"]["dayfirst"],
                                   datetime_mode=INPUT["DataBase"]["datetime_mode"],
                                   low_memory=INPUT["DataBase"]["low_memory"],
                                   drop_rows=INPUT["DataBase"]["drop_rows"])

else:
    db_path = INPUT["DataBase"]["path_DataBase"]

if INPUT["DataBase"]["colnames_preset"] == 'MetOcean':
    COLNAMES = {
        "datetime": r"datetime (ISO 8601) [UTC]",
        "v_m": r"Wind Speed at 10m (WS_{10}) [m/s]",
        "dir_v_m": r"Wind Direction at 10m (WD_{10}) [\Deg.N-from]",
        "H_s": r"Sign. Wave Height (H_{m0}) [m]",
        "T_p": r"Peak Wave Period (T_{p}) [s]",
        "dir_T_p": r"Peak Wave Direction (PWD) [\Deg.N-from]",
        "H_s_wind": r"Sign. Wave Height - Wind-Sea (H_{m0-Sea}) [m]",
        "T_p_wind": r"Peak Wave Period - Wind-Sea (T_{p-Sea}) [s]",
        "dir_T_mean": r"Mean Wave Direction - Wind-Sea (MWD_{Sea}) [\Deg.N-from]",
        "H_s_swell": r"Sign. Wave Height - Swell (H_{m0-Swell}) [m]",
        "T_p_swell": r"Peak Wave Period - Swell (T_{p-Swell}) [s]",
        "dir_T_mean_Swell": r"Mean Wave Direction - Swell (MWD_{Swell}) [\Deg.N-from]",
        "dir_T_mean_Wind": r"Mean Wave Direction - Wind-Sea (MWD_{Sea}) [\Deg.N-from]",
        "WL": "Water Level (WL) [mMSL]",
        "dir_curr": r"Current Direction (CD) [\Deg.N-to]",
        "WL_tide": r"Water Level - Tide (WL_{Tide}) [mMSL]",
        "v_curr": r"Current Speed (CS) [m/s]",
        "v_curr_tot_5": r"Current Speed - Total at 5% of Water Column (CS_{5%}) [m/s]",
        "v_curr_tide_5": r"Current Speed - Tide at 5% of Water Column (CS_{Tide-5%}) [m/s]",
        "v_curr_tot_50": r"Current Speed - Total at 50% of Water Column (CS_{50%}) [m/s]",
        "v_curr_tide_50": r"Current Speed - Tide at 50% of Water Column (CS_{Tide-50%}) [m/s]",
        "v_curr_tot_100": r"Current Speed - Total at 100% of Water Column (CS_{100%}) [m/s]",
        "v_curr_tide_100": r"Current Speed - Tide at 100% of Water Column (CS_{Tide-100%}) [m/s]"}

elif INPUT["DataBase"]["colnames_preset"] == 'APGMer':

    COLNAMES = {
        "datetime": r"Datetime",
        "v_m": r"Ws,10min,150m",
        "dir_v_m": r"Wdir,1hr,10m",
        "H_max": r"Hmax",
        "H_s": r"Hs",
        "T_p": r"Tp",
        "dir_T_p": r"Pdir",
        "H_s_wind": r"Hs,Sea",
        "T_p_wind": r"Tp,Sea",
        "dir_T_mean": r"Mdir",
        "H_s_swell": r"Hs,Swell",
        "T_p_swell": r"Tp,Swell",
        "dir_T_mean_Swell": r"Mdir,Swell",
        "dir_T_mean_Wind": r"Mdir,Sea",
        "dir_curr": r"Cdir,DA",
        "dir_curr_surf": r"Cdir,Surf",
        "dir_curr_bed": "",
        "dir_curr_tide": r"Cdir,Tide",
        "dir_curr_tide_surf": r"Cdir,Tid,Surf",
        "dir_curr_tide_bed": r"Cdir,Tid,Bed",
        "WL": r"WL",
        "WL_tide": r"WL,Tid",
        "v_curr": r"Cs,DA",
        "v_curr_surf": r"Cs,Surf",
        "v_curr_bed": r"Cs,Bed",
        "v_curr_tide": r"Cs,Tide",
        "v_curr_tide_surf": r"Cs,Tid,Surf",
        "v_curr_tide_bed": r"Cs,Tid,Bed",
    }

elif INPUT["DataBase"]["colnames_preset"] is None:
    COLNAMES = INPUT["ColumNames"]
else:
    print("please choose colnames_preset from 'MetOcen', 'APGMer' or None to import Colnames from Input file")

# %% Calculation
angle_grid, angle_grid_mod = hc_calc.angles(INPUT["AngleSection"]["mode_angle"], INPUT["AngleSection"]["N_angle"], INPUT["AngleSection"]["angle_start"],
                                            width=INPUT["AngleSection"]["width_angle"])

DATA_OUT["VMHS"] = {}
DATA_OUT["HSTP"] = {}
DATA_OUT["VMTP"] = {}
DATA_OUT["RWI"] = {}
DATA_OUT["WaveBreak_Steep"] = {}
DATA_OUT["table_vmhs"] = {}
DATA_OUT["table_vmtp"] = {}
DATA_OUT["AngleDeviation"] = {}
DATA_OUT["Validation"] = {}
DATA_OUT["SensorEval"] = {}
DATA_OUT["Weibull"] = {}
DATA_OUT["ExtremeConture"] = {}


# VMHS
if (('wind' in INPUT["Toggle_Modules"].get("calc_VMHS", {}))
        or ('wind' in INPUT["Toggle_Modules"].get("calc_VMTP", {}))
        or ('wind' in INPUT["Toggle_Modules"].get("calc_Tables", {}))
        or ('wind' in INPUT["Toggle_Modules"].get("calc_Validation", {}))):
    print("calculating VMHS Wind Sea...")

    Input = INPUT["VMHS_wind"]
    table_name = 'Hindcast_combined'
    column_names = [COLNAMES["dir_v_m"], COLNAMES["v_m"], COLNAMES["H_s_wind"]]

    Calc = hc_calc.Calculation()
    Calc.anglecol = COLNAMES["dir_v_m"]

    df = Calc.initilize_from_db(db_path, table_name, column_names, timeframe=timeframe)

    # filter for nans
    indizes_in = Calc.initilize_filter(None, mode='nans')
    df = df.loc[indizes_in]

    directional = hc_calc.calc_VMHS(df[COLNAMES["v_m"]], df[COLNAMES["H_s_wind"]], df[COLNAMES["dir_v_m"]], angle_grid_mod,
                                    N_grid=Input["N_grid"],
                                    deg_reg=Input["deg_reg"],
                                    model_reg=Input["model_reg"],
                                    cut_reg=Input["cut_reg"],
                                    weighting_reg=Input["weighting_reg"],
                                    zone_reg=Input["zone_reg"],
                                    zone_line=Input["zone_line"],
                                    bin_min=Input["bin_min"],
                                    average_correction=Input["average_correction"],
                                    avrg_method=Input["avrg_method"],
                                    make_monotone=Input["make_monotone"])

    omni = hc_calc.calc_VMHS(df[COLNAMES["v_m"]], df[COLNAMES["H_s_wind"]], df[COLNAMES["dir_v_m"]], None,
                             N_grid=Input["N_grid"],
                             deg_reg=Input["deg_reg"],
                             model_reg=Input["model_reg"],
                             cut_reg=Input["cut_reg"],
                             weighting_reg=Input["weighting_reg"],
                             zone_reg=Input["zone_reg"],
                             zone_line=Input["zone_line"],
                             bin_min=Input["bin_min"],
                             average_correction=Input["average_correction"],
                             avrg_method=Input["avrg_method"],
                             make_monotone=Input["make_monotone"])

    Calc.result = omni + directional

    DATA_OUT["VMHS"]["wind"] = Calc

if (('swell' in INPUT["Toggle_Modules"].get("calc_VMHS", {}))
        or ('swell' in INPUT["Toggle_Modules"].get("calc_VMTP", {}))
        or ('swell' in INPUT["Toggle_Modules"].get("calc_Tables", {}))
        or ('swell' in INPUT["Toggle_Modules"].get("calc_Validation", {}))):
    print("calculating VMHS Swell Sea...")

    Input = INPUT["VMHS_swell"]
    table_name = 'Hindcast_combined'
    column_names = [COLNAMES["dir_T_mean_Swell"], COLNAMES["v_m"], COLNAMES["H_s_swell"]]

    Calc = hc_calc.Calculation()
    df = Calc.initilize_from_db(db_path, table_name, column_names, timeframe=timeframe)

    # filter for nans
    indizes_in = Calc.initilize_filter(None, mode='nans')
    df = df.loc[indizes_in]

    directional = hc_calc.calc_VMHS(df[COLNAMES["v_m"]], df[COLNAMES["H_s_swell"]], df[COLNAMES["dir_T_mean_Swell"]], angle_grid_mod,
                                    N_grid=Input["N_grid"],
                                    deg_reg=Input["deg_reg"],
                                    model_reg=Input["model_reg"],
                                    cut_reg=Input["cut_reg"],
                                    weighting_reg=Input["weighting_reg"],
                                    zone_reg=Input["zone_reg"],
                                    zone_line=Input["zone_line"],
                                    bin_min=Input["bin_min"],
                                    average_correction=Input["average_correction"],
                                    avrg_method=Input["avrg_method"],
                                    make_monotone=Input["make_monotone"])

    omni = hc_calc.calc_VMHS(df[COLNAMES["v_m"]], df[COLNAMES["H_s_swell"]], df[COLNAMES["dir_T_mean_Swell"]], None,
                             N_grid=Input["N_grid"],
                             deg_reg=Input["deg_reg"],
                             model_reg=Input["model_reg"],
                             cut_reg=Input["cut_reg"],
                             weighting_reg=Input["weighting_reg"],
                             zone_reg=Input["zone_reg"],
                             zone_line=Input["zone_line"],
                             bin_min=Input["bin_min"],
                             average_correction=Input["average_correction"],
                             avrg_method=Input["avrg_method"],
                             make_monotone=Input["make_monotone"])

    Calc.result = omni + directional
    DATA_OUT["VMHS"]["swell"] = Calc

if (('total' in INPUT["Toggle_Modules"].get("calc_VMHS", {}))
        or ('total' in INPUT["Toggle_Modules"].get("calc_VMTP", {}))
        or ('total' in INPUT["Toggle_Modules"].get("calc_Tables", {}))
        or ('total' in INPUT["Toggle_Modules"].get("calc_Validation", {}))):
    print("calculating VMHS Total Sea...")

    Input = INPUT["VMHS_total"]
    table_name = 'Hindcast_combined'
    column_names = [COLNAMES["dir_T_mean"], COLNAMES["v_m"], COLNAMES["H_s"]]

    Calc = hc_calc.Calculation()
    df = Calc.initilize_from_db(db_path, table_name, column_names, timeframe=timeframe)

    # filter for nans
    indizes_in = Calc.initilize_filter(None, mode='nans')
    df = df.loc[indizes_in]

    directional = hc_calc.calc_VMHS(df[COLNAMES["v_m"]], df[COLNAMES["H_s"]], df[COLNAMES["dir_T_mean"]], angle_grid_mod,
                                    N_grid=Input["N_grid"],
                                    deg_reg=Input["deg_reg"],
                                    model_reg=Input["model_reg"],
                                    cut_reg=Input["cut_reg"],
                                    weighting_reg=Input["weighting_reg"],
                                    zone_reg=Input["zone_reg"],
                                    zone_line=Input["zone_line"],
                                    bin_min=Input["bin_min"],
                                    average_correction=Input["average_correction"],
                                    avrg_method=Input["avrg_method"],
                                    make_monotone=Input["make_monotone"])

    omni = hc_calc.calc_VMHS(df[COLNAMES["v_m"]], df[COLNAMES["H_s"]], df[COLNAMES["dir_T_mean"]], None,
                             N_grid=Input["N_grid"],
                             deg_reg=Input["deg_reg"],
                             model_reg=Input["model_reg"],
                             cut_reg=Input["cut_reg"],
                             weighting_reg=Input["weighting_reg"],
                             zone_reg=Input["zone_reg"],
                             zone_line=Input["zone_line"],
                             bin_min=Input["bin_min"],
                             average_correction=Input["average_correction"],
                             avrg_method=Input["avrg_method"],
                             make_monotone=Input["make_monotone"])

    Calc.result = omni + directional
    DATA_OUT["VMHS"]["total"] = Calc

# HSTP
if (('wind' in INPUT["Toggle_Modules"].get("calc_HSTP", {}))
        or ('wind' in INPUT["Toggle_Modules"].get("calc_VMTP", {}))
        or ('wind' in INPUT["Toggle_Modules"].get("calc_Tables", {}))
        or ('wind' in INPUT["Toggle_Modules"].get("calc_Validation", {}))):
    print("calculating HSTP Wind Sea...")

    Input = INPUT["HSTP_wind"]
    table_name = 'Hindcast_combined'
    column_names = [COLNAMES["H_s_wind"], COLNAMES["T_p_wind"], COLNAMES["dir_v_m"]]

    Calc = hc_calc.Calculation()
    df = Calc.initilize_from_db(db_path, table_name, column_names, timeframe=timeframe)

    # filter for nans
    indizes_in = Calc.initilize_filter(None, mode='nans')
    df = df.loc[indizes_in]

    omni = hc_calc.calc_HSTP(df[COLNAMES["H_s_wind"]], df[COLNAMES["T_p_wind"]], df[COLNAMES["dir_v_m"]], None,
                             N_grid=Input["N_grid"],
                             deg_reg=Input["deg_reg"],
                             model_reg=Input["model_reg"],
                             cut_reg=Input["cut_reg"],
                             weighting_reg=Input["weighting_reg"],
                             zone_reg=Input["zone_reg"],
                             zone_line=Input["zone_line"],
                             bin_min=Input["bin_min"],
                             quantile=Input["quantile"],
                             quant_up=Input["quant_up"],
                             quant_low=Input["quant_low"],
                             percentiles=Input["percentiles"],
                             avrg_method=Input["avrg_method"])

    directional = hc_calc.calc_HSTP(df[COLNAMES["H_s_wind"]], df[COLNAMES["T_p_wind"]], df[COLNAMES["dir_v_m"]], angle_grid_mod,
                                    N_grid=Input["N_grid"],
                                    deg_reg=Input["deg_reg"],
                                    model_reg=Input["model_reg"],
                                    cut_reg=Input["cut_reg"],
                                    weighting_reg=Input["weighting_reg"],
                                    zone_reg=Input["zone_reg"],
                                    zone_line=Input["zone_line"],
                                    bin_min=Input["bin_min"],
                                    quantile=Input["quantile"],
                                    quant_up=Input["quant_up"],
                                    quant_low=Input["quant_low"],
                                    percentiles=Input["percentiles"],
                                    avrg_method=Input["avrg_method"])

    Calc.result = omni + directional

    DATA_OUT["HSTP"]["wind"] = Calc

if (('swell' in INPUT["Toggle_Modules"].get("calc_HSTP", {}))
        or ('swell' in INPUT["Toggle_Modules"].get("calc_VMTP", {}))
        or ('swell' in INPUT["Toggle_Modules"].get("calc_Tables", {}))
        or ('swell' in INPUT["Toggle_Modules"].get("calc_Validation", {}))):
    print("calculating HSTP Swell Sea...")

    Input = INPUT["HSTP_swell"]
    table_name = 'Hindcast_combined'
    column_names = [COLNAMES["H_s_swell"], COLNAMES["T_p_swell"], COLNAMES["dir_T_mean_Swell"]]

    Calc = hc_calc.Calculation()
    df = Calc.initilize_from_db(db_path, table_name, column_names, timeframe=timeframe)

    # filter for nans
    indizes_in = Calc.initilize_filter(None, mode='nans')
    df = df.loc[indizes_in]

    directional = hc_calc.calc_HSTP(df[COLNAMES["H_s_swell"]], df[COLNAMES["T_p_swell"]], df[COLNAMES["dir_T_mean_Swell"]], angle_grid_mod,
                                    N_grid=Input["N_grid"],
                                    deg_reg=Input["deg_reg"],
                                    model_reg=Input["model_reg"],
                                    cut_reg=Input["cut_reg"],
                                    weighting_reg=Input["weighting_reg"],
                                    zone_reg=Input["zone_reg"],
                                    zone_line=Input["zone_line"],
                                    bin_min=Input["bin_min"],
                                    quantile=Input["quantile"],
                                    quant_up=Input["quant_up"],
                                    quant_low=Input["quant_low"],
                                    percentiles=Input["percentiles"],
                                    avrg_method=Input["avrg_method"])

    omni = hc_calc.calc_HSTP(df[COLNAMES["H_s_swell"]], df[COLNAMES["T_p_swell"]], df[COLNAMES["dir_T_mean_Swell"]], None,
                             N_grid=Input["N_grid"],
                             deg_reg=Input["deg_reg"],
                             model_reg=Input["model_reg"],
                             cut_reg=Input["cut_reg"],
                             weighting_reg=Input["weighting_reg"],
                             zone_reg=Input["zone_reg"],
                             zone_line=Input["zone_line"],
                             bin_min=Input["bin_min"],
                             quantile=Input["quantile"],
                             quant_up=Input["quant_up"],
                             quant_low=Input["quant_low"],
                             percentiles=Input["percentiles"],
                             avrg_method=Input["avrg_method"])

    Calc.result = omni + directional
    DATA_OUT["HSTP"]["swell"] = Calc

if (('total' in INPUT["Toggle_Modules"].get("calc_HSTP", {}))
        or ('total' in INPUT["Toggle_Modules"].get("calc_VMTP", {}))
        or ('total' in INPUT["Toggle_Modules"].get("calc_Tables", {}))
        or ('total' in INPUT["Toggle_Modules"].get("calc_Validation", {}))):
    print("calculating HSTP total Sea...")

    Input = INPUT["HSTP_total"]
    table_name = 'Hindcast_combined'
    column_names = [COLNAMES["H_s"], COLNAMES["T_p"], COLNAMES["dir_T_mean"]]

    Calc = hc_calc.Calculation()
    df = Calc.initilize_from_db(db_path, table_name, column_names, timeframe=timeframe)

    # filter for nans
    indizes_in = Calc.initilize_filter(None, mode='nans')
    df = df.loc[indizes_in]

    directional = hc_calc.calc_HSTP(df[COLNAMES["H_s"]], df[COLNAMES["T_p"]], df[COLNAMES["dir_T_mean"]], angle_grid_mod,
                                    N_grid=Input["N_grid"],
                                    deg_reg=Input["deg_reg"],
                                    model_reg=Input["model_reg"],
                                    cut_reg=Input["cut_reg"],
                                    weighting_reg=Input["weighting_reg"],
                                    zone_reg=Input["zone_reg"],
                                    zone_line=Input["zone_line"],
                                    bin_min=Input["bin_min"],
                                    quantile=Input["quantile"],
                                    quant_up=Input["quant_up"],
                                    quant_low=Input["quant_low"],
                                    percentiles=Input["percentiles"],
                                    avrg_method=Input["avrg_method"])

    omni = hc_calc.calc_HSTP(df[COLNAMES["H_s"]], df[COLNAMES["T_p"]], df[COLNAMES["dir_T_mean"]], None,
                             N_grid=Input["N_grid"],
                             deg_reg=Input["deg_reg"],
                             model_reg=Input["model_reg"],
                             cut_reg=Input["cut_reg"],
                             weighting_reg=Input["weighting_reg"],
                             zone_reg=Input["zone_reg"],
                             zone_line=Input["zone_line"],
                             bin_min=Input["bin_min"],
                             quantile=Input["quantile"],
                             quant_up=Input["quant_up"],
                             quant_low=Input["quant_low"],
                             percentiles=Input["percentiles"],
                             avrg_method=Input["avrg_method"])

    Calc.result = omni + directional
    DATA_OUT["HSTP"]["total"] = Calc

# VMTP
if (('wind' in INPUT["Toggle_Modules"].get("calc_VMTP", {}))
        or ('wind' in INPUT["Toggle_Modules"].get("calc_Tables", {}))
        or ('wind' in INPUT["Toggle_Modules"].get("calc_Validation", {}))):
    print("calculating VMTP Wind Sea...")

    table_name = 'Hindcast_combined'
    column_names = [COLNAMES["T_p_wind"], COLNAMES["v_m"], COLNAMES["dir_v_m"]]

    Calc = hc_calc.Calculation()
    Calc.initilize_from_db(db_path, table_name, column_names, timeframe=timeframe)

    Calc.result = hc_calc.calc_VMTP(DATA_OUT["VMHS"]["wind"].result, DATA_OUT["HSTP"]["wind"].result, fill_range=False)
    DATA_OUT["VMTP"]["wind"] = Calc

if (('swell' in INPUT["Toggle_Modules"].get("calc_VMTP", {}))
        or ('swell' in INPUT["Toggle_Modules"].get("calc_Tables", {}))
        or ('swell' in INPUT["Toggle_Modules"].get("calc_Validation", {}))):
    print("calculating VMTP Swell Sea...")

    table_name = 'Hindcast_combined'
    column_names = [COLNAMES["T_p_swell"], COLNAMES["v_m"], COLNAMES["dir_T_mean_Swell"]]

    Calc = hc_calc.Calculation()
    df = Calc.initilize_from_db(db_path, table_name, column_names, timeframe=timeframe)

    Calc.result = hc_calc.calc_VMTP(DATA_OUT["VMHS"]["swell"].result, DATA_OUT["HSTP"]["swell"].result, vm_points=df[COLNAMES["v_m"]], fill_range=True)
    DATA_OUT["VMTP"]["swell"] = Calc

if (('total' in INPUT["Toggle_Modules"].get("calc_VMTP", {}))
        or ('total' in INPUT["Toggle_Modules"].get("calc_Tables", {}))
        or ('total' in INPUT["Toggle_Modules"].get("calc_Validation", {}))):
    print("calculating VMTP total Sea...")

    table_name = 'Hindcast_combined'
    column_names = [COLNAMES["T_p"], COLNAMES["v_m"], COLNAMES["dir_T_mean"]]

    Calc = hc_calc.Calculation()
    Calc.initilize_from_db(db_path, table_name, column_names, timeframe=timeframe)

    Calc.result = hc_calc.calc_VMTP(DATA_OUT["VMHS"]["total"].result, DATA_OUT["HSTP"]["total"].result, fill_range=False)
    DATA_OUT["VMTP"]["total"] = Calc

# Table
if (('wind' in INPUT["Toggle_Modules"].get("calc_Tables", {}))
        or ('wind' in INPUT["Toggle_Modules"].get("calc_Validation", {}))):

    print("calculating Tables Wind Sea...")

    vmhs = DATA_OUT["VMHS"]["wind"]
    vmtp = DATA_OUT["VMTP"]["wind"]

    vm_data = vmhs.load_from_db([COLNAMES["v_m"]])
    vm_data = vm_data[vm_data.keys()[0]]

    Input = INPUT["Tables"]

    vm_zone = Input["vm_zone"]
    if Input["vm_zone"][1] is None:
        vm_zone[1] = max(vm_data.values)
    vm_grid = gl.range_stepfix(Input["vm_step"], vm_zone)

    Calc = hc_calc.Calculation()

    Calc.result = hc_calc.calc_tables(DATA_OUT["VMHS"]["wind"].result, vm_grid, vm_data)
    DATA_OUT["table_vmhs"]["wind"] = Calc

    Calc = hc_calc.Calculation()
    Calc.result = hc_calc.calc_tables(DATA_OUT["VMTP"]["wind"].result, vm_grid, vm_data)
    DATA_OUT["table_vmtp"]["wind"] = Calc

if (('swell' in INPUT["Toggle_Modules"].get("calc_Tables", {}))
        or ('swell' in INPUT["Toggle_Modules"].get("calc_Validation", {}))):

    print("calculating Tables Swell Sea...")

    vmhs = DATA_OUT["VMHS"]["swell"]
    vmtp = DATA_OUT["VMTP"]["swell"]

    vm_data = vmhs.load_from_db([COLNAMES["v_m"]])
    vm_data = vm_data[vm_data.keys()[0]]

    Input = INPUT["Tables"]

    vm_zone = Input["vm_zone"]
    if Input["vm_zone"][1] is None:
        vm_zone[1] = max(vm_data.values)

    vm_grid = gl.range_stepfix(Input["vm_step"], vm_zone)

    Calc = hc_calc.Calculation()
    Calc.result = hc_calc.calc_tables(DATA_OUT["VMHS"]["swell"].result, vm_grid, vm_data)
    DATA_OUT["table_vmhs"]["swell"] = Calc

    Calc = hc_calc.Calculation()
    Calc.result = hc_calc.calc_tables(DATA_OUT["VMTP"]["swell"].result, vm_grid, vm_data)
    DATA_OUT["table_vmtp"]["swell"] = Calc
    # print("calculating Tables Swell Sea...")
    #
    # vmhs = DATA_OUT["VMHS"]["swell"]
    # vmtp = DATA_OUT["VMTP"]["swell"]
    #
    # vm_data = vmhs.load_from_db([COLNAMES["v_m"]])
    # vm_data = vm_data[vm_data.keys()[0]]
    #
    # Input = INPUT["Tables"]
    #
    # vm_zone = Input["vm_zone"]
    # if Input["vm_zone"][1] is None:
    #     vm_zone[1] = max(vm_data.values)
    # vm_grid = gl.range_stepfix(Input["vm_step"], vm_zone)
    #
    # Calc = hc_calc.Calculation()
    #
    # Calc.result = hc_calc.calc_tables(DATA_OUT["VMHS"]["swell"].result, vm_grid, vm_data)
    # DATA_OUT["table_vmhs"]["swell"] = Calc
    #
    # Calc = hc_calc.Calculation()
    # Calc.result = hc_calc.calc_tables(DATA_OUT["VMTP"]["swell"].result, vm_grid, vm_data)
    # DATA_OUT["table_vmtp"]["swell"] = Calc

if (('total' in INPUT["Toggle_Modules"].get("calc_Tables", {}))
        or ('total' in INPUT["Toggle_Modules"].get("calc_Validation", {}))):
    print("calculating Tables Total Sea...")

    vmhs = DATA_OUT["VMHS"]["total"]
    vmtp = DATA_OUT["VMTP"]["total"]

    vm_data = vmhs.load_from_db([COLNAMES["v_m"]])
    vm_data = vm_data[vm_data.keys()[0]]

    Input = INPUT["Tables"]

    vm_zone = Input["vm_zone"]
    if Input["vm_zone"][1] is None:
        vm_zone[1] = max(vm_data.values)

    vm_grid = gl.range_stepfix(Input["vm_step"], vm_zone)

    Calc = hc_calc.Calculation()
    Calc.result = hc_calc.calc_tables(DATA_OUT["VMHS"]["total"].result, vm_grid, vm_data)
    DATA_OUT["table_vmhs"]["total"] = Calc

    Calc = hc_calc.Calculation()
    Calc.result = hc_calc.calc_tables(DATA_OUT["VMTP"]["total"].result, vm_grid, vm_data)
    DATA_OUT["table_vmtp"]["total"] = Calc
    # print("calculating Tables total Sea...")
    #
    # vmhs = DATA_OUT["VMHS"]["total"]
    # vmtp = DATA_OUT["VMTP"]["total"]
    #
    # vm_data = vmhs.load_from_db([COLNAMES["v_m"]])
    # vm_data = vm_data[vm_data.keys()[0]]
    #
    # Input = INPUT["Tables"]
    #
    # vm_zone = Input["vm_zone"]
    # if Input["vm_zone"][1] is None:
    #     vm_zone[1] = max(vm_data.values)
    # vm_grid = gl.range_stepfix(Input["vm_step"], vm_zone)
    #
    # Calc = hc_calc.Calculation()
    #
    # Calc.result = hc_calc.calc_tables(DATA_OUT["VMHS"]["total"].result, vm_grid, vm_data)
    # DATA_OUT["table_vmhs"]["total"] = Calc
    #
    # Calc = hc_calc.Calculation()
    # Calc.result = hc_calc.calc_tables(DATA_OUT["VMTP"]["total"].result, vm_grid, vm_data)
    # DATA_OUT["table_vmtp"]["total"] = Calc

# RWI
if 'wind' in INPUT["Toggle_Modules"].get("calc_RWI", {}):
    print("calculating RWI Wind Sea...")

    table_name = 'Hindcast_combined'
    column_names = [COLNAMES["H_s_wind"], COLNAMES["T_p_wind"], COLNAMES["dir_v_m"]]

    Calc = hc_calc.Calculation()
    df = Calc.initilize_from_db(db_path, table_name, column_names, timeframe=timeframe)

    # filter for nans
    indizes_in = Calc.initilize_filter(None, mode='nans')
    df = df.loc[indizes_in]

    directional, RWI_max = hc_calc.calc_RWI(df[COLNAMES["H_s_wind"]],
                                            df[COLNAMES["T_p_wind"]],
                                            df[COLNAMES["dir_v_m"]],
                                            angle_grid_mod,
                                            INPUT["Structure"]["f_0"])

    omni, RWI_max = hc_calc.calc_RWI(df[COLNAMES["H_s_wind"]],
                                     df[COLNAMES["T_p_wind"]],
                                     df[COLNAMES["dir_v_m"]],
                                     None,
                                     INPUT["Structure"]["f_0"])

    Calc.result = omni + directional

    DATA_OUT["RWI"]["wind"] = Calc

if 'total' in INPUT["Toggle_Modules"].get("calc_RWI", {}):
    print("calculating RWI total Sea...")

    table_name = 'Hindcast_combined'
    column_names = [COLNAMES["H_s"], COLNAMES["T_p"], COLNAMES["dir_v_m"]]

    Calc = hc_calc.Calculation()
    df = Calc.initilize_from_db(db_path, table_name, column_names, timeframe=timeframe)

    # filter for nans
    indizes_in = Calc.initilize_filter(None, mode='nans')
    df = df.loc[indizes_in]

    directional, RWI_max = hc_calc.calc_RWI(df[COLNAMES["H_s"]],
                                            df[COLNAMES["T_p"]],
                                            df[COLNAMES["dir_v_m"]],
                                            angle_grid_mod,
                                            INPUT["Structure"]["f_0"])

    omni, RWI_max = hc_calc.calc_RWI(df[COLNAMES["H_s"]],
                                     df[COLNAMES["T_p"]],
                                     df[COLNAMES["dir_v_m"]],
                                     None,
                                     INPUT["Structure"]["f_0"])

    Calc.result = omni + directional

    DATA_OUT["RWI"]["total"] = Calc

# WaveBreak
if 'wind' in INPUT["Toggle_Modules"].get("calc_WaveBreak_Steep", {}):
    print("calculating WaveBreakSteep Wind Sea...")

    Input = INPUT["Structure"]
    table_name = 'Hindcast_combined'
    column_names = [COLNAMES["dir_v_m"], COLNAMES["H_s_wind"], COLNAMES["T_p_wind"]]

    Calc = hc_calc.Calculation()
    df = Calc.initilize_from_db(db_path, table_name, column_names, timeframe=timeframe)

    # filter for nans
    indizes_in = Calc.initilize_filter(None, mode='nans')
    df = df.loc[indizes_in]

    directional = hc_calc.calc_WaveBreak_Steep(df[COLNAMES["H_s_wind"]],
                                               df[COLNAMES["T_p_wind"]],
                                               df[COLNAMES["dir_v_m"]],
                                               angle_grid_mod,
                                               Input["steep_crit"],
                                               Input["d"])

    omni = hc_calc.calc_WaveBreak_Steep(df[COLNAMES["H_s_wind"]],
                                        df[COLNAMES["T_p_wind"]],
                                        df[COLNAMES["dir_v_m"]],
                                        None,
                                        Input["steep_crit"],
                                        Input["d"])

    Calc.result = omni + directional

    DATA_OUT["WaveBreak_Steep"]["wind"] = Calc

if 'total' in INPUT["Toggle_Modules"].get("calc_WaveBreak_Steep", {}):
    print("calculating WaveBreakSteep total Sea...")

    Input = INPUT["Structure"]
    table_name = 'Hindcast_combined'
    column_names = [COLNAMES["dir_v_m"], COLNAMES["H_s"], COLNAMES["T_p"]]

    Calc = hc_calc.Calculation()
    df = Calc.initilize_from_db(db_path, table_name, column_names, timeframe=timeframe)

    # filter for nans
    indizes_in = Calc.initilize_filter(None, mode='nans')
    df = df.loc[indizes_in]

    directional = hc_calc.calc_WaveBreak_Steep(df[COLNAMES["H_s"]],
                                               df[COLNAMES["T_p"]],
                                               df[COLNAMES["dir_v_m"]],
                                               angle_grid_mod,
                                               Input["steep_crit"],
                                               Input["d"])

    omni = hc_calc.calc_WaveBreak_Steep(df[COLNAMES["H_s"]],
                                        df[COLNAMES["T_p"]],
                                        df[COLNAMES["dir_v_m"]],
                                        None,
                                        Input["steep_crit"],
                                        Input["d"])

    Calc.result = omni + directional

    DATA_OUT["WaveBreak_Steep"]["total"] = Calc

# Angle deviation
if INPUT["Toggle_Modules"].get("calc_AngleDeviation", {}):
    print("calculating Angle Deviation...")

    Input = INPUT["AngleDeviation"]
    table_name = 'Hindcast_combined'
    column_names = [COLNAMES[INPUT["Toggle_Modules"]["calc_AngleDeviation"][0]],
                    COLNAMES[INPUT["Toggle_Modules"]["calc_AngleDeviation"][1]],
                    COLNAMES[INPUT["Toggle_Modules"]["calc_AngleDeviation"][2]]]

    Calc = hc_calc.Calculation()
    df = Calc.initilize_from_db(db_path, table_name, column_names, timeframe=timeframe)

    # filter for nans
    indizes_in = Calc.initilize_filter(None, mode='nans')
    df = df.loc[indizes_in]

    if INPUT["AngleDeviation"]["filter_by"] is not None:
        cols = [COLNAMES[curr] for curr in INPUT["AngleDeviation"]["filter_by"]]
        indizes_in = Calc.initilize_filter(cols, INPUT["AngleDeviation"]["margin"])
        df = df[df.index.isin(indizes_in)]

    # omni, only comparison of global angles
    diff, diff_abs = gl.angle_deviation(df[column_names[0]], df[column_names[1]])

    diff_rolling_mean, bin_edges, _ = sc.stats.binned_statistic(df[column_names[0]], diff_abs, statistic='mean', bins=np.linspace(0, 360, 200))
    diff_rolling_mean = pd.Series(diff_rolling_mean, index=(bin_edges[1:] + bin_edges[:-1]) / 2)

    omni = [hc_calc.Segment(0,
                            angles=None,
                            result={"points": pd.DataFrame({"diff": diff, "diff_abs": diff_abs}, index=df.index), "mean": diff_rolling_mean},
                            colnames={'ang_orig': column_names[0], 'ang_comp': column_names[1]},
                            angle_name=None,
                            indizes=list(df.index))]

    directional = hc_calc.calc_angle_deviation_tables(df[column_names[0]],
                                                      df[column_names[1]],
                                                      df[column_names[2]],
                                                      angle_grid_mod,
                                                      v_m_zone=Input["v_m_zone"],
                                                      v_m_step=Input["v_m_step"],
                                                      N_angle_comp_sec=Input["N_angle_comp_sec"])

    Calc.result = omni + directional

    DATA_OUT["AngleDeviation"] = Calc

# Roseplot
if len(INPUT["Toggle_Modules"].get("calc_Roseplots", {})) > 0:
    print("calculating Roseplots...")
    DATA_OUT["Roseplot"] = {}
    for Roseplot_cols in INPUT["Toggle_Modules"]["calc_Roseplots"]:
        table_name = 'Hindcast_combined'
        column_names = [COLNAMES[col] for col in Roseplot_cols]

        Calc = hc_calc.Calculation()
        df = Calc.initilize_from_db(db_path, table_name, column_names, timeframe=timeframe)

        # filter for nans
        indizes_in = Calc.initilize_filter(None, mode='nans')
        df = df.loc[indizes_in]

        temp, bins = hc_calc.calc_Roseplot(df[column_names[0]], df[column_names[1]], angle_grid_mod)

        Calc.result = {"table": temp, "r_bins": bins, "r_max": np.max(df[column_names[1]])}

        DATA_OUT["Roseplot"][f"{Roseplot_cols[0]} over {Roseplot_cols[1]}"] = Calc

# Extreme Values
if len(INPUT["Toggle_Modules"].get("calc_ExtremeValues", {})) > 0:
    print("calculating ExtremeValues...")

    DATA_OUT["ExtremeValues"] = {}
    Input = INPUT["ExtremeValues"]
    for cols in INPUT["Toggle_Modules"]["calc_ExtremeValues"]:
        table_name = 'Hindcast_combined'
        column_names = [COLNAMES[col] for col in cols]

        Calc = hc_calc.Calculation()
        df = Calc.initilize_from_db(db_path, table_name, column_names, timeframe=timeframe)

        directional = hc_calc.calc_ExtemeValues(df[column_names[1]],
                                                df[column_names[0]],
                                                angle_grid_mod,
                                                T_return_single=Input["T_return"],
                                                conf_inter_mode=Input["conf_inter_mode"],
                                                conf_inter_algorithm=Input["conf_inter_algorithm"],
                                                N_itter=Input["N_itter"],
                                                freq_samp=Input["freq_samp"],
                                                perc_up=Input["perc_up"],
                                                perc_down=Input["perc_down"],
                                                time_window_offset=Input["time_window_offset"]
                                                )

        omni = hc_calc.calc_ExtemeValues(df[column_names[1]],
                                         df[column_names[0]],
                                         None,
                                         T_return_single=Input["T_return"],
                                         conf_inter_mode=Input["conf_inter_mode"],
                                         conf_inter_algorithm=Input["conf_inter_algorithm"],
                                         N_itter=Input["N_itter"],
                                         freq_samp=Input["freq_samp"],
                                         perc_up=Input["perc_up"],
                                         perc_down=Input["perc_down"],
                                         time_window_offset=Input["time_window_offset"]
                                         )

        Calc.result = omni + directional

        DATA_OUT["ExtremeValues"][f"{cols[0]} over {cols[1]}"] = Calc

# Extreme Conture Plots
if len(INPUT["Toggle_Modules"].get("calc_ExtremeConture", {})) > 0:
    print("calculating Extreme Conture Plots...")
    DATA_OUT["ExtremeConturePlots"] = {}
    for cols in INPUT["Toggle_Modules"]["calc_ExtremeConture"]:
        table_name = 'Hindcast_combined'
        column_names = [COLNAMES[col] for col in cols]

        Calc = hc_calc.Calculation()
        df = Calc.initilize_from_db(db_path, table_name, column_names, timeframe=timeframe)

        # filter for nans
        indizes_in = Calc.initilize_filter(None, mode='nans')
        df = df.loc[indizes_in]

        out_direc = hc_calc.calc_extreme_contures(df[column_names[0]],
                                            df[column_names[1]],
                                            df[column_names[2]],
                                            angle_grid,
                                            INPUT["ExtremeValues"]["T_return"])

        out_omni = hc_calc.calc_extreme_contures(df[column_names[0]],
                                            df[column_names[1]],
                                            df[column_names[2]],
                                            None,
                                            INPUT["ExtremeValues"]["T_return"])

        Calc.result = out_omni + out_direc
        DATA_OUT["ExtremeConture"][f"{cols[0]} over {cols[1]}"] = Calc



# Validation
if 'wind' in INPUT["Toggle_Modules"].get("calc_Validation", {}):
    print("calculating Validation wind...")

    Input = INPUT["Validation_wind"]
    table_name = 'Hindcast_combined'
    column_names = [COLNAMES["H_s_wind"], COLNAMES["T_p_wind"], COLNAMES["dir_v_m"], COLNAMES["v_m"]]

    Calc = hc_calc.Calculation()

    df_data = gl.export_df_from_sql(db_path, table_name, column_names=column_names, timeframe=timeframe)
    df_data['gamma'] = 3.3

    df_data = df_data.dropna(how='any')

    JBOOST_proj_Path = INPUT['DataBase']['JBOOST_proj']
    JBOOST_proj_input_path = INPUT['DataBase']['JBOOST_input']
    JBOOST_exe_path = path_main + '\\JBOOST\\'

    if not Input["from_DB"]:
        print(f"   looking in database: {db_path} for current nodes and timeframe")
        table_name_DEL, colnames_DEL = hc_calc.update_DEL_db(db_path,
                                                             df_data[column_names[0]],
                                                             df_data[column_names[1]],
                                                             df_data['gamma'],
                                                             proj_path=JBOOST_proj_Path,
                                                             input_path=JBOOST_proj_input_path,
                                                             exe_path=JBOOST_exe_path)

        colnames_DEL = gl.find_string_with_substrings(colnames_DEL, Input['nodes_to_load'])

        df = Calc.initilize_from_db(db_path, table_name_DEL, colnames_DEL, timeframe=timeframe, indizes=df_data.index)

    else:
        print(f"   loading from database: {db_path} in table {Input['table_name']} for current nodes and timeframe")

        df = gl.export_df_from_sql(db_path, Input['table_name'], timeframe=timeframe, indizes=df_data.index)
        df = gl.filter_df_cols_by_keywords(df, Input['nodes_to_load'])

        Calc.basedata = {"dbname": db_path,
                         "tablename": Input['table_name'],
                         "colnames_ini": df.keys,
                         "db_timeframe": [df.index[0], df.index[1]],
                         "N_rows": len(df),
                         "sample_rate": gl.median_sample_rate(df.index),
                         "indizes": df.index}

    Calc.initilize_filter(None, mode='nans')

    print(f"   processing calculated/loaded DEL data and comparing to condensed data in tables")
    result = hc_calc.calc_Validation(df,
                                     df_data[column_names[3]],
                                     df_data[column_names[2]],
                                     DATA_OUT["table_vmhs"]["wind"].result,
                                     DATA_OUT["table_vmtp"]["wind"].result,
                                     INPUT["DataBase"]["JBOOST_proj"],
                                     INPUT["DataBase"]["JBOOST_input"],
                                     r".\\JBOOST\\")

    Calc.result = result

    DATA_OUT["Validation"]["wind"] = Calc

if 'swell' in INPUT["Toggle_Modules"].get("calc_Validation", {}):
    print("calculating Validation swell...")

    Input = INPUT["Validation_swell"]
    table_name = 'Hindcast_combined'
    column_names = [COLNAMES["H_s_swell"], COLNAMES["T_p_swell"], COLNAMES["dir_T_mean_Swell"], COLNAMES["v_m"]]

    Calc = hc_calc.Calculation()
    df_data = gl.export_df_from_sql(db_path, table_name, column_names=column_names, timeframe=timeframe)
    df_data['gamma'] = 3.3
    df_data = df_data.dropna(how='any')

    JBOOST_proj_Path = INPUT['DataBase']['JBOOST_proj']
    JBOOST_proj_input_path = INPUT['DataBase']['JBOOST_input']
    JBOOST_exe_path = path_main + '\\JBOOST\\'

    if not Input["from_DB"]:
        print(f"   looking in database: {db_path} for current nodes and timeframe")
        table_name_DEL, colnames_DEL = hc_calc.update_DEL_db(db_path,
                                                             df_data[column_names[0]],
                                                             df_data[column_names[1]],
                                                             df_data['gamma'],
                                                             proj_path=JBOOST_proj_Path,
                                                             input_path=JBOOST_proj_input_path,
                                                             exe_path=JBOOST_exe_path)

        colnames_DEL = gl.find_string_with_substrings(colnames_DEL, Input['nodes_to_load'])

        df = Calc.initilize_from_db(db_path, table_name_DEL, colnames_DEL, timeframe=timeframe, indizes=df_data.index)

    else:
        print(f"   loading from database: {db_path} in table {Input['table_name']} for current nodes and timeframe")

        df = gl.export_df_from_sql(db_path, Input['table_name'], timeframe=timeframe)

        indizes_in = df.index.intersection(df_data.index)
        indizes_in_2 = df_data.index.intersection(df.index)

        if len(df.index.difference(df_data.index)) != 0:
            print(f"   loaded data from DEL table and from Hindcast table have different indizes, "
                  f"{len(df.index.difference(df_data.index))} points dropped from DEL data, {len(df_data.index.difference(df.index))} dropped from Hincast data")

        df = df.loc[indizes_in]
        df_data = df_data.loc[indizes_in]

        df = gl.filter_df_cols_by_keywords(df, Input['nodes_to_load'])

        Calc.basedata = {"dbname": db_path,
                         "tablename": Input['table_name'],
                         "colnames_ini": df.keys,
                         "db_timeframe": [df.index[0], df.index[1]],
                         "N_rows": len(df),
                         "sample_rate": df.index[1] - df.index[0],
                         "indizes": df.index}

    print(f"   processing calculated/loaded DEL data and comparing to condensed data in tables")
    result = hc_calc.calc_Validation(df,
                                     df_data[column_names[3]],
                                     df_data[column_names[2]],
                                     DATA_OUT["table_vmhs"]["swell"].result,
                                     DATA_OUT["table_vmtp"]["swell"].result,
                                     INPUT["DataBase"]["JBOOST_proj"],
                                     INPUT["DataBase"]["JBOOST_input"],
                                     r".\\JBOOST\\")

    Calc.result = result

    DATA_OUT["Validation"]["swell"] = Calc

# SensorEval
if len(INPUT["Toggle_Modules"].get("calc_SensorEval", {})) > 0:
    print("calculating Sensor Evaluation...")

    for colname in INPUT["Toggle_Modules"]["calc_SensorEval"]:
        table_name = 'Hindcast_combined'
        colname_data = COLNAMES[colname]

        Calc = hc_calc.Calculation()
        df = Calc.initilize_from_db(db_path, table_name, [colname_data], timeframe=timeframe)

        # filter for nans
        indizes_in = Calc.initilize_filter(None, mode='nans')
        df = df.loc[indizes_in]

        #  directional = hc_calc.calc_histogram(df[column_names[0]],
        #                                      df[column_names[1]],
        #                                        angle_grid_mod)

        omni = hc_calc.calc_histogram(df[colname_data],
                                      None,
                                      None)

        Calc.result = omni

        DATA_OUT["SensorEval"][f"{colname}"] = Calc

# Weibull
if len(INPUT["Toggle_Modules"].get("calc_Weibull", {})) > 0:
    print("calculating Weibull fit...")

    for colnames in INPUT["Toggle_Modules"]["calc_Weibull"]:
        table_name = 'Hindcast_combined'
        column_names = [COLNAMES[colnames[0]], COLNAMES[colnames[1]]]

        Calc = hc_calc.Calculation()
        df = Calc.initilize_from_db(db_path, table_name, column_names, timeframe=timeframe)

        # filter for nans
        indizes_in = Calc.initilize_filter(None, mode='nans')
        df = df.loc[indizes_in]

        directional = hc_calc.calc_weibull(df[column_names[1]],
                                           df[column_names[0]],
                                           angle_grid_mod)

        omni = hc_calc.calc_weibull(df[column_names[1]],
                                    None,
                                    None)

        Calc.result = directional + omni

        DATA_OUT["Weibull"][f"{colnames[1]} over {colnames[0]}"] = Calc
# %% Plot
figsize_fullpage = [size * 0.39370079 for size in INPUT["Toggle_Modules"].get("writing_box", {})]
figsize_halfpage = [figsize_fullpage[0], figsize_fullpage[1] / 2]

if 'wind' in INPUT["Toggle_Modules"].get("plot_VMHS", {}):
    print('plotting VMHS wind...')

    if not ('wind' in INPUT["Toggle_Modules"].get("calc_VMHS", {})):
        print("   please toggle calculation to plot")

    else:

        Calc = DATA_OUT["VMHS"]["wind"]

        Tiles = []
        Tiles_omni = []

        df = Calc.load_from_db(colnames_ini=True)
        titels = Calc.create_segment_title()
        titels = gl.alias(titels, COLNAMES, INPUT["Aliase"])

        for i, Seg in enumerate(Calc.result):

            Seg.indizes = pd.to_datetime(Seg.indizes)
            point_data = df[df.index.isin(Seg.indizes)]

            tile_curr = hc_plt.Tile(i,
                                    x_label=gl.alias(Seg.colnames['x'], COLNAMES, INPUT["Aliase"]),
                                    y_label=gl.alias(Seg.colnames['y'], COLNAMES, INPUT["Aliase"]),
                                    title=titels[i])

            Line_mean = hc_plt.Line(x=Seg.result["x"],
                                    y=Seg.result["mean result plot"],
                                    label='mean',
                                    color='black')

            tile_curr.add_line(Line_mean)

            scatter = hc_plt.Scatter(x=point_data[Seg.colnames["x"]],
                                     y=point_data[Seg.colnames["y"]],
                                     cmap='cool',
                                     size=2,
                                     cmap_norm='sqrt')

            tile_curr.add_scatter(scatter)

            if Seg.angles is not None:
                Tiles.append(tile_curr)
            else:
                Tiles_omni.append(tile_curr)

        FIG_direc = hc_plt.plot_tiled(Tiles, global_max=['auto', 'auto'], global_min=[0, 0], grid=[3, 2], figsize=figsize_fullpage)

        FIG_omni = hc_plt.plot_tiled(Tiles_omni, global_max=['auto', 'auto'], global_min=[0, 0], grid=[1, 1], figsize=figsize_halfpage)

        if 'png' in INPUT["Toggle_Modules"]["plot_as"]:
            gl.save_figs_as_png(FIG_direc + FIG_omni, path_out + 'VMHS_wind', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

        if 'pdf' in INPUT["Toggle_Modules"]["plot_as"]:
            gl.save_figs_as_pdf(FIG_direc + FIG_omni, path_out + 'VMHS_wind', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

if 'swell' in INPUT["Toggle_Modules"].get("plot_VMHS", {}):
    print('plotting VMHS swell...')

    if not ('swell' in INPUT["Toggle_Modules"].get("calc_VMHS", {})):
        print("   please toggle calculation to plot")

    else:
        Calc = DATA_OUT["VMHS"]["swell"]

        Tiles = []
        Tiles_omni = []

        df = Calc.load_from_db(colnames_ini=True)
        titels = Calc.create_segment_title()
        titels = gl.alias(titels, COLNAMES, INPUT["Aliase"])

        for i, Seg in enumerate(Calc.result):

            Seg.indizes = pd.to_datetime(Seg.indizes)
            point_data = df[df.index.isin(Seg.indizes)]

            tile_curr = hc_plt.Tile(i,
                                    x_label=gl.alias(Seg.colnames['x'], COLNAMES, INPUT["Aliase"]),
                                    y_label=gl.alias(Seg.colnames['y'], COLNAMES, INPUT["Aliase"]),
                                    title=titels[i])

            Line_mean = hc_plt.Line(x=Seg.result["x"],
                                    y=Seg.result["mean result plot"],
                                    label='mean',
                                    color='black')

            tile_curr.add_line(Line_mean)

            scatter = hc_plt.Scatter(x=point_data[Seg.colnames["x"]],
                                     y=point_data[Seg.colnames["y"]],
                                     cmap='cool',
                                     size=2,
                                     cmap_norm='sqrt')

            tile_curr.add_scatter(scatter)

            if Seg.angles is not None:
                Tiles.append(tile_curr)
            else:
                Tiles_omni.append(tile_curr)

        FIG_direc = hc_plt.plot_tiled(Tiles, global_max=['auto', 'auto'], global_min=[0, 0], grid=[3, 2], figsize=figsize_fullpage)

        FIG_omni = hc_plt.plot_tiled(Tiles_omni, global_max=['auto', 'auto'], global_min=[0, 0], grid=[1, 1], figsize=figsize_halfpage)

        if 'png' in INPUT["Toggle_Modules"]["plot_as"]:
            gl.save_figs_as_png(FIG_direc + FIG_omni, path_out + 'VMHS_swell', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

        if 'pdf' in INPUT["Toggle_Modules"]["plot_as"]:
            gl.save_figs_as_pdf(FIG_direc + FIG_omni, path_out + 'VMHS_swell', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

if 'total' in INPUT["Toggle_Modules"].get("plot_VMHS", {}):
    print('plotting VMHS total...')

    if not ('total' in INPUT["Toggle_Modules"].get("calc_VMHS", {})):
        print("   please toggle calculation to plot")

    else:
        Calc = DATA_OUT["VMHS"]["total"]

        Tiles = []
        Tiles_omni = []

        df = Calc.load_from_db(colnames_ini=True)
        titels = Calc.create_segment_title()
        titels = gl.alias(titels, COLNAMES, INPUT["Aliase"])

        for i, Seg in enumerate(Calc.result):

            Seg.indizes = pd.to_datetime(Seg.indizes)
            point_data = df[df.index.isin(Seg.indizes)]

            tile_curr = hc_plt.Tile(i,
                                    x_label=gl.alias(Seg.colnames['x'], COLNAMES, INPUT["Aliase"]),
                                    y_label=gl.alias(Seg.colnames['y'], COLNAMES, INPUT["Aliase"]),
                                    title=titels[i])

            Line_mean = hc_plt.Line(x=Seg.result["x"],
                                    y=Seg.result["mean result plot"],
                                    label='mean',
                                    color='black')

            tile_curr.add_line(Line_mean)

            scatter = hc_plt.Scatter(x=point_data[Seg.colnames["x"]],
                                     y=point_data[Seg.colnames["y"]],
                                     cmap='cool',
                                     size=2,
                                     cmap_norm='sqrt')

            tile_curr.add_scatter(scatter)

            if Seg.angles is not None:
                Tiles.append(tile_curr)
            else:
                Tiles_omni.append(tile_curr)

        FIG_direc = hc_plt.plot_tiled(Tiles, global_max=['auto', 'auto'], global_min=[0, 0], grid=[3, 2], figsize=figsize_fullpage)

        FIG_omni = hc_plt.plot_tiled(Tiles_omni, global_max=['auto', 'auto'], global_min=[0, 0], grid=[1, 1], figsize=figsize_halfpage)

        if 'png' in INPUT["Toggle_Modules"]["plot_as"]:
            gl.save_figs_as_png(FIG_direc + FIG_omni, path_out + 'VMHS_total', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

        if 'pdf' in INPUT["Toggle_Modules"]["plot_as"]:
            gl.save_figs_as_pdf(FIG_direc + FIG_omni, path_out + 'VMHS_total', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

if 'wind' in INPUT["Toggle_Modules"].get("plot_HSTP", {}):
    print('plotting HSTP wind...')

    if not ('wind' in INPUT["Toggle_Modules"].get("calc_HSTP", {})):
        print("   please toggle calculation to plot")

    else:
        Calc = DATA_OUT["HSTP"]["wind"]
        Input = INPUT["HSTP_wind"]
        Tiles = []
        Tiles_omni = []

        df = Calc.load_from_db(colnames_ini=True)
        titels = Calc.create_segment_title()
        titels = gl.alias(titels, COLNAMES, INPUT["Aliase"])

        for i, Seg in enumerate(Calc.result):

            Seg.indizes = pd.to_datetime(Seg.indizes)
            point_data = df[df.index.isin(Seg.indizes)]

            tile_curr = hc_plt.Tile(i,
                                    x_label=gl.alias(Seg.colnames['x'], COLNAMES, INPUT["Aliase"]),
                                    y_label=gl.alias(Seg.colnames['y'], COLNAMES, INPUT["Aliase"]),
                                    title=titels[i])

            key_plot = [key for key in Seg.result.columns if 'plot' in key]
            key_mean = [key for key in key_plot if 'mean' in key]
            key_perc = [key for key in key_plot if 'percentile' in key]
            key_quantile = [key for key in Seg.result.columns if 'quantile' in key]

            Line_mean = hc_plt.Line(x=Seg.result["x"],
                                    y=Seg.result["mean result plot"],
                                    label='mean',
                                    color='black',
                                    linestyle='-')
            tile_curr.add_line(Line_mean)

            Line_perc_low = hc_plt.Line(x=Seg.result["x"],
                                        y=Seg.result[key_perc[0]],
                                        label=key_perc[0].replace('result', '').replace('plot', ''),
                                        color='black',
                                        linestyle='--')
            tile_curr.add_line(Line_perc_low)

            Line_perc_up = hc_plt.Line(x=Seg.result["x"],
                                       y=Seg.result[key_perc[1]],
                                       label=key_perc[1].replace('result', '').replace('plot', ''),
                                       color='black',
                                       linestyle=':')
            tile_curr.add_line(Line_perc_up)

            if len(key_quantile) > 0:
                Line_quant = hc_plt.Line(x=Seg.result["x"],
                                         y=Seg.result[key_quantile[0]],
                                         label='selected correlation',
                                         color='red',
                                         linestyle='-')
                tile_curr.add_line(Line_quant)

                Line_quant_up = hc_plt.Line(x=None,
                                            y=[1 / Input["quant_up"]],
                                            label=f'upper quantile, f={Input["quant_up"]} Hz',
                                            color='green',
                                            linestyle=':')
                tile_curr.add_line(Line_quant_up)

                Line_quant_low = hc_plt.Line(x=None,
                                             y=[1 / Input["quant_low"]],
                                             label=f'lower quantile, f={Input["quant_low"]} Hz',
                                             color='green',
                                             linestyle='--')
                tile_curr.add_line(Line_quant_low)

            scatter = hc_plt.Scatter(x=point_data[Seg.colnames["x"]],
                                     y=point_data[Seg.colnames["y"]],
                                     cmap='cool',
                                     size=2,
                                     cmap_norm='sqrt')

            tile_curr.add_scatter(scatter)

            if Seg.angles is not None:
                tile_curr.legend_fontsize = 6
                Tiles.append(tile_curr)

            else:
                Tiles_omni.append(tile_curr)

        FIG_direc = hc_plt.plot_tiled(Tiles, global_max=['auto', 'auto'], global_min=[0, 0], grid=[3, 2], figsize=figsize_fullpage)

        FIG_omni = hc_plt.plot_tiled(Tiles_omni, global_max=['auto', 'auto'], global_min=[0, 0], grid=[1, 1], figsize=figsize_halfpage)

        if 'png' in INPUT["Toggle_Modules"]["plot_as"]:
            gl.save_figs_as_png(FIG_direc + FIG_omni, path_out + 'HSTP_wind', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

        if 'pdf' in INPUT["Toggle_Modules"]["plot_as"]:
            gl.save_figs_as_pdf(FIG_direc + FIG_omni, path_out + 'HSTP_wind', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

if 'swell' in INPUT["Toggle_Modules"].get("plot_HSTP", {}):
    print('plotting HSTP swell...')

    if not ('swell' in INPUT["Toggle_Modules"].get("calc_HSTP", {})):
        print("   please toggle calculation to plot")

    else:

        Calc = DATA_OUT["HSTP"]["swell"]
        Input = INPUT["HSTP_swell"]
        Tiles = []
        Tiles_omni = []

        df = Calc.load_from_db(colnames_ini=True)
        titels = Calc.create_segment_title()
        titels = gl.alias(titels, COLNAMES, INPUT["Aliase"])

        for i, Seg in enumerate(Calc.result):

            Seg.indizes = pd.to_datetime(Seg.indizes)
            point_data = df[df.index.isin(Seg.indizes)]

            tile_curr = hc_plt.Tile(i,
                                    x_label=gl.alias(Seg.colnames['x'], COLNAMES, INPUT["Aliase"]),
                                    y_label=gl.alias(Seg.colnames['y'], COLNAMES, INPUT["Aliase"]),
                                    title=titels[i])

            key_plot = [key for key in Seg.result.columns if 'plot' in key]
            key_mean = [key for key in key_plot if 'mean' in key]
            key_perc = [key for key in key_plot if 'percentile' in key]
            key_quantile = [key for key in Seg.result.columns if 'quantile' in key]

            Line_mean = hc_plt.Line(x=Seg.result["x"],
                                    y=Seg.result["mean result plot"],
                                    label='mean',
                                    color='black',
                                    linestyle='-')
            tile_curr.add_line(Line_mean)

            # Line_perc_low = hc_plt.Line(x=Seg.result["x"],
            #                             y=Seg.result[key_perc[0]],
            #                             label=key_perc[0].replace('result', '').replace('plot', ''),
            #                             color='black',
            #                             linestyle='--')
            # tile_curr.add_line(Line_perc_low)

            # Line_perc_up = hc_plt.Line(x=Seg.result["x"],
            #                            y=Seg.result[key_perc[1]],
            #                            label='mean',
            #                            color='black',
            #                            linestyle=':')
            # tile_curr.add_line(Line_perc_up)

            if len(key_quantile) > 0:
                Line_quant = hc_plt.Line(x=Seg.result["x"],
                                         y=Seg.result[key_quantile[0]],
                                         label='selected correlation',
                                         color='red',
                                         linestyle='-')
                tile_curr.add_line(Line_quant)

                Line_quant_up = hc_plt.Line(x=None,
                                            y=[1 / Input["quant_up"]],
                                            label=f'upper quantile, f={Input["quant_up"]} Hz',
                                            color='green',
                                            linestyle=':')
                tile_curr.add_line(Line_quant_up)

                Line_quant_low = hc_plt.Line(x=None,
                                             y=[1 / Input["quant_low"]],
                                             label=f'lower quantile, f={Input["quant_low"]} Hz',
                                             color='green',
                                             linestyle='--')
                tile_curr.add_line(Line_quant_low)

            scatter = hc_plt.Scatter(x=point_data[Seg.colnames["x"]],
                                     y=point_data[Seg.colnames["y"]],
                                     cmap='cool',
                                     size=2,
                                     cmap_norm='sqrt')

            tile_curr.add_scatter(scatter)

            if Seg.angles is not None:
                Tiles.append(tile_curr)

            else:
                Tiles_omni.append(tile_curr)

        FIG_direc = hc_plt.plot_tiled(Tiles, global_max=['auto', 'auto'], global_min=[0, 0], grid=[3, 2], figsize=figsize_fullpage)

        FIG_omni = hc_plt.plot_tiled(Tiles_omni, global_max=['auto', 'auto'], global_min=[0, 0], grid=[1, 1], figsize=figsize_halfpage)

        if 'png' in INPUT["Toggle_Modules"]["plot_as"]:
            gl.save_figs_as_png(FIG_direc + FIG_omni, path_out + 'HSTP_swell', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

        if 'pdf' in INPUT["Toggle_Modules"]["plot_as"]:
            gl.save_figs_as_pdf(FIG_direc + FIG_omni, path_out + 'HSTP_swell', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

if 'total' in INPUT["Toggle_Modules"].get("plot_HSTP", {}):
    print('plotting HSTP total...')

    if not ('total' in INPUT["Toggle_Modules"].get("calc_HSTP", {})):
        print("   please toggle calculation to plot")

    else:

        Calc = DATA_OUT["HSTP"]["total"]
        Input = INPUT["HSTP_total"]
        Tiles = []
        Tiles_omni = []

        df = Calc.load_from_db(colnames_ini=True)
        titels = Calc.create_segment_title()
        titels = gl.alias(titels, COLNAMES, INPUT["Aliase"])

        for i, Seg in enumerate(Calc.result):

            Seg.indizes = pd.to_datetime(Seg.indizes)
            point_data = df[df.index.isin(Seg.indizes)]

            tile_curr = hc_plt.Tile(i,
                                    x_label=gl.alias(Seg.colnames['x'], COLNAMES, INPUT["Aliase"]),
                                    y_label=gl.alias(Seg.colnames['y'], COLNAMES, INPUT["Aliase"]),
                                    title=titels[i])

            key_plot = [key for key in Seg.result.columns if 'plot' in key]
            key_mean = [key for key in key_plot if 'mean' in key]
            key_perc = [key for key in key_plot if 'percentile' in key]
            key_quantile = [key for key in Seg.result.columns if 'quantile' in key]

            Line_mean = hc_plt.Line(x=Seg.result["x"],
                                    y=Seg.result["mean result plot"],
                                    label='mean',
                                    color='black',
                                    linestyle='-')
            tile_curr.add_line(Line_mean)

            # Line_perc_low = hc_plt.Line(x=Seg.result["x"],
            #                             y=Seg.result[key_perc[0]],
            #                             label=key_perc[0].replace('result', '').replace('plot', ''),
            #                             color='black',
            #                             linestyle='--')
            # tile_curr.add_line(Line_perc_low)

            # Line_perc_up = hc_plt.Line(x=Seg.result["x"],
            #                            y=Seg.result[key_perc[1]],
            #                            label='mean',
            #                            color='black',
            #                            linestyle=':')
            # tile_curr.add_line(Line_perc_up)

            if len(key_quantile) > 0:
                Line_quant = hc_plt.Line(x=Seg.result["x"],
                                         y=Seg.result[key_quantile[0]],
                                         label='selected correlation',
                                         color='red',
                                         linestyle='-')
                tile_curr.add_line(Line_quant)

                Line_quant_up = hc_plt.Line(x=None,
                                            y=[1 / Input["quant_up"]],
                                            label=f'upper quantile, f={Input["quant_up"]} Hz',
                                            color='green',
                                            linestyle=':')
                tile_curr.add_line(Line_quant_up)

                Line_quant_low = hc_plt.Line(x=None,
                                             y=[1 / Input["quant_low"]],
                                             label=f'lower quantile, f={Input["quant_low"]} Hz',
                                             color='green',
                                             linestyle='--')
                tile_curr.add_line(Line_quant_low)

            scatter = hc_plt.Scatter(x=point_data[Seg.colnames["x"]],
                                     y=point_data[Seg.colnames["y"]],
                                     cmap='cool',
                                     size=2,
                                     cmap_norm='sqrt')

            tile_curr.add_scatter(scatter)

            if Seg.angles is not None:
                Tiles.append(tile_curr)

            else:
                Tiles_omni.append(tile_curr)

        FIG_direc = hc_plt.plot_tiled(Tiles, global_max=['auto', 'auto'], global_min=[0, 0], grid=[3, 2], figsize=figsize_fullpage)

        FIG_omni = hc_plt.plot_tiled(Tiles_omni, global_max=['auto', 'auto'], global_min=[0, 0], grid=[1, 1], figsize=figsize_halfpage)

        if 'png' in INPUT["Toggle_Modules"]["plot_as"]:
            gl.save_figs_as_png(FIG_direc + FIG_omni, path_out + 'HSTP_total', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

        if 'pdf' in INPUT["Toggle_Modules"]["plot_as"]:
            gl.save_figs_as_pdf(FIG_direc + FIG_omni, path_out + 'HSTP_total', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

if 'wind' in INPUT["Toggle_Modules"].get("plot_VMTP", {}):
    print('plotting VMTP wind...')

    if not ('wind' in INPUT["Toggle_Modules"].get("calc_VMTP", {})):
        print("   please toggle calculation to plot")
    else:
        Calc = DATA_OUT["VMTP"]["wind"]

        Tiles = []
        Tiles_omni = []

        df = Calc.load_from_db(colnames_ini=True)

        titels = Calc.create_segment_title()
        titels = gl.alias(titels, COLNAMES, INPUT["Aliase"])

        for i, Seg in enumerate(Calc.result):

            Seg.indizes = pd.to_datetime(Seg.indizes)
            point_data = df[df.index.isin(Seg.indizes)]

            tile_curr = hc_plt.Tile(i,
                                    x_label=gl.alias(Seg.colnames['x'], COLNAMES, INPUT["Aliase"]),
                                    y_label=gl.alias(Seg.colnames['y'], COLNAMES, INPUT["Aliase"]),
                                    title=titels[i])

            Line_mean = hc_plt.Line(x=Seg.result["x"],
                                    y=Seg.result["mean result plot"],
                                    label='mean',
                                    color='black')

            tile_curr.add_line(Line_mean)

            scatter = hc_plt.Scatter(x=point_data[Seg.colnames["x"]],
                                     y=point_data[Seg.colnames["y"]],
                                     cmap='cool',
                                     size=2,
                                     cmap_norm='sqrt')

            tile_curr.add_scatter(scatter)

            if Seg.angles is not None:
                Tiles.append(tile_curr)
            else:
                Tiles_omni.append(tile_curr)

        FIG_direc = hc_plt.plot_tiled(Tiles, global_max=['auto', 'auto'], global_min=[0, 0], grid=[3, 2], figsize=figsize_fullpage)

        FIG_omni = hc_plt.plot_tiled(Tiles_omni, global_max=['auto', 'auto'], global_min=[0, 0], grid=[1, 1], figsize=figsize_halfpage)

        if 'png' in INPUT["Toggle_Modules"]["plot_as"]:
            gl.save_figs_as_png(FIG_direc + FIG_omni, path_out + 'VMTP_wind', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

        if 'pdf' in INPUT["Toggle_Modules"]["plot_as"]:
            gl.save_figs_as_pdf(FIG_direc + FIG_omni, path_out + 'VMTP_wind', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

if 'swell' in INPUT["Toggle_Modules"].get("plot_VMTP", {}):

    print('plotting VMTP swell...')

    if not ('swell' in INPUT["Toggle_Modules"].get("calc_VMTP", {})):
        print("   please toggle calculation to plot")
    else:

        Calc = DATA_OUT["VMTP"]["swell"]

        Tiles = []
        Tiles_omni = []

        df = Calc.load_from_db(colnames_ini=True)

        titels = Calc.create_segment_title()
        titels = gl.alias(titels, COLNAMES, INPUT["Aliase"])

        for i, Seg in enumerate(Calc.result):

            Seg.indizes = pd.to_datetime(Seg.indizes)
            point_data = df[df.index.isin(Seg.indizes)]

            tile_curr = hc_plt.Tile(i,
                                    x_label=gl.alias(Seg.colnames['x'], COLNAMES, INPUT["Aliase"]),
                                    y_label=gl.alias(Seg.colnames['y'], COLNAMES, INPUT["Aliase"]),
                                    title=titels[i])

            Line_mean = hc_plt.Line(x=Seg.result["x"],
                                    y=Seg.result["mean result plot"],
                                    label='mean',
                                    color='black')

            tile_curr.add_line(Line_mean)

            scatter = hc_plt.Scatter(x=point_data[Seg.colnames["x"]],
                                     y=point_data[Seg.colnames["y"]],
                                     cmap='cool',
                                     size=2,
                                     cmap_norm='sqrt')

            tile_curr.add_scatter(scatter)

            if Seg.angles is not None:
                Tiles.append(tile_curr)
            else:
                Tiles_omni.append(tile_curr)

        FIG_direc = hc_plt.plot_tiled(Tiles, global_max=['auto', 'auto'], global_min=[0, 0], grid=[3, 2], figsize=figsize_fullpage)

        FIG_omni = hc_plt.plot_tiled(Tiles_omni, global_max=['auto', 'auto'], global_min=[0, 0], grid=[1, 1], figsize=figsize_halfpage)

        if 'png' in INPUT["Toggle_Modules"]["plot_as"]:
            gl.save_figs_as_png(FIG_direc + FIG_omni, path_out + 'VMTP_swell', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

        if 'pdf' in INPUT["Toggle_Modules"]["plot_as"]:
            gl.save_figs_as_pdf(FIG_direc + FIG_omni, path_out + 'VMTP_swell', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

if 'total' in INPUT["Toggle_Modules"].get("plot_VMTP", {}):

    print('plotting VMTP total...')

    if not ('total' in INPUT["Toggle_Modules"].get("calc_VMTP", {})):
        print("   please toggle calculation to plot")
    else:

        Calc = DATA_OUT["VMTP"]["total"]

        Tiles = []
        Tiles_omni = []

        df = Calc.load_from_db(colnames_ini=True)

        titels = Calc.create_segment_title()
        titels = gl.alias(titels, COLNAMES, INPUT["Aliase"])

        for i, Seg in enumerate(Calc.result):

            Seg.indizes = pd.to_datetime(Seg.indizes)
            point_data = df[df.index.isin(Seg.indizes)]

            tile_curr = hc_plt.Tile(i,
                                    x_label=gl.alias(Seg.colnames['x'], COLNAMES, INPUT["Aliase"]),
                                    y_label=gl.alias(Seg.colnames['y'], COLNAMES, INPUT["Aliase"]),
                                    title=titels[i])

            Line_mean = hc_plt.Line(x=Seg.result["x"],
                                    y=Seg.result["mean result plot"],
                                    label='mean',
                                    color='black')

            tile_curr.add_line(Line_mean)

            scatter = hc_plt.Scatter(x=point_data[Seg.colnames["x"]],
                                     y=point_data[Seg.colnames["y"]],
                                     cmap='cool',
                                     size=2,
                                     cmap_norm='sqrt')

            tile_curr.add_scatter(scatter)

            if Seg.angles is not None:
                Tiles.append(tile_curr)
            else:
                Tiles_omni.append(tile_curr)

        FIG_direc = hc_plt.plot_tiled(Tiles, global_max=['auto', 'auto'], global_min=[0, 0], grid=[3, 2], figsize=figsize_fullpage)

        FIG_omni = hc_plt.plot_tiled(Tiles_omni, global_max=['auto', 'auto'], global_min=[0, 0], grid=[1, 1], figsize=figsize_halfpage)

        if 'png' in INPUT["Toggle_Modules"]["plot_as"]:
            gl.save_figs_as_png(FIG_direc + FIG_omni, path_out + 'VMTP_total', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

        if 'pdf' in INPUT["Toggle_Modules"]["plot_as"]:
            gl.save_figs_as_pdf(FIG_direc + FIG_omni, path_out + 'VMTP_total', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

if 'wind' in INPUT["Toggle_Modules"].get("plot_RWI", {}):

    print('plotting RWI wind...')

    if not ('wind' in INPUT["Toggle_Modules"].get("calc_RWI", {})):
        print("   please toggle calculation to plot")
    else:

        Calc = DATA_OUT["RWI"]["wind"]

        Tiles = []
        Tiles_omni = []

        df = Calc.load_from_db(colnames_ini=True)

        titels = Calc.create_segment_title()
        titels = gl.alias(titels, COLNAMES, INPUT["Aliase"])

        for i, Seg in enumerate(Calc.result):

            Seg.indizes = pd.to_datetime(Seg.indizes)
            point_data = df[df.index.isin(Seg.indizes)]

            tile_curr = hc_plt.Tile(i,
                                    x_label=gl.alias(Seg.colnames['x'], COLNAMES, INPUT["Aliase"]),
                                    y_label=gl.alias(Seg.colnames['y'], COLNAMES, INPUT["Aliase"]),
                                    title=titels[i])

            scatter = hc_plt.Scatter(x=point_data[Seg.colnames["x"]],
                                     y=point_data[Seg.colnames["y"]],
                                     cmap='cool',
                                     cmap_mode='manual',
                                     c=Seg.result.values,
                                     size=2,
                                     cbar=True,
                                     cbar_label="RWI = $\\sqrt{S(f_0)}$ (Resonance Wave Intesity) $[\\sqrt{m^2/Hz}]$",
                                     cbar_label_fontsize=6)

            tile_curr.add_scatter(scatter)

            Line_f0 = hc_plt.Line(x=None,
                                  y=[1 / INPUT["Structure"]["f_0"]],
                                  label=f'f_0, f={INPUT["Structure"]["f_0"]} Hz',
                                  color='green',
                                  linestyle=':')
            tile_curr.add_line(Line_f0)

            if Seg.angles is not None:
                Tiles.append(tile_curr)
            else:
                Tiles_omni.append(tile_curr)

        FIG_direc = hc_plt.plot_tiled(Tiles, global_max=['auto', 'auto'], global_min=[0, 0], grid=[3, 2], scatter_max='auto', figsize=figsize_fullpage)

        FIG_omni = hc_plt.plot_tiled(Tiles_omni, global_max=['auto', 'auto'], global_min=[0, 0], grid=[1, 1], figsize=figsize_halfpage)

        if 'png' in INPUT["Toggle_Modules"]["plot_as"]:
            gl.save_figs_as_png(FIG_direc + FIG_omni, path_out + 'RWI_wind', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

        if 'pdf' in INPUT["Toggle_Modules"]["plot_as"]:
            gl.save_figs_as_pdf(FIG_direc + FIG_omni, path_out + 'RWI_wind', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

if 'total' in INPUT["Toggle_Modules"].get("plot_RWI", {}):

    print('plotting RWI total...')

    if not ('total' in INPUT["Toggle_Modules"].get("calc_RWI", {})):
        print("   please toggle calculation to plot")
    else:

        Calc = DATA_OUT["RWI"]["total"]

        Tiles = []
        Tiles_omni = []

        df = Calc.load_from_db(colnames_ini=True)

        titels = Calc.create_segment_title()
        titels = gl.alias(titels, COLNAMES, INPUT["Aliase"])

        for i, Seg in enumerate(Calc.result):

            Seg.indizes = pd.to_datetime(Seg.indizes)
            point_data = df[df.index.isin(Seg.indizes)]

            tile_curr = hc_plt.Tile(i,
                                    x_label=gl.alias(Seg.colnames['x'], COLNAMES, INPUT["Aliase"]),
                                    y_label=gl.alias(Seg.colnames['y'], COLNAMES, INPUT["Aliase"]),
                                    title=titels[i])

            scatter = hc_plt.Scatter(x=point_data[Seg.colnames["x"]],
                                     y=point_data[Seg.colnames["y"]],
                                     cmap='cool',
                                     cmap_mode='manual',
                                     c=Seg.result.values,
                                     size=2,
                                     cbar=True,
                                     cbar_label="RWI = $\\sqrt{S(f_0)}$ (Resonance Wave Intesity) $[\\sqrt{m^2/Hz}]$",
                                     cbar_label_fontsize=6)

            tile_curr.add_scatter(scatter)

            Line_f0 = hc_plt.Line(x=None,
                                  y=[1 / INPUT["Structure"]["f_0"]],
                                  label=f'f_0, f={INPUT["Structure"]["f_0"]} Hz',
                                  color='green',
                                  linestyle=':')
            tile_curr.add_line(Line_f0)

            if Seg.angles is not None:
                Tiles.append(tile_curr)
            else:
                Tiles_omni.append(tile_curr)

        FIG_direc = hc_plt.plot_tiled(Tiles, global_max=['auto', 'auto'], global_min=[0, 0], grid=[3, 2], scatter_max='auto', figsize=figsize_fullpage)

        FIG_omni = hc_plt.plot_tiled(Tiles_omni, global_max=['auto', 'auto'], global_min=[0, 0], grid=[1, 1], figsize=figsize_halfpage)

        if 'png' in INPUT["Toggle_Modules"]["plot_as"]:
            gl.save_figs_as_png(FIG_direc + FIG_omni, path_out + 'RWI_total', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

        if 'pdf' in INPUT["Toggle_Modules"]["plot_as"]:
            gl.save_figs_as_pdf(FIG_direc + FIG_omni, path_out + 'RWI_total', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

if 'wind' in INPUT["Toggle_Modules"].get("plot_BreakSteep", {}):

    print('plotting WaveBreak_Steep wind...')

    if not ('wind' in INPUT["Toggle_Modules"].get("calc_WaveBreak_Steep", {})):
        print("   please toggle calculation to plot")
    else:

        Calc = DATA_OUT["WaveBreak_Steep"]["wind"]

        Tiles = []
        Tiles_omni = []

        df = Calc.load_from_db(colnames_ini=True)

        titels = Calc.create_segment_title()
        titels = gl.alias(titels, COLNAMES, INPUT["Aliase"])

        for i, Seg in enumerate(Calc.result):

            Seg.indizes = pd.to_datetime(Seg.indizes)
            point_data = df[df.index.isin(Seg.indizes)]

            tile_curr = hc_plt.Tile(i,
                                    x_label=gl.alias(Seg.colnames['x'], COLNAMES, INPUT["Aliase"]),
                                    y_label=gl.alias(Seg.colnames['y'], COLNAMES, INPUT["Aliase"]),
                                    title=titels[i])

            c_krit = Seg.result["steepness"]

            c_krit.iloc[Seg.result["bool_break"] == False] = float('nan')

            scatter = hc_plt.Scatter(x=point_data[Seg.colnames["x"]],
                                     y=point_data[Seg.colnames["y"]],
                                     cmap='cool',
                                     cmap_mode='manual',
                                     c=c_krit.values,
                                     size=2,
                                     cbar=True,
                                     cbar_label="steepness [-]")

            tile_curr.add_scatter(scatter)

            if Seg.angles is not None:
                Tiles.append(tile_curr)
            else:
                Tiles_omni.append(tile_curr)

        FIG_direc = hc_plt.plot_tiled(Tiles, global_max=['auto', 'auto'], global_min=[0, 0], grid=[3, 2], scatter_max='auto', scatter_min=INPUT["Structure"]["steep_crit"],
                                      figsize=figsize_fullpage)

        FIG_omni = hc_plt.plot_tiled(Tiles_omni, global_max=['auto', 'auto'], global_min=[0, 0], grid=[1, 1], figsize=figsize_halfpage)

        if 'png' in INPUT["Toggle_Modules"]["plot_as"]:
            gl.save_figs_as_png(FIG_direc + FIG_omni, path_out + 'WaveBreak_wind', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

        if 'pdf' in INPUT["Toggle_Modules"]["plot_as"]:
            gl.save_figs_as_pdf(FIG_direc + FIG_omni, path_out + 'WaveBreak_wind', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

if 'swell' in INPUT["Toggle_Modules"].get("plot_BreakSteep", {}):

    print('plotting WaveBreak_Steep swell...')

    if not ('swell' in INPUT["Toggle_Modules"].get("calc_WaveBreak_Steep", {})):
        print("   please toggle calculation to plot")
    else:

        Calc = DATA_OUT["WaveBreak_Steep"]["swell"]

        Tiles = []
        Tiles_omni = []

        df = Calc.load_from_db(colnames_ini=True)

        titels = Calc.create_segment_title()
        titels = gl.alias(titels, COLNAMES, INPUT["Aliase"])

        for i, Seg in enumerate(Calc.result):

            Seg.indizes = pd.to_datetime(Seg.indizes)
            point_data = df[df.index.isin(Seg.indizes)]

            tile_curr = hc_plt.Tile(i,
                                    x_label=gl.alias(Seg.colnames['x'], COLNAMES, INPUT["Aliase"]),
                                    y_label=gl.alias(Seg.colnames['y'], COLNAMES, INPUT["Aliase"]),
                                    title=titels[i])

            c_krit = Seg.result["steepness"]

            c_krit.iloc[Seg.result["bool_break"] == False] = float('nan')

            scatter = hc_plt.Scatter(x=point_data[Seg.colnames["x"]],
                                     y=point_data[Seg.colnames["y"]],
                                     cmap='cool',
                                     cmap_mode='manual',
                                     c=c_krit.values,
                                     size=2,
                                     cbar=True,
                                     cbar_label="steepness [-]")

            tile_curr.add_scatter(scatter)

            if Seg.angles is not None:
                Tiles.append(tile_curr)
            else:
                Tiles_omni.append(tile_curr)

        FIG_direc = hc_plt.plot_tiled(Tiles, global_max=['auto', 'auto'], global_min=[0, 0], grid=[3, 2], scatter_max='auto', scatter_min=INPUT["Structure"]["steep_crit"],
                                      figsize=figsize_fullpage)

        FIG_omni = hc_plt.plot_tiled(Tiles_omni, global_max=['auto', 'auto'], global_min=[0, 0], grid=[1, 1], figsize=figsize_halfpage)

        if 'png' in INPUT["Toggle_Modules"]["plot_as"]:
            gl.save_figs_as_png(FIG_direc + FIG_omni, path_out + 'WaveBreak_swell', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

        if 'pdf' in INPUT["Toggle_Modules"]["plot_as"]:
            gl.save_figs_as_pdf(FIG_direc + FIG_omni, path_out + 'WaveBreak_swell', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

if 'wind' in INPUT["Toggle_Modules"].get("plot_Tables", {}):
    print('plotting Tables wind...')

    if not ('wind' in INPUT["Toggle_Modules"].get("calc_Tables", {})):
        print("   please toggle calculation to plot")
    else:
        # VMHS
        Calc = DATA_OUT["table_vmhs"]["wind"]

        titel = f"'{gl.alias(Calc.result[0].colnames['y'], COLNAMES, INPUT['Aliase'])}'" + "\n " + f"in '{gl.alias(Calc.result[1].angle_name, COLNAMES, INPUT['Aliase'])}' directional sections" + "\n" + r"\small " + f"with v_m = '{gl.alias(Calc.result[0].colnames['x'], COLNAMES, INPUT['Aliase'])}'"

        FIG = [hc_plt.plot_table_condesation(Calc, figsize=figsize_fullpage, titel=titel)]

        if 'png' in INPUT["Toggle_Modules"]["plot_as"]:
            gl.save_figs_as_png(FIG, path_out + 'table_vmhs_wind', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

        if 'pdf' in INPUT["Toggle_Modules"]["plot_as"]:
            gl.save_figs_as_pdf(FIG, path_out + 'table_vmhs_wind', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

        # VMTP
        Calc = DATA_OUT["table_vmtp"]["wind"]
        titel = f"'{gl.alias(Calc.result[0].colnames['y'], COLNAMES, INPUT['Aliase'])}'" + "\n " + f"in '{gl.alias(Calc.result[1].angle_name, COLNAMES, INPUT['Aliase'])}' directional sections" + "\n" + r"\small " + f"with v_m = '{gl.alias(Calc.result[0].colnames['x'], COLNAMES, INPUT['Aliase'])}'"

        FIG = [hc_plt.plot_table_condesation(Calc, figsize=figsize_fullpage, titel=titel)]

        if 'png' in INPUT["Toggle_Modules"]["plot_as"]:
            gl.save_figs_as_png(FIG, path_out + 'table_vmtp_wind', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

        if 'pdf' in INPUT["Toggle_Modules"]["plot_as"]:
            gl.save_figs_as_pdf(FIG, path_out + 'table_vmtp_wind', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

if 'swell' in INPUT["Toggle_Modules"].get("plot_Tables", {}):
    print('plotting Tables swell...')

    if not ('swell' in INPUT["Toggle_Modules"].get("calc_Tables", {})):
        print("   please toggle calculation to plot")
    else:
        # VMHS
        Calc = DATA_OUT["table_vmhs"]["swell"]
        titel = f"'{gl.alias(Calc.result[0].colnames['y'], COLNAMES, INPUT['Aliase'])}'" + "\n " + f"in '{gl.alias(Calc.result[1].angle_name, COLNAMES, INPUT['Aliase'])}' directional sections" + "\n" + r"\small " + f"with v_m = '{gl.alias(Calc.result[0].colnames['x'], COLNAMES, INPUT['Aliase'])}'"
        FIG = [hc_plt.plot_table_condesation(Calc, figsize=figsize_fullpage, titel=titel)]

        if 'png' in INPUT["Toggle_Modules"]["plot_as"]:
            gl.save_figs_as_png(FIG, path_out + 'table_vmhs_swell', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

        if 'pdf' in INPUT["Toggle_Modules"]["plot_as"]:
            gl.save_figs_as_pdf(FIG, path_out + 'table_vmhs_swell', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

        # VMTP
        Calc = DATA_OUT["table_vmtp"]["swell"]
        titel = f"'{gl.alias(Calc.result[0].colnames['y'], COLNAMES, INPUT['Aliase'])}'" + "\n " + f"in '{gl.alias(Calc.result[1].angle_name, COLNAMES, INPUT['Aliase'])}' directional sections" + "\n" + r"\small " + f"with v_m = '{gl.alias(Calc.result[0].colnames['x'], COLNAMES, INPUT['Aliase'])}'"

        FIG = [hc_plt.plot_table_condesation(Calc, figsize=figsize_fullpage, titel=titel)]

        if 'png' in INPUT["Toggle_Modules"]["plot_as"]:
            gl.save_figs_as_png(FIG, path_out + 'table_vmtp_swell', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

        if 'pdf' in INPUT["Toggle_Modules"]["plot_as"]:
            gl.save_figs_as_pdf(FIG, path_out + 'table_vmtp_swell', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

if 'total' in INPUT["Toggle_Modules"].get("plot_Tables", {}):
    print('plotting Tables total...')

    if not ('total' in INPUT["Toggle_Modules"].get("calc_Tables", {})):
        print("   please toggle calculation to plot")
    else:
        # VMHS
        Calc = DATA_OUT["table_vmhs"]["total"]
        titel = f"'{gl.alias(Calc.result[0].colnames['y'], COLNAMES, INPUT['Aliase'])}'" + "\n " + f"in '{gl.alias(Calc.result[1].angle_name, COLNAMES, INPUT['Aliase'])}' directional sections" + "\n" + r"\small " + f"with v_m = '{gl.alias(Calc.result[0].colnames['x'], COLNAMES, INPUT['Aliase'])}'"
        FIG = [hc_plt.plot_table_condesation(Calc, figsize=figsize_fullpage, titel=titel)]

        if 'png' in INPUT["Toggle_Modules"]["plot_as"]:
            gl.save_figs_as_png(FIG, path_out + 'table_vmhs_total', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

        if 'pdf' in INPUT["Toggle_Modules"]["plot_as"]:
            gl.save_figs_as_pdf(FIG, path_out + 'table_vmhs_total', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

        # VMTP
        Calc = DATA_OUT["table_vmtp"]["total"]
        titel = f"'{gl.alias(Calc.result[0].colnames['y'], COLNAMES, INPUT['Aliase'])}'" + "\n " + f"in '{gl.alias(Calc.result[1].angle_name, COLNAMES, INPUT['Aliase'])}' directional sections" + "\n" + r"\small " + f"with v_m = '{gl.alias(Calc.result[0].colnames['x'], COLNAMES, INPUT['Aliase'])}'"

        FIG = [hc_plt.plot_table_condesation(Calc, figsize=figsize_fullpage, titel=titel)]

        if 'png' in INPUT["Toggle_Modules"]["plot_as"]:
            gl.save_figs_as_png(FIG, path_out + 'table_vmtp_total', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

        if 'pdf' in INPUT["Toggle_Modules"]["plot_as"]:
            gl.save_figs_as_pdf(FIG, path_out + 'table_vmtp_total', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

if INPUT["Toggle_Modules"].get("plot_AngleDeviation", {}):
    print('plotting AngleDeviation...')

    if not (len(INPUT["Toggle_Modules"].get("calc_AngleDeviation", {})) > 2):
        print("   please toggle calculation to plot")
    else:
        # scatter_plot omni
        Calc = DATA_OUT["AngleDeviation"]

        FIG_Tables = []
        FIG_scatter = []
        df = Calc.load_from_db(colnames_ini=True)
        # title management
        title = r"\small Occurrence probability of missalinment " + "\n" + f" ({Calc.result[0].colnames['ang_comp']} - {Calc.result[0].colnames['ang_orig']})"
        subtitle = Calc.create_segment_title(mode='sparse')
        subsubtitle = f"with v_m = {Calc.result[1].colnames['v_m']}"
        titles = [title + "\n " + subtitle_curr + "\n " + subsubtitle for subtitle_curr in subtitle]
        titles = gl.alias(titles, COLNAMES, INPUT["Aliase"])

        for i, Seg in enumerate(Calc.result):

            if Seg.angles is not None:
                temp = hc_plt.table(data=Seg.result.values,
                                    collabels=[str(deg_curr) + "°" for deg_curr in list(Seg.result.columns)],
                                    rowlabels=list(Seg.result.index),
                                    row_label_name="v_m",
                                    titel=titles[i],
                                    formater='.2e',
                                    heatmap=True,
                                    figsize=figsize_fullpage)

                FIG_Tables.append(temp)
            else:
                x_label = gl.alias(Seg.colnames['ang_orig'], COLNAMES, INPUT["Aliase"])
                tile_scatter = hc_plt.Tile(i, x_label=x_label, y_label='deviation [°]',
                                           title=f"Missalignment to {Calc.result[0].colnames['ang_comp']}" + "\n" + f"over {Calc.result[1].colnames['ang_comp']}")

                point_data = df[df.index.isin(Seg.indizes)]

                scatter = hc_plt.Scatter(x=point_data[Seg.colnames['ang_orig']],
                                         y=Seg.result["points"]["diff"],
                                         cmap='cool',
                                         size=2
                                         )

                tile_scatter.add_scatter(scatter)

                line = hc_plt.Line(x=Seg.result["mean"].index,
                                   y=Seg.result["mean"].values,
                                   label=f"rolling absolute mean (global mean = {round(np.mean(Seg.result['mean'].values), 2)} deg)",
                                   color='black')

                tile_scatter.add_line(line)

                temp = hc_plt.plot_tiled([tile_scatter], grid=[1, 1], figsize=figsize_halfpage)
                FIG_scatter = FIG_scatter + temp

        if 'png' in INPUT["Toggle_Modules"]["plot_as"]:
            gl.save_figs_as_png(FIG_scatter, path_out + 'angle_deviation_scatter', dpi=INPUT["Toggle_Modules"]["dpi_figures"])
            gl.save_figs_as_png(FIG_Tables, path_out + 'angle_deviation_table', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

        if 'pdf' in INPUT["Toggle_Modules"]["plot_as"]:
            gl.save_figs_as_pdf(FIG_scatter, path_out + 'angle_deviation_scatter', dpi=INPUT["Toggle_Modules"]["dpi_figures"])
            gl.save_figs_as_pdf(FIG_Tables, path_out + 'angle_deviation_table', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

if INPUT["Toggle_Modules"].get("plot_Roseplots", {}) and len(INPUT["Toggle_Modules"].get("calc_Roseplots", {})) > 0:
    print("plotting Roseplots...")
    i = 0
    Tiles = []
    for Roseplot_name, Calc in DATA_OUT["Roseplot"].items():
        titel = f'Roseplot {Calc.basedata["colnames_ini"][1]} over' + "\n" + f'{Calc.basedata["colnames_ini"][0]}'

        titel = gl.alias(titel, COLNAMES, INPUT["Aliase"])
        radial = Calc.result["table"].div(Calc.basedata['N_rows'] / 100)
        radial = [radial[col].tolist() for col in radial]

        Tile = hc_plt.PolarTile(i, title=titel)

        temp = hc_plt.RoseBar(angles=Calc.result["table"].columns.values,
                              r_bins=Calc.result["r_bins"],
                              radial_data=radial,
                              radial_mode='summed',
                              radial_datatype='percent',
                              cbar=True,
                              cbar_label=gl.alias(Calc.basedata['colnames_ini'][1], COLNAMES, INPUT["Aliase"]),
                              r_max=Calc.result["r_max"]
                              )

        Tile.add_RoseBar(temp)

        Tiles.append(Tile)

        i = i + 1

    FIG = hc_plt.plot_tiled(Tiles, figsize=figsize_fullpage)

    if 'png' in INPUT["Toggle_Modules"]["plot_as"]:
        gl.save_figs_as_png(FIG, path_out + 'Roseplots', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

    if 'pdf' in INPUT["Toggle_Modules"]["plot_as"]:
        gl.save_figs_as_pdf(FIG, path_out + 'Roseplots', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

if INPUT["Toggle_Modules"].get("plot_ExtremeValues", {}) and len(INPUT["Toggle_Modules"].get("calc_ExtremeValues", {})) > 0:
    print("plotting ExtremeValues...")

    for Calc_name, Calc in DATA_OUT["ExtremeValues"].items():

        # Timeseries
        i = 0
        Tiles = []
        Tiles_omni = []
        df = Calc.load_from_db(colnames_ini=True)
        titels = Calc.create_segment_title()
        titels = gl.alias(titels, COLNAMES, INPUT["Aliase"])

        for i, Seg in enumerate(Calc.result):

            Seg.indizes = pd.to_datetime(Seg.indizes)
            df_seg = df[df.index.isin(Seg.indizes)]

            tile_curr = hc_plt.Tile(i,
                                    x_label='date',
                                    y_label=gl.alias(Seg.colnames['x'], COLNAMES, INPUT["Aliase"]),
                                    title=titels[i])

            Line = hc_plt.Line(x=df_seg[Seg.colnames['x']].index,
                               y=df_seg[Seg.colnames['x']].values,
                               color='black',
                               linewidth=1)

            tile_curr.add_line(Line)

            scatter = hc_plt.Scatter(x=Seg.result["points"]["x_max"].index,
                                     y=Seg.result["points"]["x_max"].values,
                                     color='red',
                                     label='Extreme Values',
                                     size=1)

            tile_curr.add_scatter(scatter)

            if Seg.angles is not None:
                Tiles.append(tile_curr)
            else:
                Tiles_omni.append(tile_curr)

        FIG_direc = hc_plt.plot_tiled(Tiles, grid=[3, 2], figsize=figsize_fullpage)

        FIG_omni = hc_plt.plot_tiled(Tiles_omni, grid=[1, 1], figsize=figsize_halfpage)

        if 'png' in INPUT["Toggle_Modules"]["plot_as"]:
            gl.save_figs_as_png(FIG_direc + FIG_omni, path_out + f'Timeseries_{Calc_name}', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

        if 'pdf' in INPUT["Toggle_Modules"]["plot_as"]:
            gl.save_figs_as_pdf(FIG_direc + FIG_omni, path_out + f'Timeseries_{Calc_name}', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

        # qq
        i = 0
        Tiles = []
        Tiles_Omni = []

        for i, Seg in enumerate(Calc.result):
            title = titels[i] + "\n" + r"\scriptsize " + (f'intervall mode: {Seg.result["meta"]["intervall_mode"]}, '
                                                          f'intervall algorithm: {Seg.result["meta"]["intervall_algorithm"]}, '
                                                          f'itterations: {Seg.result["meta"]["N_itter"]}')
            tile_curr = hc_plt.Tile(i,
                                    x_label='Anual Maxima of' + '\n' + gl.alias(Seg.colnames['x'], COLNAMES, INPUT["Aliase"]),
                                    y_label='Theoretical Maxima (gumbel)',
                                    title=title)

            scatter = hc_plt.Scatter(x=Seg.result["points"]["x_max"].values,
                                     y=Seg.result["points"]["x_theorie"].values,
                                     color='red',
                                     size=2,
                                     label='Maximal Values')

            tile_curr.add_scatter(scatter)

            Line = hc_plt.Line(x=[min(Seg.result["points"]["x_theorie"].values), max(Seg.result["points"]["x_theorie"].values)],
                               y=[min(Seg.result["points"]["x_theorie"].values), max(Seg.result["points"]["x_theorie"].values)],
                               color='black')

            tile_curr.add_line(Line)

            if INPUT["ExtremeValues"]["conf_inter_mode"] == "percentile":
                label_up = f'upper bound ({INPUT["ExtremeValues"]["perc_up"]}th percentile)'
                label_down = f'lower bound ({INPUT["ExtremeValues"]["perc_down"]}th percentile)'

            if INPUT["ExtremeValues"]["conf_inter_mode"] == "std":
                label_up = f'upper bound (mean + std_deviation)'
                label_down = f'lower bound (mean - std_deviation)'

            Line = hc_plt.Line(x=Seg.result["points"]["x_theorie"].values,
                               y=Seg.result["points"]["band_up"].values,
                               color='black',
                               linestyle='--',
                               label=label_up)

            tile_curr.add_line(Line)

            Line = hc_plt.Line(x=Seg.result["points"]["x_theorie"].values,
                               y=Seg.result["points"]["band_down"].values,
                               color='black',
                               linestyle=':',
                               label=label_down)

            tile_curr.add_line(Line)

            if Seg.angles is not None:
                Tiles.append(tile_curr)
            else:
                Tiles_omni.append(tile_curr)

        FIG_direc = hc_plt.plot_tiled(Tiles, global_max=[None, None], global_min=[None, None], grid=[3, 2], figsize=figsize_fullpage)

        FIG_omni = hc_plt.plot_tiled(Tiles_omni, global_max=[None, None], global_min=[None, None], grid=[1, 1], figsize=figsize_halfpage)

        if 'png' in INPUT["Toggle_Modules"]["plot_as"]:
            gl.save_figs_as_png(FIG_direc + FIG_omni, path_out + f'qq_{Calc_name}', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

        if 'pdf' in INPUT["Toggle_Modules"]["plot_as"]:
            gl.save_figs_as_pdf(FIG_direc + FIG_omni, path_out + f'qq_{Calc_name}', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

        # T_Return
        i = 0
        Tiles = []
        Tiles_Omni = []
        titels = Calc.create_segment_title()
        titels = gl.alias(titels, COLNAMES, INPUT["Aliase"])

        for i, Seg in enumerate(Calc.result):

            title = titels[i] + "\n" + r"\scriptsize " + (f'intervall mode: {Seg.result["meta"]["intervall_mode"]}, '
                                                          f'intervall algorithm: {Seg.result["meta"]["intervall_algorithm"]}, '
                                                          f'itterations: {Seg.result["meta"]["N_itter"]}')

            tile_curr = hc_plt.Tile(i, x_label='Returnperiod [years]',
                                    y_label=gl.alias(Seg.colnames['x'], COLNAMES, INPUT["Aliase"]),
                                    title=title,
                                    x_norm='log')

            scatter = hc_plt.Scatter(x=Seg.result["points"]["T_R_x_max"].values,
                                     y=Seg.result["points"]["x_max"].values,
                                     color='red',
                                     size=1)

            tile_curr.add_scatter(scatter)

            Line = hc_plt.Line(x=Seg.result["T_return"]["T_R_grid"].values,
                               y=Seg.result["T_return"]["band_up"].values,
                               color='black',
                               linestyle=':')

            tile_curr.add_line(Line)

            Line = hc_plt.Line(x=Seg.result["T_return"]["T_R_grid"].values,
                               y=Seg.result["T_return"]["band_down"].values,
                               color='black',
                               linestyle='--')

            tile_curr.add_line(Line)

            Line = hc_plt.Line(x=Seg.result["T_return"]["T_R_grid"].values,
                               y=Seg.result["T_return"]["middle"].values,
                               color='black',
                               linestyle='-')

            tile_curr.add_line(Line)

            errorlims = np.vstack(
                (np.array(Seg.result['T_return_single']['middle'] - Seg.result['T_return_single']['down']),
                 np.array(Seg.result['T_return_single']['up'] - Seg.result['T_return_single']['middle'])))
            Errorbar = hc_plt.ErrorBar(x=Seg.result['T_return_single']['T_Return'].values,
                                       y=Seg.result['T_return_single']['middle'].values,
                                       errorlims=errorlims,
                                       color='blue',
                                       fmt='None')

            tile_curr.add_errorbar(Errorbar)

            scatter = hc_plt.Scatter(x=Seg.result['T_return_single']['T_Return'].values,
                                     y=Seg.result['T_return_single']['middle'].values,
                                     color='blue')

            tile_curr.add_scatter(scatter)

            T_R_text = Seg.result["T_return_single"]
            T_R_text.iloc[:, 1:] = gl.significant_digits(T_R_text.values[:, 1:], 3).astype(float)
            new_order = ['T_Return', 'down', 'middle', 'up']
            T_R_text = T_R_text[new_order]

            Textbox = hc_plt.Textbox(data=T_R_text,
                                     fontsize=7,
                                     corner1=(0.4, 0.4),
                                     corner2=(0.9, 0.1))

            tile_curr.add_textbox(Textbox)

            if Seg.angles is not None:
                Tiles.append(tile_curr)
            else:
                Tiles_omni.append(tile_curr)

        FIG_direc = hc_plt.plot_tiled(Tiles, global_max=[None, None], global_min=[None, None], grid=[3, 2], figsize=figsize_fullpage)

        FIG_omni = hc_plt.plot_tiled(Tiles_omni, global_max=[None, None], global_min=[None, None], grid=[1, 1], figsize=figsize_halfpage)

        if 'png' in INPUT["Toggle_Modules"]["plot_as"]:
            gl.save_figs_as_png(FIG_direc + FIG_omni, path_out + f'T_return_{Calc_name}', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

        if 'pdf' in INPUT["Toggle_Modules"]["plot_as"]:
            gl.save_figs_as_pdf(FIG_direc + FIG_omni, path_out + f'T_return_{Calc_name}', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

if 'wind' in INPUT["Toggle_Modules"].get("plot_Validation", {}):
    print('plotting Validation wind...')

    if not ('wind' in INPUT["Toggle_Modules"].get("calc_Validation", {})):
        print("   please toggle calculation to plot")
    else:
        # lines
        Calc = DATA_OUT["Validation"]["wind"]

        Tiles = []
        Tiles_omni = []
        titels = Calc.create_segment_title()

        titels = gl.alias(titels, COLNAMES, INPUT["Aliase"])

        for i, Seg in enumerate(Calc.result):

            cmap_lines = LinearSegmentedColormap.from_list(
                "custom_colormap", ['#ff0000', '#00ff00'])
            range_colors = np.linspace(0, 1, len(Seg.result["condensed"]['vm_vise'].columns))

            barcolor = [0.5, 0.5, 0.5]

            Seg.indizes = pd.to_datetime(Seg.indizes)

            Meta = Seg.result["meta"]

            Legend = [{"color": 'grey', "label": "hindcast", "linestyle": '-'},
                      {"color": 'grey', "label": "condensed", "linestyle": '--'}]

            tile_curr = hc_plt.Tile(i,
                                    x_label=gl.alias(Seg.colnames["Hindcast"]["v_m"], COLNAMES, INPUT["Aliase"]),
                                    y_label=f'Bending DEL [Nm]' + r" $\vert$ " + f'm={Meta["SN_slope"]}' + "\n" + f'N_ref={Meta["N_ref"]:.2e}' + r" $\vert$ " + f'lifetime={Meta["design_life"]}y',
                                    y_label_right='number of datapoints',
                                    title=titels[i],
                                    spinecolor_right=barcolor,
                                    legend=Legend)

            configs = list(Seg.result["condensed"]['vm_vise'].columns)
            configs = [config for config in configs if config != 'count']

            textbox = []
            colors = []

            for i, config in enumerate(configs):
                Line_condensed = hc_plt.Line(x=Seg.result["condensed"]['vm_vise'][config].index,
                                             y=Seg.result["condensed"]['vm_vise'][config].values,
                                             color=cmap_lines(range_colors[i]),
                                             linestyle='--')
                tile_curr.add_line(Line_condensed)

                Line_hindcast = hc_plt.Line(x=Seg.result["hindcast"]['vm_vise'][config].index,
                                            y=Seg.result["hindcast"]['vm_vise'][config].values,
                                            color=cmap_lines(range_colors[i]),
                                            linestyle='-')

                tile_curr.add_line(Line_hindcast)

                # textbox
                textbox.append(str(config))
                textbox.append(f"DEL Hindcast: {Seg.result['hindcast']['added'][config].values[0]:.3e}")
                textbox.append(f"DEL Condensed: {Seg.result['condensed']['added'][config].values[0]:.3e}")
                textbox.append(
                    f"Condensed/Hindcast: {round(Seg.result['condensed']['added'][config].values[0] / Seg.result['hindcast']['added'][config].values[0] * 100, 1)}" + r"\%")
                colors.append([cmap_lines(range_colors[i])])
                colors.append(['black'])
                colors.append(['black'])
                colors.append(['black'])

            Textbox_data = pd.DataFrame(data=textbox)
            colors_data = pd.DataFrame(data=colors)

            Bar_count = hc_plt.Bar(x=Seg.result["condensed"]['vm_vise'][configs[0]].index,
                                   y=Seg.result["condensed"]['vm_vise']["count"],
                                   width=1,
                                   bottom=None,
                                   align='center',
                                   color=barcolor,
                                   alpha=0.5,
                                   spinecolor=barcolor,
                                   yy_side='right')

            tile_curr.add_bar(Bar_count)

            if Seg.angles is not None:
                Textbox_DEL = hc_plt.Textbox(Textbox_data,
                                             fontsize=7,
                                             corner1=[0.4, 1],
                                             corner2=[1, 0.4],
                                             colors=colors_data,
                                             orientation_h='left',
                                             orientation_v='center',
                                             header=False)

                tile_curr.add_textbox(Textbox_DEL)
                Tiles.append(tile_curr)

            else:
                Textbox_DEL = hc_plt.Textbox(Textbox_data,
                                             fontsize=8,
                                             corner1=[0.6, 1],
                                             corner2=[1, 0.4],
                                             colors=colors_data,
                                             orientation_h='left',
                                             orientation_v='center',
                                             header=False)

                tile_curr.add_textbox(Textbox_DEL)
                Tiles_omni.append(tile_curr)

        FIG_direc = hc_plt.plot_tiled(Tiles, global_max=['auto', 'auto'], global_min=[0, 0], grid=[3, 2], figsize=figsize_fullpage)
        FIG_omni = hc_plt.plot_tiled(Tiles_omni, global_max=['auto', 'auto'], global_min=[0, 0], grid=[1, 1], figsize=figsize_halfpage)

        if 'png' in INPUT["Toggle_Modules"]["plot_as"]:
            gl.save_figs_as_png(FIG_direc + FIG_omni, path_out + 'Valid_line_wind', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

        if 'pdf' in INPUT["Toggle_Modules"]["plot_as"]:
            gl.save_figs_as_pdf(FIG_direc + FIG_omni, path_out + 'Valid_line_wind', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

        # scatter
        df_DEL = Calc.load_from_db()

        for config in INPUT["Validation_wind"]["scatter_configs"]:
            Tiles = []
            Tiles_omni = []
            titels = Calc.create_segment_title()

            titels = [f'Wind Sea, config: {config} ' + "\n" + title + "\n" for title in titels]

            titels = gl.alias(titels, COLNAMES, INPUT["Aliase"])

            for i, Seg in enumerate(Calc.result):
                Meta = Seg.result["meta"]

                df_data = gl.export_df_from_sql(Calc.basedata["dbname"],
                                                "Hindcast_combined",
                                                column_names=[Seg.colnames["Hindcast"]["H_s"], Seg.colnames["Hindcast"]["T_p"]],
                                                indizes=Seg.indizes)

                tile_curr = hc_plt.Tile(i,
                                        x_label=gl.alias(Seg.colnames["Hindcast"]["H_s"], COLNAMES, INPUT["Aliase"]),
                                        y_label=gl.alias(Seg.colnames["Hindcast"]["T_p"], COLNAMES, INPUT["Aliase"]),
                                        title=titels[i])

                scatter = hc_plt.Scatter(x=df_data[Seg.colnames["Hindcast"]["H_s"]],
                                         y=df_data[Seg.colnames["Hindcast"]["T_p"]],
                                         cmap='cool',
                                         cmap_mode='manual',
                                         c=df_DEL.loc[Seg.indizes, config].values,
                                         size=2,
                                         cbar=True,
                                         cbar_label=f'Bending DEL [Nm]' + r" $\vert$ " + f'm={Meta["SN_slope"]}' + "\n" + f'N_ref={Meta["N_ref"]:.2e}' + r" $\vert$ " + f'lifetime={Meta["design_life"]}y')

                tile_curr.add_scatter(scatter)

                Line_f0 = hc_plt.Line(x=None,
                                      y=[1 / INPUT["Structure"]["f_0"]],
                                      label=f'f_0, f={INPUT["Structure"]["f_0"]} Hz',
                                      color='green',
                                      linestyle=':')
                tile_curr.add_line(Line_f0)

                if Seg.angles is not None:
                    Tiles.append(tile_curr)
                else:
                    Tiles_omni.append(tile_curr)

            FIG_direc = hc_plt.plot_tiled(Tiles, global_max=['auto', 'auto'], global_min=[0, 0], grid=[3, 2], figsize=figsize_fullpage)
            FIG_omni = hc_plt.plot_tiled(Tiles_omni, global_max=['auto', 'auto'], global_min=[0, 0], grid=[1, 1], figsize=figsize_halfpage)

            if 'png' in INPUT["Toggle_Modules"]["plot_as"]:
                gl.save_figs_as_png(FIG_direc + FIG_omni, path_out + f'Valid_scatter_wind_{config}', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

            if 'pdf' in INPUT["Toggle_Modules"]["plot_as"]:
                gl.save_figs_as_pdf(FIG_direc + FIG_omni, path_out + f'Valid_scatter_wind_{config}', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

if 'swell' in INPUT["Toggle_Modules"].get("plot_Validation", {}):
    print('plotting Validation swell...')

    if not ('swell' in INPUT["Toggle_Modules"].get("calc_Validation", {})):
        print("   please toggle calculation to plot")
    else:
        # lines
        Calc = DATA_OUT["Validation"]["swell"]

        Tiles = []
        Tiles_omni = []
        titels = Calc.create_segment_title()

        titels = gl.alias(titels, COLNAMES, INPUT["Aliase"])

        for i, Seg in enumerate(Calc.result):

            cmap_lines = LinearSegmentedColormap.from_list(
                "custom_colormap", ['#ff0000', '#00ff00'])
            range_colors = np.linspace(0, 1, len(Seg.result["condensed"]['vm_vise'].columns))

            barcolor = [0.5, 0.5, 0.5]

            Seg.indizes = pd.to_datetime(Seg.indizes)

            Meta = Seg.result["meta"]

            Legend = [{"color": 'grey', "label": "hindcast", "linestyle": '-'},
                      {"color": 'grey', "label": "condensed", "linestyle": '--'}]

            tile_curr = hc_plt.Tile(i,
                                    x_label=gl.alias(Seg.colnames["Hindcast"]["v_m"], COLNAMES, INPUT["Aliase"]),
                                    y_label=f'Bending DEL [Nm]' + r" $\vert$ " + f'm={Meta["SN_slope"]}' + "\n" + f'N_ref={Meta["N_ref"]:.2e}' + r" $\vert$ " + f'lifetime={Meta["design_life"]}y',
                                    y_label_right='number of datapoints',
                                    title=titels[i],
                                    spinecolor_right=barcolor,
                                    legend=Legend)

            configs = list(Seg.result["condensed"]['vm_vise'].columns)
            configs = [config for config in configs if config != 'count']

            textbox = []
            colors = []

            for i, config in enumerate(configs):
                Line_condensed = hc_plt.Line(x=Seg.result["condensed"]['vm_vise'][config].index,
                                             y=Seg.result["condensed"]['vm_vise'][config].values,
                                             color=cmap_lines(range_colors[i]),
                                             linestyle='--')
                tile_curr.add_line(Line_condensed)

                Line_hindcast = hc_plt.Line(x=Seg.result["hindcast"]['vm_vise'][config].index,
                                            y=Seg.result["hindcast"]['vm_vise'][config].values,
                                            color=cmap_lines(range_colors[i]),
                                            linestyle='-')

                tile_curr.add_line(Line_hindcast)

                # textbox
                textbox.append(str(config))
                textbox.append(f"DEL Hindcast: {Seg.result['hindcast']['added'][config].values[0]:.3e}")
                textbox.append(f"DEL Condensed: {Seg.result['condensed']['added'][config].values[0]:.3e}")
                textbox.append(
                    f"Condensed/Hindcast: {round(Seg.result['condensed']['added'][config].values[0] / Seg.result['hindcast']['added'][config].values[0] * 100, 1)}" + r"\%")
                colors.append([cmap_lines(range_colors[i])])
                colors.append(['black'])
                colors.append(['black'])
                colors.append(['black'])

            Textbox_data = pd.DataFrame(data=textbox)
            colors_data = pd.DataFrame(data=colors)

            Bar_count = hc_plt.Bar(x=Seg.result["condensed"]['vm_vise'][configs[0]].index,
                                   y=Seg.result["condensed"]['vm_vise']["count"],
                                   width=1,
                                   bottom=None,
                                   align='center',
                                   color=barcolor,
                                   alpha=0.5,
                                   spinecolor=barcolor,
                                   yy_side='right')

            tile_curr.add_bar(Bar_count)

            if Seg.angles is not None:
                Textbox_DEL = hc_plt.Textbox(Textbox_data,
                                             fontsize=6,
                                             corner1=[0.4, 1],
                                             corner2=[1, 0.4],
                                             colors=colors_data,
                                             orientation_h='left',
                                             orientation_v='center',
                                             header=False)

                tile_curr.add_textbox(Textbox_DEL)
                Tiles.append(tile_curr)

            else:
                Textbox_DEL = hc_plt.Textbox(Textbox_data,
                                             fontsize=8,
                                             corner1=[0.6, 1],
                                             corner2=[1, 0.4],
                                             colors=colors_data,
                                             orientation_h='left',
                                             orientation_v='center',
                                             header=False)

                tile_curr.add_textbox(Textbox_DEL)
                Tiles_omni.append(tile_curr)

        FIG_direc = hc_plt.plot_tiled(Tiles, global_max=['auto', 'auto'], global_min=[0, 0], grid=[3, 2], figsize=figsize_fullpage)
        FIG_omni = hc_plt.plot_tiled(Tiles_omni, global_max=['auto', 'auto'], global_min=[0, 0], grid=[1, 1], figsize=figsize_halfpage)

        if 'png' in INPUT["Toggle_Modules"]["plot_as"]:
            gl.save_figs_as_png(FIG_direc + FIG_omni, path_out + 'Valid_line_swell', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

        if 'pdf' in INPUT["Toggle_Modules"]["plot_as"]:
            gl.save_figs_as_pdf(FIG_direc + FIG_omni, path_out + 'Valid_line_swell', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

        # scatter
        df_DEL = Calc.load_from_db()

        for config in INPUT["Validation_swell"]["scatter_configs"]:
            Tiles = []
            Tiles_omni = []
            titels = Calc.create_segment_title()

            titels = [f'Swell Sea, config: {config} ' + "\n" + title + "\n" for title in titels]

            titels = gl.alias(titels, COLNAMES, INPUT["Aliase"])

            for i, Seg in enumerate(Calc.result):
                Meta = Seg.result["meta"]

                df_data = gl.export_df_from_sql(Calc.basedata["dbname"],
                                                "Hindcast_combined",
                                                column_names=[Seg.colnames["Hindcast"]["H_s"], Seg.colnames["Hindcast"]["T_p"]],
                                                indizes=Seg.indizes)

                tile_curr = hc_plt.Tile(i,
                                        x_label=gl.alias(Seg.colnames["Hindcast"]["H_s"], COLNAMES, INPUT["Aliase"]),
                                        y_label=gl.alias(Seg.colnames["Hindcast"]["T_p"], COLNAMES, INPUT["Aliase"]),
                                        title=titels[i])

                scatter = hc_plt.Scatter(x=df_data[Seg.colnames["Hindcast"]["H_s"]],
                                         y=df_data[Seg.colnames["Hindcast"]["T_p"]],
                                         cmap='cool',
                                         cmap_mode='manual',
                                         c=df_DEL.loc[Seg.indizes, config].values,
                                         size=2,
                                         cbar=True,
                                         cbar_label=f'Bending DEL [Nm]' + r" $\vert$ " + f'm={Meta["SN_slope"]}' + "\n " + f'N_ref={Meta["N_ref"]:.2e}' + r" $\vert$ " + f'lifetime={Meta["design_life"]}y')

                tile_curr.add_scatter(scatter)

                Line_f0 = hc_plt.Line(x=None,
                                      y=[1 / INPUT["Structure"]["f_0"]],
                                      label=f'f_0, f={INPUT["Structure"]["f_0"]} Hz',
                                      color='green',
                                      linestyle=':')
                tile_curr.add_line(Line_f0)

                if Seg.angles is not None:
                    Tiles.append(tile_curr)
                else:
                    Tiles_omni.append(tile_curr)

            FIG_direc = hc_plt.plot_tiled(Tiles, global_max=['auto', 'auto'], global_min=[0, 0], grid=[3, 2], figsize=figsize_fullpage)
            FIG_omni = hc_plt.plot_tiled(Tiles_omni, global_max=['auto', 'auto'], global_min=[0, 0], grid=[1, 1], figsize=figsize_halfpage)

            if 'png' in INPUT["Toggle_Modules"]["plot_as"]:
                gl.save_figs_as_png(FIG_direc + FIG_omni, path_out + f'Valid_scatter_swell_{config}', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

            if 'pdf' in INPUT["Toggle_Modules"]["plot_as"]:
                gl.save_figs_as_pdf(FIG_direc + FIG_omni, path_out + f'Valid_scatter_swell_{config}', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

if INPUT["Toggle_Modules"].get("plot_SensorEval", {}):
    print("plotting SensorEval...")

    for Calc_name, Calc in DATA_OUT["SensorEval"].items():

        Tiles = []
        Tiles_omni = []

        titels = Calc.create_segment_title(mode='sparse')

        titels = gl.alias(titels, COLNAMES, INPUT["Aliase"])

        for i, Seg in enumerate(Calc.result):

            titel = f'Histogramm with binsize={Seg.result["bin_size"]}, ' + titels[i]

            Seg.indizes = pd.to_datetime(Seg.indizes)

            # tile histo
            tile_histo = hc_plt.Tile(i,
                                     x_label=gl.alias(Seg.colnames['x'], COLNAMES, INPUT["Aliase"]),
                                     y_label='number of datapoints [-]',
                                     title=titel)

            bars_histo = hc_plt.Bar(x=Seg.result["center"],
                                    y=Seg.result["count"],
                                    color='grey',
                                    width=Seg.result["bin_size"])

            tile_histo.add_bar(bars_histo)

            # tile timeseries
            df = Calc.load_from_db(colnames_ini=True, indizes=Seg.indizes)
            x = df[Seg.colnames['x']].values
            titel = f'Timeseries with min={min(x)}' + r" $\vert$ " + f'max={max(x)}' + r" $\vert$ " + f'standard deviation={round(np.std(x), 4)}, ' + titels[i]

            tile_time = hc_plt.Tile(i, x_label='date', y_label=gl.alias(Seg.colnames['x'], COLNAMES, INPUT["Aliase"]), title=titel)

            timeseries = hc_plt.Line(x=df[Seg.colnames['x']].index,
                                     y=df[Seg.colnames['x']].values,
                                     color='grey',
                                     linewidth=1)

            tile_time.add_line(timeseries)

            if Seg.angles is not None:
                Tiles.append(tile_histo)
                Tiles.append(tile_time)
            else:
                Tiles_omni.append(tile_histo)
                Tiles_omni.append(tile_time)

        # FIG_direc = hc_plt.plot_tiled(Tiles, global_max=['auto', 'auto'], global_min=[0, 0], grid=[3, 2])

        FIG_omni = hc_plt.plot_tiled(Tiles_omni, global_max=[None, None], global_min=[None, None], grid=[2, 1], figsize=figsize_halfpage)

        if 'png' in INPUT["Toggle_Modules"]["plot_as"]:
            gl.save_figs_as_png(FIG_omni, path_out + f'SensorEval_{Calc_name}', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

        if 'pdf' in INPUT["Toggle_Modules"]["plot_as"]:
            gl.save_figs_as_pdf(FIG_omni, path_out + f'SensorEval_{Calc_name}', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

if INPUT["Toggle_Modules"].get("plot_Weibull", {}):
    print("plotting Weibull fit...")

    for Calc_name, Calc in DATA_OUT["Weibull"].items():

        Tiles = []
        Tiles_omni = []

        titels = Calc.create_segment_title(mode='sparse')

        titels = gl.alias(titels, COLNAMES, INPUT["Aliase"])

        for i, Seg in enumerate(Calc.result):

            titel = f'Weibull fit with binsize={Seg.result["bin_size"]}, ' + "\n" + titels[i]

            Seg.indizes = pd.to_datetime(Seg.indizes)

            # tile histo
            tile_curr = hc_plt.Tile(i,
                                    x_label=gl.alias(Seg.colnames['x'], COLNAMES, INPUT["Aliase"]),
                                    y_label='probability density [-]',
                                    title=titel)

            bars_histo = hc_plt.Bar(x=Seg.result["center"],
                                    y=Seg.result["prob"],
                                    color='grey',
                                    width=Seg.result["bin_size"])

            tile_curr.add_bar(bars_histo)

            line_weibull = hc_plt.Line(x=Seg.result["weibull"].index,
                                       y=Seg.result["weibull"].values,
                                       label=f'Weibull fit')

            tile_curr.add_line(line_weibull)

            textbox_data = pd.DataFrame(data=[[key, round(Seg.result["weibull_params"][key], 2)] for key in Seg.result["weibull_params"].keys()])

            textbox_params = hc_plt.Textbox(data=textbox_data, corner1=[0.6, 0.95], corner2=[0.9, 0.45], header=False, orientation_h='left')
            tile_curr.add_textbox(textbox_params)

            if Seg.angles is not None:
                Tiles.append(tile_curr)
            else:
                Tiles_omni.append(tile_curr)

        FIG_direc = hc_plt.plot_tiled(Tiles, global_max=[None, None], global_min=[0, 0], grid=[3, 2], figsize=figsize_fullpage)

        FIG_omni = hc_plt.plot_tiled(Tiles_omni, global_max=[None, None], global_min=[None, None], grid=[1, 1], figsize=figsize_halfpage)

        if 'png' in INPUT["Toggle_Modules"]["plot_as"]:
            gl.save_figs_as_png(FIG_direc + FIG_omni, path_out + f'Weibull_{Calc_name}', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

        if 'pdf' in INPUT["Toggle_Modules"]["plot_as"]:
            gl.save_figs_as_pdf(FIG_direc + FIG_omni, path_out + f'Weibull_{Calc_name}', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

if INPUT["Toggle_Modules"].get("plot_ExtremeConture", {}) and len(INPUT["Toggle_Modules"].get("calc_ExtremeConture", {})) > 0:
    print('plotting ExtremeConture...')

    for name, Calc in DATA_OUT["ExtremeConture"].items():

        Input = INPUT["ExtremeValues"]
        Tiles = []
        Tiles_omni = []

        df = Calc.load_from_db(colnames_ini=True)
        titels = Calc.create_segment_title()
        titels = gl.alias(titels, COLNAMES, INPUT["Aliase"])

        for i, Seg in enumerate(Calc.result):

            Seg.indizes = pd.to_datetime(Seg.indizes)
            point_data = df[df.index.isin(Seg.indizes)]

            tile_curr = hc_plt.Tile(i,
                                    x_label=gl.alias(Seg.colnames['x'], COLNAMES, INPUT["Aliase"]),
                                    y_label=gl.alias(Seg.colnames['y'], COLNAMES, INPUT["Aliase"]),
                                    title=titels[i])

            scatter = hc_plt.Scatter(x=point_data[Seg.colnames["x"]],
                                     y=point_data[Seg.colnames["y"]],
                                     cmap='cool',
                                     size=2,
                                     cmap_norm='sqrt')

            tile_curr.add_scatter(scatter)

            # add contures
            # Define the colors for the colormap: red to dark green
            colors = [(1, 0, 0),  # Red
                      (0, 0.4, 0)]  # Dark green
            # Create the colormap
            cmap_name = 'red_to_darkgreen'
            red_to_darkgreen = LinearSegmentedColormap.from_list(cmap_name, colors)

            color = np.linspace(1, 0, len(Seg.result))
            i = 0
            for name, data in Seg.result.items():

                contour = hc_plt.Line(x=data["x"],
                                      y=data["y"],
                                      label=name,
                                      color=red_to_darkgreen(color[i]),
                                      linewidth=0.8)

                tile_curr.add_line(contour)

                i = i + 1
            if Seg.angles is not None:
                Tiles.append(tile_curr)

            else:
                Tiles_omni.append(tile_curr)


        FIG_direc = hc_plt.plot_tiled(Tiles, global_max=['auto', 'auto'], global_min=[0, 0], grid=[3, 2], figsize=figsize_fullpage)

        FIG_omni = hc_plt.plot_tiled(Tiles_omni, global_max=['auto', 'auto'], global_min=[0, 0], grid=[1, 1], figsize=figsize_halfpage)

        if 'png' in INPUT["Toggle_Modules"]["plot_as"]:
            gl.save_figs_as_png(FIG_direc + FIG_omni, path_out + 'ExtremeConture', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

        if 'pdf' in INPUT["Toggle_Modules"]["plot_as"]:
            gl.save_figs_as_pdf(FIG_direc + FIG_omni, path_out + 'ExtremeConture', dpi=INPUT["Toggle_Modules"]["dpi_figures"])


# %% Data Out
if INPUT["DataOut"]["CSV_out"]:
    print("saving CSV Data...")

    path_csv = os.path.join(path_out, 'csv_data')

    try:
        # Create the new folder
        os.makedirs(path_csv, exist_ok=True)  # exist_ok=True prevents an error if the folder already exists
        print(f"Folder created successfully at: {path_csv}")
    except Exception as e:
        print(f"An error occurred: {e}")

    unpack_funcs = {"flat_data": lambda x: [segment.result for segment in x.result],
                    "flat_angles": lambda x: [segment.angles for segment in x.result],
                    "flat_exclude_omni": lambda x: [segment.result for segment in x.result if segment.angles is not None],
                    }

    if "wind" in DATA_OUT["VMHS"]:
        print("   VMHS wind")
        calc = DATA_OUT["VMHS"]["wind"]
        data = unpack_funcs["flat_data"](calc)
        table_names = unpack_funcs["flat_angles"](calc)
        table_names = ["omnidirectional" if name is None else f"{name[0]} to {name[1]}" for name in table_names]

        gl.save_df_list_to_excel(path_csv + r'//VMHS_wind', data, sheet_names=table_names)

    if "swell" in DATA_OUT["VMHS"]:
        print("   VMHS swell")
        calc = DATA_OUT["VMHS"]["swell"]
        data = unpack_funcs["flat_data"](calc)
        table_names = unpack_funcs["flat_angles"](calc)
        table_names = ["omnidirectional" if name is None else f"{name[0]} to {name[1]}" for name in table_names]

        gl.save_df_list_to_excel(path_csv + r'//VMHS_swell', data, sheet_names=table_names)

    if "wind" in DATA_OUT["HSTP"]:
        print("   HSTP wind")
        calc = DATA_OUT["HSTP"]["wind"]
        data = unpack_funcs["flat_data"](calc)
        table_names = unpack_funcs["flat_angles"](calc)
        table_names = ["omnidirectional" if name is None else f"{name[0]} to {name[1]}" for name in table_names]

        gl.save_df_list_to_excel(path_csv + r'//HSTP_wind', data, sheet_names=table_names)

    if "swell" in DATA_OUT["HSTP"]:
        print("   HSTP swell")
        calc = DATA_OUT["HSTP"]["swell"]
        data = unpack_funcs["flat_data"](calc)
        table_names = unpack_funcs["flat_angles"](calc)
        table_names = ["omnidirectional" if name is None else f"{name[0]} to {name[1]}" for name in table_names]

        gl.save_df_list_to_excel(path_csv + r'//HSTP_swell', data, sheet_names=table_names)

    if "wind" in DATA_OUT["VMTP"]:
        print("   VMTP wind")
        calc = DATA_OUT["VMTP"]["wind"]
        data = unpack_funcs["flat_data"](calc)
        table_names = unpack_funcs["flat_angles"](calc)
        table_names = ["omnidirectional" if name is None else f"{name[0]} to {name[1]}" for name in table_names]

        gl.save_df_list_to_excel(path_csv + r'//VMTP_wind', data, sheet_names=table_names)

    if "swell" in DATA_OUT["VMTP"]:
        print("   VMTP swell")
        calc = DATA_OUT["VMTP"]["swell"]
        data = unpack_funcs["flat_data"](calc)
        table_names = unpack_funcs["flat_angles"](calc)
        table_names = ["omnidirectional" if name is None else f"{name[0]} to {name[1]}" for name in table_names]

        gl.save_df_list_to_excel(path_csv + r'//VMTP_swell', data, sheet_names=table_names)

    if "wind" in DATA_OUT["RWI"]:
        print("   RWI wind")
        calc = DATA_OUT["RWI"]["wind"]
        df = calc.load_from_db(colnames_ini=True)
        data = unpack_funcs["flat_data"](calc)

        for i, segment in enumerate(calc.result):
            df_temp = df.loc[segment.indizes]

            data[i] = pd.concat((data[i], df_temp), axis=1)

        table_names = unpack_funcs["flat_angles"](calc)
        table_names = ["omnidirectional" if name is None else f"{name[0]} to {name[1]}" for name in table_names]

        gl.save_df_list_to_excel(path_csv + r'//RWI_wind', data, sheet_names=table_names)

    if "total" in DATA_OUT["RWI"]:
        print("   RWI total")
        calc = DATA_OUT["RWI"]["total"]
        df = calc.load_from_db(colnames_ini=True)
        data = unpack_funcs["flat_data"](calc)

        for i, segment in enumerate(calc.result):
            df_temp = df.loc[segment.indizes]

            data[i] = pd.concat((data[i], df_temp), axis=1)

        table_names = unpack_funcs["flat_angles"](calc)
        table_names = ["omnidirectional" if name is None else f"{name[0]} to {name[1]}" for name in table_names]

        gl.save_df_list_to_excel(path_csv + r'//RWI_swell', data, sheet_names=table_names)

    if "wind" in DATA_OUT["WaveBreak_Steep"]:
        print("   WaveBreak_Steep wind")
        calc = DATA_OUT["WaveBreak_Steep"]["wind"]
        data = unpack_funcs["flat_data"](calc)
        table_names = unpack_funcs["flat_angles"](calc)
        table_names = ["omnidirectional" if name is None else f"{name[0]} to {name[1]}" for name in table_names]

        gl.save_df_list_to_excel(path_csv + r'//WaveBreak_Steep_wind', data, sheet_names=table_names)

    if "wind" in DATA_OUT["table_vmhs"]:
        print("   table_vmhs wind")
        calc = DATA_OUT["table_vmhs"]["wind"]
        data = unpack_funcs["flat_data"](calc)
        table_names = unpack_funcs["flat_angles"](calc)
        table_names = ["omnidirectional" if name is None else f"{name[0]} to {name[1]}" for name in table_names]

        # combined
        combined = []
        for table in data:
            combined.append(table["value"])

        combined = np.array(combined).T

        df_combined = pd.DataFrame(combined, columns=table_names, index=data[0]["vm_edges"])

        data.insert(0, df_combined)
        table_names.insert(0, 'combined')
        gl.save_df_list_to_excel(path_csv + r'//table_vmhs_wind', data, sheet_names=table_names)

    if "swell" in DATA_OUT["table_vmhs"]:
        print("   table_vmh sswell")
        calc = DATA_OUT["table_vmhs"]["swell"]
        data = unpack_funcs["flat_data"](calc)
        table_names = unpack_funcs["flat_angles"](calc)
        table_names = ["omnidirectional" if name is None else f"{name[0]} to {name[1]}" for name in table_names]

        # combined
        combined = []
        for table in data:
            combined.append(table["value"])

        combined = np.array(combined).T

        df_combined = pd.DataFrame(combined, columns=table_names, index=data[0]["vm_edges"])

        data.insert(0, df_combined)
        table_names.insert(0, 'combined')
        gl.save_df_list_to_excel(path_csv + r'//table_vmhs_swell', data, sheet_names=table_names)

    if "wind" in DATA_OUT["table_vmtp"]:
        print("   table_vmtp wind")
        calc = DATA_OUT["table_vmtp"]["wind"]
        data = unpack_funcs["flat_data"](calc)
        table_names = unpack_funcs["flat_angles"](calc)
        table_names = ["omnidirectional" if name is None else f"{name[0]} to {name[1]}" for name in table_names]

        # combined
        combined = []
        for table in data:
            combined.append(table["value"])

        combined = np.array(combined).T

        df_combined = pd.DataFrame(combined, columns=table_names, index=data[0]["vm_edges"])

        data.insert(0, df_combined)
        table_names.insert(0, 'combined')
        gl.save_df_list_to_excel(path_csv + r'//table_vmtp_wind', data, sheet_names=table_names)

    if "swell" in DATA_OUT["table_vmtp"]:
        print("   table_vmtp swell")
        calc = DATA_OUT["table_vmtp"]["swell"]
        data = unpack_funcs["flat_data"](calc)
        table_names = unpack_funcs["flat_angles"](calc)
        table_names = ["omnidirectional" if name is None else f"{name[0]} to {name[1]}" for name in table_names]

        # combined
        combined = []
        for table in data:
            combined.append(table["value"])

        combined = np.array(combined).T

        df_combined = pd.DataFrame(combined, columns=table_names, index=data[0]["vm_edges"])

        data.insert(0, df_combined)
        table_names.insert(0, 'combined')
        gl.save_df_list_to_excel(path_csv + r'//table_vmtp_swell', data, sheet_names=table_names)

    if DATA_OUT["AngleDeviation"]:
        print("   AngleDeviation")
        calc = DATA_OUT["AngleDeviation"]
        data = unpack_funcs["flat_exclude_omni"](calc)
        table_names = unpack_funcs["flat_angles"](calc)
        table_names = table_names[1:]
        table_names = [f"{name[0]} to {name[1]}" for name in table_names]
        gl.save_df_list_to_excel(path_csv + r'//AngleDeviation', data, sheet_names=table_names)

    if "wind" in DATA_OUT["Validation"]:
        print("   Validation wind")
        calc = DATA_OUT["Validation"]["wind"]

        unpac_func_hindcast_vm_vise = lambda x: [df.result["hindcast"]["vm_vise"] for df in x.result]
        unpac_func_condensed_vm_vise = lambda x: [df.result["condensed"]["vm_vise"] for df in x.result]
        unpac_func_condensed_added = lambda x: [df.result["condensed"]["added"] for df in x.result]
        unpac_func_hindcast_added = lambda x: [df.result["hindcast"]["added"] for df in x.result]

        data_hindcast_vm_vise = unpac_func_hindcast_vm_vise(calc)
        data_condensed_vm_vise = unpac_func_condensed_vm_vise(calc)
        data_condensed_added = unpac_func_condensed_added(calc)
        data_hindcast_added = unpac_func_hindcast_added(calc)

        table_names = unpack_funcs["flat_angles"](calc)
        table_names = ["omnidirectional" if name is None else f"{name[0]} to {name[1]}" for name in table_names]

        table_hindcast_added_combined = [list(table.values[0]) for table in data_hindcast_added]
        table_hindcast_added_combined = pd.DataFrame(table_hindcast_added_combined, columns=data_hindcast_added[0].columns)
        table_hindcast_added_combined["angles"] = table_names

        table_condensed_added_combined = [list(table.values[0]) for table in data_condensed_added]
        table_condensed_added_combined = pd.DataFrame(table_condensed_added_combined, columns=data_condensed_added[0].columns)
        table_condensed_added_combined["angles"] = table_names

        table_wind_added = [table_hindcast_added_combined, table_condensed_added_combined]

        gl.save_df_list_to_excel(path_csv + r'//Validation_wind_hindcast_vm_vise', data_hindcast_vm_vise, sheet_names=table_names)
        gl.save_df_list_to_excel(path_csv + r'//Validation_wind_condensed_vm_vise', data_condensed_vm_vise, sheet_names=table_names)
        gl.save_df_list_to_excel(path_csv + r'//Validation_wind_added', table_wind_added, sheet_names=["hindcast", "condensed"])

    if "swell" in DATA_OUT["Validation"]:
        print("   Validation swell")
        calc = DATA_OUT["Validation"]["swell"]

        unpac_func_hindcast_vm_vise = lambda x: [df.result["hindcast"]["vm_vise"] for df in x.result]
        unpac_func_condensed_vm_vise = lambda x: [df.result["condensed"]["vm_vise"] for df in x.result]
        unpac_func_condensed_added = lambda x: [df.result["condensed"]["added"] for df in x.result]
        unpac_func_hindcast_added = lambda x: [df.result["hindcast"]["added"] for df in x.result]

        data_hindcast_vm_vise = unpac_func_hindcast_vm_vise(calc)
        data_condensed_vm_vise = unpac_func_condensed_vm_vise(calc)
        data_hindcast_added = unpac_func_hindcast_added(calc)
        data_condensed_added = unpac_func_condensed_added(calc)

        table_names = unpack_funcs["flat_angles"](calc)
        table_names = ["omnidirectional" if name is None else f"{name[0]} to {name[1]}" for name in table_names]

        table_hindcast_added_combined = [list(table.values[0]) for table in data_hindcast_added]
        table_hindcast_added_combined = pd.DataFrame(table_hindcast_added_combined, columns=data_hindcast_added[0].columns)
        table_hindcast_added_combined["angles"] = table_names

        table_condensed_added_combined = [list(table.values[0]) for table in data_condensed_added]
        table_condensed_added_combined = pd.DataFrame(table_condensed_added_combined, columns=data_condensed_added[0].columns)
        table_condensed_added_combined["angles"] = table_names

        table_wind_added = [table_hindcast_added_combined, table_condensed_added_combined]

        gl.save_df_list_to_excel(path_csv + r'//Validation_swell_hindcast_vm_vise', data_hindcast_vm_vise, sheet_names=table_names)
        gl.save_df_list_to_excel(path_csv + r'//Validation_swell_condensed_vm_vise', data_condensed_vm_vise, sheet_names=table_names)
        gl.save_df_list_to_excel(path_csv + r'//Validation_swell_added', table_wind_added, sheet_names=["hindcast", "condensed"])

# %% MAIN - Save Infolog
print("saving Log...")

with open(path_out + 'Info_LOG.txt', "w") as log_text:
    log_text.write(INFO_LOG)
del log_text

print(f"{script_name} finished!")

#
# def Weibull_method_of_moment(X):
#     import scipy.stats as stats
#     X = X + 0.0001;
#     n = len(X);
#     m1 = np.mean(X);
#     cm1 = np.mean((X - np.mean(X)) ** 1);
#     m2 = np.var(X);
#     cm2 = np.mean((X - np.mean(X)) ** 2);
#     m3 = stats.skew(X);
#     cm3 = np.mean((X - np.mean(X)) ** 3);
#
#     from scipy.special import gamma
#     def m1fun(a, b, c):
#         return a + b * gamma(1 + 1 / c)
#
#     def cm2fun(b, c):
#         return b ** 2 * (gamma(1 + 2 / c) - gamma(1 + 1 / c) ** 2)
#
#     def cm3fun(b, c):
#         return b ** 3 * (gamma(1 + 3 / c) - 3 * gamma(1 + 1 / c) * gamma(1 + 2 / c) + 2 * gamma(1 + 1 / c) ** 3)
#
#     def cfun(c):
#         return abs(np.sqrt(cm3fun(1, c) ** 2 / cm2fun(1, c) ** 3) - np.sqrt(cm3 ** 2 / cm2 ** 3))
#
#     from scipy import optimize
#     cHat = optimize.fminbound(cfun, -2, 5)  # shape
#
#     def bfun(b):
#         return abs(cm2fun(b, cHat) - cm2)
#
#     bHat = optimize.fminbound(bfun, -5, 30)  # scale
#
#     def afun(a):
#         return abs(m1fun(a, bHat, cHat) - m1)
#
#     aHat = optimize.fminbound(afun, -5, 30)  # location
#
#     return cHat, aHat, bHat  # shape, location, scale
#
#
# def joint_distribution_Hs_Tp(data, var_hs='hs', var_tp='tp', periods=None, adjustment=None):
#     """
#     This fuction will plot Hs-Tp joint distribution using LogNoWe model (the Lognormal + Weibull distribution)
#     df : dataframe,
#     var1 : Hs: significant wave height,
#     var2 : Tp: Peak period
#     file_out: Hs-Tp joint distribution, optional
#     """
#     if periods is None:
#         periods = [1, 10, 100, 10000]
#
#     if adjustment == 'NORSOK':
#         periods_adj = np.array([x * 6 for x in periods])
#     else:
#         periods_adj = periods
#
#     df = data
#     pd.options.mode.chained_assignment = None  # default='warn'
#     df.loc[:, 'hs'] = df[var_hs].values
#     df.loc[:, 'tp'] = Tp_correction(df[var_tp].values)
#
#     import scipy.stats as stats
#     from scipy.optimize import curve_fit
#     from scipy.signal import find_peaks
#
#     # calculate lognormal and weibull parameters and plot the PDFs
#     mu = np.mean(np.log(df.hs.values))  # mean of ln(Hs)
#     std = np.std(np.log(df.hs.values))  # standard deviation of ln(Hs)
#     alpha = mu
#     sigma = std
#
#     h = np.linspace(start=0.01, stop=30, num=1500)
#
#     if 0 < mu < 5:
#         pdf_Hs1 = 1 / (np.sqrt(2 * np.pi) * alpha * h) * np.exp(-(np.log(h) - sigma) ** 2 / (2 * alpha ** 2))
#     else:
#         param = stats.lognorm.fit(df.hs.values, )  # shape, loc, scale
#         pdf_lognorm = stats.lognorm.pdf(h, param[0], loc=param[1], scale=param[2])
#         pdf_Hs1 = pdf_lognorm
#
#     param = Weibull_method_of_moment(df.hs.values)  # stats.weibull_min.fit(df.hs.values) # shape, loc, scale
#     pdf_Hs2 = stats.weibull_min.pdf(h, param[0], loc=param[1], scale=param[2])
#
#     # Find the index where two PDF cut, between P60 and P99
#     for i in range(len(h)):
#         if abs(h[i] - np.percentile(df.hs.values, 60)) < 0.1:
#             i1 = i
#
#         if abs(h[i] - np.percentile(df.hs.values, 99)) < 0.1:
#             i2 = i
#
#     epsilon = abs(pdf_Hs1[i1:i2] - pdf_Hs2[i1:i2])
#     param = find_peaks(1 / epsilon)
#     try:
#         index = param[0][1]
#     except:
#         try:
#             index = param[0][0]
#         except:
#             index = np.where(epsilon == epsilon.min())[0]
#     index = index + i1
#
#     # Merge two functions and do smoothing around the cut
#     eta = h[index]
#     pdf_Hs = h * 0
#     for i in range(len(h)):
#         if h[i] < eta:
#             pdf_Hs[i] = pdf_Hs1[i]
#         else:
#             pdf_Hs[i] = pdf_Hs2[i]
#
#     for i in range(len(h)):
#         if eta - 0.5 < h[i] < eta + 0.5:
#             pdf_Hs[i] = np.mean(pdf_Hs[i - 10:i + 10])
#
#     #####################################################
#     # calcualte a1, a2, a3, b1, b2, b3
#     # firstly calcualte mean_hs, mean_lnTp, variance_lnTp
#     Tp = df.tp.values
#     Hs = df.hs.values
#     maxHs = max(Hs)
#     if maxHs < 2:
#         intx = 0.05
#     elif 2 <= maxHs < 3:
#         intx = 0.1
#     elif 3 <= maxHs < 4:
#         intx = 0.2
#     elif 4 <= maxHs < 10:
#         intx = 0.5
#     else:
#         intx = 1.0
#
#     mean_hs = []
#     variance_lnTp = []
#     mean_lnTp = []
#
#     hs_bin = np.arange(0, maxHs + intx, intx)
#     for i in range(len(hs_bin) - 1):
#         idxs = np.where((hs_bin[i] <= Hs) & (Hs < hs_bin[i + 1]))
#         if Hs[idxs].shape[0] > 15:
#             mean_hs.append(np.mean(Hs[idxs]))
#             mean_lnTp.append(np.mean(np.log(Tp[idxs])))
#             variance_lnTp.append(np.var(np.log(Tp[idxs])))
#
#     mean_hs = np.asarray(mean_hs)
#     mean_lnTp = np.asarray(mean_lnTp)
#     variance_lnTp = np.asarray(variance_lnTp)
#
#     # calcualte a1, a2, a3
#     parameters, covariance = curve_fit(Gauss3, mean_hs, mean_lnTp)
#     a1 = parameters[0]
#     a2 = parameters[1]
#     a3 = 0.36
#
#     # calcualte b1, b2, b3
#     start = 1
#     x = mean_hs[start:]
#     y = variance_lnTp[start:]
#     parameters, covariance = curve_fit(Gauss4, x, y)
#     b1 = 0.005
#     b2 = parameters[0]
#     b3 = parameters[1]
#
#     # calculate pdf Hs, Tp
#     t = np.linspace(start=0.01, stop=40, num=2000)
#
#     f_Hs_Tp = np.zeros((len(h), len(t)))
#     pdf_Hs_Tp = f_Hs_Tp * 0
#
#     for i in range(len(h)):
#         mu = a1 + a2 * h[i] ** a3
#         std2 = b1 + b2 * np.exp(-b3 * h[i])
#         std = np.sqrt(std2)
#
#         f_Hs_Tp[i, :] = 1 / (np.sqrt(2 * np.pi) * std * t) * np.exp(-(np.log(t) - mu) ** 2 / (2 * std2))
#         pdf_Hs_Tp[i, :] = pdf_Hs[i] * f_Hs_Tp[i, :]
#
#     interval = ((df.index[-1] - df.index[0]).days + 1) * 24 / df.shape[0]  # in hours
#
#     t3 = []
#     h3 = []
#     X = []
#     hs_tpl_tph = pd.DataFrame()
#
#     # Assuming Hs_Tp_curve() returns four values, otherwise adjust accordingly
#     for i in range(len(periods)):
#         t3_val, h3_val, X_val, hs_tpl_tph_val = Hs_Tp_curve(df.hs.values, pdf_Hs, pdf_Hs_Tp, f_Hs_Tp, h, t, interval, X=periods_adj[i])
#         t3.append(t3_val)
#         h3.append(h3_val)
#         X.append(X_val)
#         hs_tpl_tph_val.columns = [f'{col}_{periods[i]}' for col in hs_tpl_tph_val.columns]
#         hs_tpl_tph = pd.concat([hs_tpl_tph, hs_tpl_tph_val], axis=1)
#
#     # if save_rve:
#     #    hs_tpl_tph[3].to_csv(str(param[2])+'_year.csv', index=False)
#
#     return a1, a2, a3, b1, b2, b3, pdf_Hs, h, t3, h3, X, hs_tpl_tph
#
#
# def Tp_correction(Tp):
#     """
#     This function will correct the Tp from ocean model which are vertical straight lines in Hs-Tp distribution
#     """
#     new_Tp = 1 + np.log(Tp / 3.244) / 0.09525
#     index = np.where(Tp >= 3.2)  # indexes of Tp
#     r = np.random.uniform(low=-0.5, high=0.5, size=len(Tp[index]))
#     Tp[index] = np.round(3.244 * np.exp(0.09525 * (new_Tp[index] - 1 - r)), 1)
#     return Tp
#
#
# def Hs_Tp_curve(data, pdf_Hs, pdf_Hs_Tp, f_Hs_Tp, h, t, interval, X=100):
#     import scipy.stats as stats
#     from scipy.signal import find_peaks
#
#     # RVE of X years
#     shape, loc, scale = Weibull_method_of_moment(data)  # shape, loc, scale
#
#     if X == 1:
#         period = 1.5873 * 365.2422 * 24 / interval
#     else:
#         period = X * 365.2422 * 24 / interval
#     rve_X = stats.weibull_min.isf(1 / period, shape, loc, scale)
#
#     # Find index of Hs=value
#     epsilon = abs(h - rve_X)
#     param = find_peaks(1 / epsilon)  # to find the index of bottom
#     index = param[0][0]  # the  index of Hs=value
#
#     # Find peak of pdf at Hs=RVE of X year
#     pdf_Hs_Tp_X = pdf_Hs_Tp[index, :]  # Find pdf at RVE of X year
#     param = find_peaks(pdf_Hs_Tp_X)  # find the peak
#     index = param[0][0]
#     f_Hs_Tp_100 = pdf_Hs_Tp_X[index]
#
#     h1 = []
#     t1 = []
#     t2 = []
#     for i in range(len(h)):
#         f3_ = f_Hs_Tp_100 / pdf_Hs[i]
#         f3 = f_Hs_Tp[i, :]
#         epsilon = abs(f3 - f3_)  # the difference
#         para = find_peaks(1 / epsilon)  # to find the bottom
#         index = para[0]
#         if t[index].shape[0] == 2:
#             h1.append(h[i])
#             t1.append(t[index][0])
#             t2.append(t[index][1])
#
#     h1 = np.asarray(h1)
#     t1 = np.asarray(t1)
#     t2 = np.asarray(t2)
#     t3 = np.concatenate((t1, t2[::-1]))  # to get correct circle order
#     h3 = np.concatenate((h1, h1[::-1]))  # to get correct circle order
#     t3 = np.concatenate((t3, t1[0:1]))  # connect the last to the first point
#     h3 = np.concatenate((h3, h1[0:1]))  # connect the last to the first point
#
#     df = pd.DataFrame()
#     df['hs'] = h1
#     df['t1'] = t1
#     df['t2'] = t2
#
#     return t3, h3, X, df
#
#
# def Gauss3(x, a1, a2):
#     y = a1 + a2 * x ** 0.36
#     return y
#
#
# def Gauss4(x, b2, b3):
#     y = 0.005 + b2 * np.exp(-x * b3)
#     return y
#
#
# import matplotlib.pyplot as plt
#
#
# def plot_joint_distribution_Hs_Tp(data, var_hs='hs', var_tp='tp', periods=None, title='Hs-Tp joint distribution', output_file='Hs.Tp.joint.distribution.png',
#                                   density_plot=False):
#     if periods is None:
#         periods = [1, 10, 100, 10000]
#     a1, a2, a3, b1, b2, b3, pdf_Hs, h, t3, h3, X, hs_tpl_tph = joint_distribution_Hs_Tp(data=data, var_hs=var_hs, var_tp=var_tp, periods=periods)
#     df = data
#     # calculate pdf Hs, Tp
#     t = np.linspace(start=0.01, stop=40, num=2000)
#
#     f_Hs_Tp = np.zeros((len(h), len(t)))
#     pdf_Hs_Tp = f_Hs_Tp * 0
#
#     for i in range(len(h)):
#         mu = a1 + a2 * h[i] ** a3
#         std2 = b1 + b2 * np.exp(-b3 * h[i])
#         std = np.sqrt(std2)
#
#         f_Hs_Tp[i, :] = 1 / (np.sqrt(2 * np.pi) * std * t) * np.exp(-(np.log(t) - mu) ** 2 / (2 * std2))
#         pdf_Hs_Tp[i, :] = pdf_Hs[i] * f_Hs_Tp[i, :]
#
#     interval = ((df.index[-1] - df.index[0]).days + 1) * 24 / df.shape[0]  # in hours
#     t_steepness, h_steepness = DVN_steepness(df, h, t, periods, interval)
#     percentile05 = find_percentile(df.hs.values, pdf_Hs_Tp, h, t, 5, periods, interval)
#     percentile50 = find_percentile(df.hs.values, pdf_Hs_Tp, h, t, 50, periods, interval)
#     percentile95 = find_percentile(df.hs.values, pdf_Hs_Tp, h, t, 95, periods, interval)
#
#     fig, ax = plt.subplots(figsize=(8, 6))
#     df = df[df['hs'] >= 0.1]
#     if density_plot is False:
#         plt.scatter(df.tp.values, df.hs.values, c='red', label='data', s=3)
#     else:
#         from matplotlib.colors import LogNorm
#         # plt.scatter(df.tp.values,df.hs.values,c='red',label='data',s=3)
#         plt.hist2d(df['tp'].values, df['hs'].values, bins=50, cmap='hot', cmin=1)
#         plt.colorbar()
#
#     for i in range(len(periods)):
#         plt.plot(t3[i], h3[i], label=str(X[i]) + '-year')
#
#     plt.plot(t_steepness, h_steepness, 'k--', label='steepness')
#
#     plt.plot(percentile50[0], percentile50[1], 'g', label='Tp-mean', linewidth=5)
#     plt.plot(percentile05[0], percentile05[1], 'g:', label='Tp-5%', linewidth=2)
#     plt.plot(percentile95[0], percentile95[1], 'g--', label='Tp-95%', linewidth=2)
#
#     plt.xlabel('Tp - Peak Period [s]')
#     plt.suptitle(title)
#     plt.ylabel('Hs - Significant Wave Height [m]')
#     plt.grid()
#     plt.legend()
#     plt.xlim([0, np.round(hs_tpl_tph['t2_' + str(np.max(periods))].max() + 1)])
#     plt.ylim([0, np.round(hs_tpl_tph['hs_' + str(np.max(periods))].max() + 1)])
#     plt.savefig(output_file, dpi=100, facecolor='white', bbox_inches='tight')
#
#     return fig
#
#
# def DVN_steepness(df, h, t, periods, interval):
#     import scipy.stats as stats
#     ## steepness
#     max_y = max(periods)
#     X = max_y  # get max 500 year
#     period = X * 365.2422 * 24 / interval
#     shape, loc, scale = Weibull_method_of_moment(df.hs.values)  # shape, loc, scale
#     rve_X = stats.weibull_min.isf(1 / period, shape, loc, scale)
#
#     h1 = []
#     t1 = []
#     h2 = []
#     t2 = []
#     h3 = []
#     t3 = []
#     g = 9.80665
#     j15 = 10000
#     for j in range(len(t)):
#         if t[j] <= 8:
#             Sp = 1 / 15
#             temp = Sp * g * t[j] ** 2 / (2 * np.pi)
#             if temp <= rve_X:
#                 h1.append(temp)
#                 t1.append(t[j])
#
#             j8 = j  # t=8
#             h1_t8 = temp
#             t8 = t[j]
#         elif t[j] >= 15:
#             Sp = 1 / 25
#             temp = Sp * g * t[j] ** 2 / (2 * np.pi)
#             if temp <= rve_X:
#                 h3.append(temp)
#                 t3.append(t[j])
#             if j < j15:
#                 j15 = j  # t=15
#                 h3_t15 = temp
#                 t15 = t[j]
#
#     xp = [t8, t15]
#     fp = [h1_t8, h3_t15]
#     t2_ = t[j8 + 1:j15]
#     h2_ = np.interp(t2_, xp, fp)
#     for i in range(len(h2_)):
#         if h2_[i] <= rve_X:
#             h2.append(h2_[i])
#             t2.append(t2_[i])
#
#     h_steepness = np.asarray(h1 + h2 + h3)
#     t_steepness = np.asarray(t1 + t2 + t3)
#
#     return t_steepness, h_steepness
#
#
# def find_percentile(data, pdf_Hs_Tp, h, t, p, periods, interval):
#     import scipy.stats as stats
#     from scipy.signal import find_peaks
#
#     ## find pecentile
#     # RVE of X years
#     max_y = max(periods)
#     X = max_y  # get max 500 year
#     period = X * 365.2422 * 24 / interval
#     shape, loc, scale = Weibull_method_of_moment(data)  # shape, loc, scale
#     rve_X = stats.weibull_min.isf(1 / period, shape, loc, scale)
#     epsilon = abs(h - rve_X)
#     param = find_peaks(1 / epsilon)  # to find the index of bottom
#     index_X = param[0][0]  # the  index of Hs=value
#
#     h1 = []
#     t1 = []
#     # Find peak of pdf at Hs=RVE of X year
#     for i in range(index_X):
#         pdf_Hs_Tp_X = pdf_Hs_Tp[i, :]  # Find pdf at RVE of X year
#         sum_pdf = sum(pdf_Hs_Tp_X)
#         for j in range(len(pdf_Hs_Tp_X)):
#             if (sum(pdf_Hs_Tp_X[:j]) / sum_pdf <= p / 100) and (sum(pdf_Hs_Tp_X[:j + 1]) / sum_pdf >= p / 100):
#                 # print (i, h[i],j,t[j])
#                 t1.append(t[j])
#                 h1.append(h[i])
#                 break
#     h1 = np.asarray(h1)
#     t1 = np.asarray(t1)
#
#     return t1, h1
#
#
# df = gl.export_df_from_sql(db_path, table_name='Hindcast_combined', column_names=[COLNAMES["H_s"], COLNAMES["T_p"]])
# out = joint_distribution_Hs_Tp(df, var_hs=COLNAMES["H_s"], var_tp=COLNAMES["T_p"])
#
# fig = plot_joint_distribution_Hs_Tp(df, var_hs=COLNAMES["H_s"], var_tp=COLNAMES["T_p"], density_plot=True)
#
# print("1")