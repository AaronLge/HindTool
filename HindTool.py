import argparse
import datetime
import inspect
import os
import shutil
from warnings import simplefilter
import numpy as np
import pandas as pd
import scipy as sc
import warnings
import re
import subprocess
import webbrowser

from matplotlib.colors import LinearSegmentedColormap

from libaries import general as gl
from libaries import hindtoolcalc as hc_calc
from libaries import hindtoolplot as hc_plt
from libaries import latex as ltx

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
def main():
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


    if args.o is None:
        if INPUT['DataOut']['dir_name'] is None:
            path_out = os.path.abspath(INPUT['DataOut']['path_out']) + '\\HindCast_' + timestamp + '/'
        else:
            path_out = INPUT['SELECT']['path_out'] + \
                       INPUT['SELECT']['dir_name'] + '/'
    else:
        path_out = os.path.abspath(args.o) + '/'

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
        table_name = 'Hind_combined'
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

    if (INPUT["Toggle_Modules"].get("plot_condensation_example", {})):
        print("calculating VMHS Wind Sea example plot...")

        Input = INPUT["VMHS_docu"]
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

        DATA_OUT["VMHS"]["wind_example"] = Calc

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
            table_name = 'Hind_combined'
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

                Line_mean = hc_plt.Line(x=Seg.result["data"]["x"],
                                        y=Seg.result["data"]["mean result plot"],
                                        label='Selected correlation',
                                        color='black')

                scatter = hc_plt.Scatter(x=point_data[Seg.colnames["x"]],
                                         y=point_data[Seg.colnames["y"]],
                                         cmap='cool',
                                         size=2,
                                         cmap_norm='sqrt')

                tile_curr.add_scatter(scatter)
                tile_curr.add_line(Line_mean)

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

    if INPUT["Toggle_Modules"].get("plot_condensation_example", {}):

        Calc = DATA_OUT["VMHS"]["wind_example"]

        Tiles = []
        Tiles_omni = []

        df = Calc.load_from_db(colnames_ini=True)
        titels = Calc.create_segment_title()
        titels = gl.alias(titels, COLNAMES, INPUT["Aliase"])

        for i, Seg in enumerate(Calc.result):

            if Seg.angles is None:
                Seg.indizes = pd.to_datetime(Seg.indizes)
                point_data = df[df.index.isin(Seg.indizes)]
                reg_zone = np.where(Seg.result["data"]["bool_reg_zone"] == 1)[0]
                use_reg_zone = np.where(Seg.result["data"]["use_regression"] == 1)[0]
                use_mean = np.where(Seg.result["data"]["use_regression"] == 0)[0]

                tile_curr = hc_plt.Tile(i,
                                        x_label=gl.alias(Seg.colnames['x'], COLNAMES, INPUT["Aliase"]),
                                        y_label=gl.alias(Seg.colnames['y'], COLNAMES, INPUT["Aliase"]),
                                        title=titels[i])

                scatter = hc_plt.Scatter(x=point_data[Seg.colnames["x"]],
                                         y=point_data[Seg.colnames["y"]],
                                         size=5,
                                         cmap_norm='sqrt',
                                         color=[0.8, 0.8, 0.8])

                Line_mean = hc_plt.Line(x=Seg.result["data"]["x"],
                                        y=Seg.result["data"]["mean"],
                                        color='black',
                                        linewidth=2.5,
                                        alpha=0.5)

                Line_mean_inrange = hc_plt.Line(x=Seg.result["data"]["x"].iloc[use_mean[0]:use_mean[-1] + 2],
                                                y=Seg.result["data"]["mean"].iloc[use_mean[0]:use_mean[-1] + 2],
                                                color='black',
                                                linewidth=2.5)

                scatter_mean = hc_plt.Scatter(x=Seg.result["data"]["x"],
                                              y=Seg.result["data"]["mean"],
                                              label='$y$ (representative bin values)',
                                              color='black',
                                              size=30)

                Line_regression = hc_plt.Line(x=Seg.result["data"]["x"],
                                              y=Seg.result["data"]["mean regression"],
                                              color='blue',
                                              linewidth=2.5,
                                              alpha=0.3)

                Line_regression_inrange = hc_plt.Line(x=Seg.result["data"]["x"].iloc[use_reg_zone],
                                                      y=Seg.result["data"]["mean regression"].iloc[use_reg_zone],
                                                      label='Resulting regression curve',
                                                      color='blue',
                                                      linewidth=1.5)

                Line_reg_zone_left = hc_plt.Line(x=[Seg.result["data"]['x'].iloc[reg_zone[0]]],
                                                 y=None,
                                                 label=r'Range of the regression base (start, end)',
                                                 color='#20D503',
                                                 linewidth=2)

                Line_reg_zone_right = hc_plt.Line(x=[Seg.result["data"]['x'].iloc[reg_zone[-1]]],
                                                  y=None,
                                                  color='#20D503',
                                                  linewidth=2)

                Line_use_reg = hc_plt.Line(x=[Seg.result["data"]['x'].iloc[use_reg_zone[0]]],
                                           y=None,
                                           label=r'Start of datapoints evaluated with regression curve',
                                           color='#20D503',
                                           linewidth=2,
                                           linestyle='--')

                scatter_regression = hc_plt.Scatter(x=Seg.result["data"]["x"].iloc[reg_zone],
                                                    y=Seg.result["data"]["mean"].iloc[reg_zone],
                                                    label=r'Values used for regression (inside "Range of the regression base")',
                                                    color='blue',
                                                    size=40,
                                                    marker='x')

                Line_result = hc_plt.Line(x=Seg.result["data"]["x"],
                                          y=Seg.result["data"]["mean result"],
                                          color='red',
                                          linestyle='--',
                                          linewidth=1)

                Scatter_result = hc_plt.Scatter(x=Seg.result["data"]["x"],
                                                y=Seg.result["data"]["mean result"],
                                                label=r'Selected correlation',
                                                color='red',
                                                size=5)

                tile_curr.add_scatter(scatter)
                # tile_curr.add_line(Line_mean)
                tile_curr.add_line(Line_regression)

                tile_curr.add_line(Line_regression_inrange)
                #  tile_curr.add_line(Line_mean_inrange)

                tile_curr.add_scatter(scatter_mean)
                tile_curr.add_scatter(scatter_regression)

                tile_curr.add_line(Line_reg_zone_right)
                tile_curr.add_line(Line_reg_zone_left)
                tile_curr.add_line(Line_use_reg)
                tile_curr.add_scatter(Scatter_result)

                Tiles.append(tile_curr)

                Tiles_omni.append(tile_curr)

        FIG_omni = hc_plt.plot_tiled(Tiles_omni, global_max=['auto', 'auto'], global_min=[0, 0], grid=[1, 1], figsize=figsize_halfpage)

        if 'png' in INPUT["Toggle_Modules"]["plot_as"]:
            gl.save_figs_as_png(FIG_omni, path_out + 'VMHS_example', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

        if 'pdf' in INPUT["Toggle_Modules"]["plot_as"]:
            gl.save_figs_as_pdf(FIG_omni, path_out + 'VMHS_example', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

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

                Line_mean = hc_plt.Line(x=Seg.result["data"]["x"],
                                        y=Seg.result["data"]["mean result plot"],
                                        label='Selected correlation',
                                        color='black')

                scatter = hc_plt.Scatter(x=point_data[Seg.colnames["x"]],
                                         y=point_data[Seg.colnames["y"]],
                                         cmap='cool',
                                         size=2,
                                         cmap_norm='sqrt')

                tile_curr.add_scatter(scatter)
                tile_curr.add_line(Line_mean)

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

                Line_mean = hc_plt.Line(x=Seg.result["data"]["x"],
                                        y=Seg.result["data"]["mean result plot"],
                                        label='Selected correlation',
                                        color='black')

                scatter = hc_plt.Scatter(x=point_data[Seg.colnames["x"]],
                                         y=point_data[Seg.colnames["y"]],
                                         cmap='cool',
                                         size=2,
                                         cmap_norm='sqrt')

                tile_curr.add_scatter(scatter)
                tile_curr.add_line(Line_mean)

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

                key_plot = [key for key in Seg.result["data"].columns if 'plot' in key]
                key_mean = [key for key in key_plot if 'mean' in key]
                key_perc = [key for key in key_plot if 'percentile' in key]
                key_quantile = [key for key in Seg.result["data"].columns if 'quantile' in key]

                Line_perc_low = hc_plt.Line(x=Seg.result["data"]["x"],
                                            y=Seg.result["data"][key_perc[0]],
                                            label=key_perc[0].replace('result', '').replace('plot', ''),
                                            color='black',
                                            linestyle='--')

                Line_mean = hc_plt.Line(x=Seg.result["data"]["x"],
                                        y=Seg.result["data"]["mean result plot"],
                                        label='50 percentile',
                                        color='black',
                                        linestyle='-')

                Line_perc_up = hc_plt.Line(x=Seg.result["data"]["x"],
                                           y=Seg.result["data"][key_perc[1]],
                                           label=key_perc[1].replace('result', '').replace('plot', ''),
                                           color='black',
                                           linestyle=':')

                scatter = hc_plt.Scatter(x=point_data[Seg.colnames["x"]],
                                         y=point_data[Seg.colnames["y"]],
                                         cmap='cool',
                                         size=2,
                                         cmap_norm='sqrt')

                tile_curr.add_scatter(scatter)
                tile_curr.add_line(Line_perc_low)
                tile_curr.add_line(Line_perc_up)
                tile_curr.add_line(Line_mean)

                if len(key_quantile) > 0:
                    Line_quant = hc_plt.Line(x=Seg.result["data"]["x"],
                                             y=Seg.result["data"][key_quantile[0]],
                                             label='Selected correlation',
                                             color='red',
                                             linestyle='-')
                    tile_curr.add_line(Line_quant)

                    Line_quant_up = hc_plt.Line(x=None,
                                                y=[1 / Input["quant_up"]],
                                                label=r'$f_{up}$' + f"$={round(Input['quant_up'], 3)}$ Hz",
                                                color='green',
                                                linestyle=':')
                    tile_curr.add_line(Line_quant_up)

                    Line_quant_low = hc_plt.Line(x=None,
                                                 y=[1 / Input["quant_low"]],
                                                 label=r'$f_{low}$' + f"$={round(Input['quant_low'], 3)}$ Hz",
                                                 color='green',
                                                 linestyle='--')
                    tile_curr.add_line(Line_quant_low)

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

                scatter = hc_plt.Scatter(x=point_data[Seg.colnames["x"]],
                                         y=point_data[Seg.colnames["y"]],
                                         cmap='cool',
                                         size=2,
                                         cmap_norm='sqrt')

                key_plot = [key for key in Seg.result["data"].columns if 'plot' in key]
                key_mean = [key for key in key_plot if 'mean' in key]
                key_perc = [key for key in key_plot if 'percentile' in key]
                key_quantile = [key for key in Seg.result["data"].columns if 'quantile' in key]

                Line_mean = hc_plt.Line(x=Seg.result["data"]["x"],
                                        y=Seg.result["data"]["mean result plot"],
                                        label='Selected correlation',
                                        color='black',
                                        linestyle='-')

                tile_curr.add_scatter(scatter)
                tile_curr.add_line(Line_mean)

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

                key_plot = [key for key in Seg.result["data"].columns if 'plot' in key]
                key_mean = [key for key in key_plot if 'mean' in key]
                key_perc = [key for key in key_plot if 'percentile' in key]
                key_quantile = [key for key in Seg.result["data"].columns if 'quantile' in key]

                Line_mean = hc_plt.Line(x=Seg.result["data"]["x"],
                                        y=Seg.result["data"]["mean result plot"],
                                        label='Selected correlation',
                                        color='black',
                                        linestyle='-')

                scatter = hc_plt.Scatter(x=point_data[Seg.colnames["x"]],
                                         y=point_data[Seg.colnames["y"]],
                                         cmap='cool',
                                         size=2,
                                         cmap_norm='sqrt')

                tile_curr.add_scatter(scatter)
                tile_curr.add_line(Line_mean)

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

                Line_mean = hc_plt.Line(x=Seg.result["data"]["x"],
                                        y=Seg.result["data"]["mean result plot"],
                                        label='Extracted correlation',
                                        color='black')

                scatter = hc_plt.Scatter(x=point_data[Seg.colnames["x"]],
                                         y=point_data[Seg.colnames["y"]],
                                         cmap='cool',
                                         size=2,
                                         cmap_norm='sqrt')

                tile_curr.add_scatter(scatter)
                tile_curr.add_line(Line_mean)

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

                Line_mean = hc_plt.Line(x=Seg.result["data"]["x"],
                                        y=Seg.result["data"]["mean result plot"],
                                        label='Extracted correlation',
                                        color='black')

                scatter = hc_plt.Scatter(x=point_data[Seg.colnames["x"]],
                                         y=point_data[Seg.colnames["y"]],
                                         cmap='cool',
                                         size=2,
                                         cmap_norm='sqrt')

                tile_curr.add_scatter(scatter)
                tile_curr.add_line(Line_mean)

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

                Line_mean = hc_plt.Line(x=Seg.result["data"]["x"],
                                        y=Seg.result["data"]["mean result plot"],
                                        label='Extracted correlation',
                                        color='black')

                scatter = hc_plt.Scatter(x=point_data[Seg.colnames["x"]],
                                         y=point_data[Seg.colnames["y"]],
                                         cmap='cool',
                                         size=2,
                                         cmap_norm='sqrt')

                tile_curr.add_scatter(scatter)
                tile_curr.add_line(Line_mean)

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

                Line_f0 = hc_plt.Line(x=None,
                                      y=[1 / INPUT["Structure"]["f_0"]],
                                      label=f'f_0, f={INPUT["Structure"]["f_0"]} Hz',
                                      color='green',
                                      linestyle=':')

                tile_curr.add_scatter(scatter)
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

                Line_f0 = hc_plt.Line(x=None,
                                      y=[1 / INPUT["Structure"]["f_0"]],
                                      label=f'$f_0={INPUT["Structure"]["f_0"]}$ Hz',
                                      color='green',
                                      linestyle=':')

                tile_curr.add_scatter(scatter)
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
                                       label=f"Rolling absolute mean (global mean = {round(np.mean(Seg.result['mean'].values), 2)} deg)",
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
                scatter = hc_plt.Scatter(x=Seg.result["points"]["x_max"].index,
                                         y=Seg.result["points"]["x_max"].values,
                                         color='red',
                                         label='Extreme Values',
                                         size=10)

                tile_curr.add_line(Line)
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
                                         size=10,
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
                                         size=10)

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

                    Line_hindcast = hc_plt.Line(x=Seg.result["hindcast"]['vm_vise'][config].index,
                                                y=Seg.result["hindcast"]['vm_vise'][config].values,
                                                color=cmap_lines(range_colors[i]),
                                                linestyle='-')

                    tile_curr.add_line(Line_hindcast, zorder=10)
                    tile_curr.add_line(Line_condensed, zorder=10)

                    # textbox
                    textbox.append(str(config))
                    textbox.append(f"DEL Hindcast: {Seg.result['hindcast']['added'][config].values[0]:.3e}")
                    textbox.append(f"DEL Condensed: {Seg.result['condensed']['added'][config].values[0]:.3e}")
                    textbox.append(
                        f"Condensed/Hindcast: {round(Seg.result['condensed']['added'][config].values[0] / Seg.result['hindcast']['added'][config].values[0] * 100, 1)}" + "$\\%$")
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

                    tile_curr.add_textbox(Textbox_DEL, zorder=9)
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
                                          label=f'$f_0={INPUT["Structure"]["f_0"]}$ Hz',
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
                    tile_curr.add_line(Line_condensed, zorder=10)

                    Line_hindcast = hc_plt.Line(x=Seg.result["hindcast"]['vm_vise'][config].index,
                                                y=Seg.result["hindcast"]['vm_vise'][config].values,
                                                color=cmap_lines(range_colors[i]),
                                                linestyle='-')

                    tile_curr.add_line(Line_hindcast, zorder=10)

                    # textbox
                    textbox.append(str(config))
                    textbox.append(f"DEL Hindcast: {Seg.result['hindcast']['added'][config].values[0]:.3e}")
                    textbox.append(f"DEL Condensed: {Seg.result['condensed']['added'][config].values[0]:.3e}")
                    textbox.append(
                        f"Condensed/Hindcast: {round(Seg.result['condensed']['added'][config].values[0] / Seg.result['hindcast']['added'][config].values[0] * 100, 1)}" + "$\\%$")
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

                tile_curr.add_bar(Bar_count, zorder=0)

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

                    tile_curr.add_textbox(Textbox_DEL, zorder=-1)
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
                                          label=f'$f_0={INPUT["Structure"]["f_0"]}$ Hz',
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


    # %% plot report tables
    if INPUT["Toggle_Modules"].get("plot_report_tables", {}) or INPUT["DataBase"].get("create_report", {}):
        try:
            INPUT_REPORT = gl.read_input_txt(INPUT["DataBase"]["Report_Input"])
        except:
            print("please specify Report_Input")
            return

        path_report = os.path.join(path_out, 'report')
        try:
            # Create the new folder
            os.makedirs(path_report, exist_ok=True)  # exist_ok=True prevents an error if the folder already exists
            print(f"Folder created successfully at: {path_report}")
        except Exception as e:
            print(f"An error occurred: {e}")

        # crate COLNAME dataframe with symbols as master
        COLNAMES_REPORT = pd.DataFrame(index=INPUT_REPORT["Symbols"].keys())
        COLNAMES_REPORT["Symbols"] = INPUT_REPORT["Symbols"].values()
        COLNAMES_REPORT["Sensor_names"] = [COLNAMES[key] if key in COLNAMES else float('nan') for key in list(COLNAMES_REPORT.index)]
        COLNAMES_REPORT["Aliase"] = [INPUT["Aliase"][key] if key in INPUT["Aliase"] else float('nan') for key in COLNAMES_REPORT.index]

        if INPUT["Toggle_Modules"].get("plot_report_tables", {}):

            # plot databases
            Meta_data = gl.export_df_from_sql(db_path, 'Hind_MetaData')

            for dataset_name, dataset_contents in Meta_data.iterrows():
                meta_para = []
                meta_value = []
                for dataset_para, dataset_value in dataset_contents.items():
                    if dataset_value is not None:
                        meta_para.append(dataset_para)
                        meta_value.append(dataset_value)

                # plot Plot_names
                data = np.array([meta_para, meta_value])
                data = data.T
                col_labels = ["Parameter", "Value"]
                FIG = hc_plt.table(data,
                                   collabels=col_labels,
                                   figsize=figsize_fullpage,
                                   datatype='str')

                if 'png' in INPUT["Toggle_Modules"]["plot_as"]:
                    gl.save_figs_as_png([FIG], path_out + f'DataSorce_{dataset_name}', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

                if 'pdf' in INPUT["Toggle_Modules"]["plot_as"]:
                    gl.save_figs_as_pdf([FIG], path_out + f'DataSorce_{dataset_name}', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

            # add datasorce reference
            datasorce_cols = gl.export_colnames_from_db(db_path)
            datasorce_keys_raw = [col for col in datasorce_cols.keys() if "Hind_raw" in col]

            for sensor_key, sensor_name in COLNAMES_REPORT["Sensor_names"].items():
                for datasorce_key_raw in datasorce_keys_raw:
                    if sensor_name in datasorce_cols[datasorce_key_raw]:
                        datasorce_key_clean = datasorce_key_raw.replace('Hind_raw_', '').replace('.csv', '')
                        COLNAMES_REPORT.loc[sensor_key, "DataSorce"] = datasorce_key_clean

            # plot sensor names
            data = np.array([COLNAMES_REPORT["Symbols"], COLNAMES_REPORT["Sensor_names"], COLNAMES_REPORT["DataSorce"]])
            data = data.T
            col_labels = ["Symbol", "Sensor name", "Data Sorce"]
            FIG = hc_plt.table(data,
                               collabels=col_labels,
                               figsize=figsize_fullpage,
                               datatype='str')

            if 'png' in INPUT["Toggle_Modules"]["plot_as"]:
                gl.save_figs_as_png([FIG], path_out + 'Sensor_names', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

            if 'pdf' in INPUT["Toggle_Modules"]["plot_as"]:
                gl.save_figs_as_pdf([FIG], path_out + 'Sensor_names', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

            # plot Plot_names
            data = np.array([COLNAMES_REPORT["Symbols"], COLNAMES_REPORT["Aliase"]])
            data = data.T
            col_labels = ["Symbol", "Plot name"]
            FIG = hc_plt.table(data,
                               collabels=col_labels,
                               figsize=figsize_fullpage,
                               datatype='str')

            if 'png' in INPUT["Toggle_Modules"]["plot_as"]:
                gl.save_figs_as_png([FIG], path_out + 'Plot_names', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

            if 'pdf' in INPUT["Toggle_Modules"]["plot_as"]:
                gl.save_figs_as_pdf([FIG], path_out + 'Plot_names', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

            # VMHS parameter
            if INPUT["Toggle_Modules"].get("plot_VMHS", {}):
                print("   plotting VMHS parameter table")

                columns_table = []
                col_labels = []
                row_labels = ['Degree of regression',
                              'Shape function f(x)',
                              'Bin number',
                              'Portion of datapoints evaluated with mean',
                              'Regression range',
                              'Evaluated range',
                              'Applied percentile for condensation',
                              'Averageing method']

                if 'wind' in INPUT["Toggle_Modules"].get("plot_VMHS", {}):

                    Input = INPUT["VMHS_wind"]

                    new_col = [
                        Input["deg_reg"],
                        'x' if Input["model_reg"] == 'poly' else r'$\sqrt{x}$' if Input["model_reg"] == 'sqrt' else '',
                        Input["N_grid"],
                        f"{Input['cut_reg']}" + r" \%",
                        f"[{Input['zone_reg'][0]} .. {'max' if Input['zone_reg'][1] is None else Input['zone_reg'][1]}]",
                        f"[{Input['zone_line'][0]} .. {'max' if Input['zone_line'][1] is None else Input['zone_line'][1]}]",
                        f"{Input['perc_mean']}" + r" \%",
                        Input["avrg_method"],
                               ]
                    columns_table.append(new_col)
                    col_labels.append('Wind Sea')

                if 'swell' in INPUT["Toggle_Modules"].get("plot_VMHS", {}):

                    Input = INPUT["VMHS_swell"]

                    new_col = [
                        Input["deg_reg"],
                        'x' if Input["model_reg"] == 'poly' else r'$\sqrt{x}$' if Input["model_reg"] == 'sqrt' else '',
                        Input["N_grid"],
                        f"{Input['cut_reg']}" + r" \%",
                        f"[{Input['zone_reg'][0]} .. {'max' if Input['zone_reg'][1] is None else Input['zone_reg'][1]}]",
                        f"[{Input['zone_line'][0]} .. {'max' if Input['zone_line'][1] is None else Input['zone_line'][1]}]",
                        f"{Input['perc_mean']}" + r" \%",
                        Input["avrg_method"],
                               ]
                    columns_table.append(new_col)
                    col_labels.append('Swell Sea')

                if 'total' in INPUT["Toggle_Modules"].get("plot_VMHS", {}):

                    Input = INPUT["VMHS_total"]

                    new_col = [
                        Input["deg_reg"],
                        'x' if Input["model_reg"] == 'poly' else r'$\sqrt{x}$' if Input["model_reg"] == 'sqrt' else '',
                        Input["N_grid"],
                        f"{Input['cut_reg']}" + r" %",
                        f"[{Input['zone_reg'][0]} .. {'max' if Input['zone_reg'][1] is None else Input['zone_reg'][1]}]",
                        f"[{Input['zone_line'][0]} .. {'max' if Input['zone_line'][1] is None else Input['zone_line'][1]}]",
                        f"{Input['perc_mean']}" + r" %",
                        Input["avrg_method"],
                               ]
                    columns_table.append(new_col)
                    col_labels.append('Total Sea')

                if len(columns_table) > 0:
                    data = np.array(columns_table)
                    data = data.T

                    FIG = hc_plt.table(data, collabels=col_labels, rowlabels=row_labels, row_label_name='Parameters', figsize=figsize_halfpage)

                    if 'png' in INPUT["Toggle_Modules"]["plot_as"]:
                        gl.save_figs_as_png([FIG], path_out + 'Report_table_VMHS', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

                    if 'pdf' in INPUT["Toggle_Modules"]["plot_as"]:
                        gl.save_figs_as_pdf([FIG], path_out + 'Report_table_VMHS', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

            # HSTP parameter
            if INPUT["Toggle_Modules"].get("plot_HSTP", {}):
                print("   plotting HSTP parameter table")

                columns_table = []
                col_labels = []
                row_labels = ['Degree of regression',
                              'Shape function f(x)',
                              'Bin number',
                              'Portion of datapoints evaluated with mean',
                              'Regression range',
                              'Evaluated range',
                              'Quantiles for selected correlation',
                              'Frequency band for selected correlation',
                              'Applied percentile for condensation',
                              'Averageing method']

                Input = INPUT["HSTP_wind"]

                if 'wind' in INPUT["Toggle_Modules"].get("plot_HSTP", {}):
                    Input = INPUT["HSTP_wind"]

                    new_col = [
                        Input["deg_reg"],
                        'x' if Input["model_reg"] == 'poly' else r'$\sqrt{x}$' if Input["model_reg"] == 'sqrt' else '',
                        Input["N_grid"],
                        f"{Input['cut_reg']}" + r" $\%$",
                        f"[{Input['zone_reg'][0]} .. {'max' if Input['zone_reg'][1] is None else Input['zone_reg'][1]}]",
                        f"[{Input['zone_line'][0]} .. {'max' if Input['zone_line'][1] is None else Input['zone_line'][1]}]",
                        f"[{Input['percentiles'][0]}$\%$ .. {Input['percentiles'][1]}$\%$]" if Input['quantile'] else "none",
                        f"[{Input['quant_up']} .. {Input['quant_low']}]" if Input['quantile'] else "none",
                        f"{Input['perc_mean']}" + r" $\%$",
                        Input["avrg_method"],
                    ]
                    columns_table.append(new_col)
                    col_labels.append('Wind Sea')

                if 'swell' in INPUT["Toggle_Modules"].get("plot_HSTP", {}):
                    Input = INPUT["HSTP_swell"]

                    new_col = [
                        Input["deg_reg"],
                        'x' if Input["model_reg"] == 'poly' else r'$\sqrt{x}$' if Input["model_reg"] == 'sqrt' else '',
                        Input["N_grid"],
                        f"{Input['cut_reg']}" + r" $\%$",
                        f"[{Input['zone_reg'][0]} .. {'max' if Input['zone_reg'][1] is None else Input['zone_reg'][1]}]",
                        f"[{Input['zone_line'][0]} .. {'max' if Input['zone_line'][1] is None else Input['zone_line'][1]}]",
                        f"[{Input['percentiles'][0]}$\%$ .. {Input['percentiles'][1]}$\%$]" if Input['quantile'] else "none",
                        f"[{Input['quant_up']} .. {Input['quant_low']}]" if Input['quantile'] else "none",
                        f"{Input['perc_mean']}" + r" $\%$",
                        Input["avrg_method"],
                    ]
                    columns_table.append(new_col)
                    col_labels.append('Swell Sea')

                if 'total' in INPUT["Toggle_Modules"].get("plot_HSTP", {}):
                    Input = INPUT["HSTP_total"]

                    new_col = [
                        Input["deg_reg"],
                        'x' if Input["model_reg"] == 'poly' else r'$\sqrt{x}$' if Input["model_reg"] == 'sqrt' else '',
                        Input["N_grid"],
                        f"{Input['cut_reg']}" + r" $\%$",
                        f"[{Input['zone_reg'][0]} .. {'max' if Input['zone_reg'][1] is None else Input['zone_reg'][1]}]",
                        f"[{Input['zone_line'][0]} .. {'max' if Input['zone_line'][1] is None else Input['zone_line'][1]}]",
                        f"[{Input['percentiles'][0]}$\%$ .. {Input['percentiles'][1]}$\%$]" if Input['quantile'] else "none",
                        f"[{Input['quant_up']} .. {Input['quant_low']}]" if Input['quantile'] else "none",
                        f"{Input['perc_mean']}" + r" $\%$",
                        Input["avrg_method"],
                    ]
                    columns_table.append(new_col)
                    col_labels.append('Total Sea')

                if len(columns_table) > 0:
                    data = np.array(columns_table)
                    data = data.T

                    FIG = hc_plt.table(data,
                                       collabels=col_labels,
                                       rowlabels=row_labels,
                                       row_label_name='Parameters',
                                       figsize=figsize_halfpage,
                                       datatype='str')

                    if 'png' in INPUT["Toggle_Modules"]["plot_as"]:
                        gl.save_figs_as_png([FIG], path_out + 'Report_table_HSTP', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

                    if 'pdf' in INPUT["Toggle_Modules"]["plot_as"]:
                        gl.save_figs_as_pdf([FIG], path_out + 'Report_table_HSTP', dpi=INPUT["Toggle_Modules"]["dpi_figures"])

        # %% Create Latex File
        if INPUT["DataBase"].get("create_report", {}):

            FIGURES = pd.DataFrame(columns=["filename", "path", "caption" ])

            if INPUT_REPORT["General"]["fig_path"] is not None:
                path_figs = INPUT_REPORT["General"]["fig_path"]
            else:
                path_figs = path_out

            png_files = [f for f in os.listdir(path_figs) if f.endswith('.png')]
            png_paths = [os.path.join(path_figs, f) for f in os.listdir(path_figs) if f.endswith('.png')]
            png_names = [name.removesuffix('.png') for name in png_files]

            FIGURES["filename"] = png_files
            FIGURES["path"] = png_paths

            FIGURES.index = png_names

            # Crete TEX content
            TEX = {}
            TEX_Main = ltx.insertLatexVars(INPUT_REPORT["General"]["path_latex_templates"] + "/template_main.txt", INPUT_REPORT["DocumentMeta"])

            # Titlepage
            TEX["titlepage"] = ltx.insertLatexVars(INPUT_REPORT["General"]["path_latex_templates"] + "/template_titlepage.txt", INPUT_REPORT["DocumentMeta"])
            TEX_Main, last_idx = ltx.include_include(TEX_Main, 'titlepage')

            TEX_Main, last_idx = ltx.include_str(TEX_Main, '\\pagestyle{fancy}', last_idx+1)

            # Introduction
            TEX["introduction"] = ltx.insertLatexVars(INPUT_REPORT["General"]["path_latex_templates"] + "/template_introduction.txt", INPUT_REPORT["introduction"])
            TEX_Main, last_idx = ltx.include_include(TEX_Main, 'introduction', line=last_idx+1)

            #Sensors
            if INPUT["Toggle_Modules"].get("plot_SensorEval", {}) and INPUT["Toggle_Modules"].get("plot_SensorEval", {}):
                print("   including Sensor sheets")

                for sensor_key in INPUT["Toggle_Modules"]["calc_SensorEval"]:

                    sensor_alias = COLNAMES_REPORT["Aliase"][sensor_key]

                    FIGURES.loc[f"SensorEval_{sensor_key}_page_1", "caption"] = f'Timeseries and Hisogram of Sensor: {sensor_alias}'

                    temp = ltx.insertLatexVars(INPUT_REPORT["General"]["path_latex_templates"] + "/template_sensorsheet.txt",
                                                        {'sensorName': sensor_alias} )

                    TEX[f"Sensor_{sensor_key}"] = ltx.insertLatexVars(INPUT_REPORT["General"]["path_latex_templates"] + "/template_sensorsheet.txt",
                                                        {'sensorName': sensor_alias} )

                    TEX[f"Sensor_{sensor_key}"] = ltx.include_Fig(TEX[f"Sensor_{sensor_key}"], FIGURES.loc[f"SensorEval_{sensor_key}_page_1"])
                    TEX_Main, _ = ltx.include_include(TEX_Main, f'Sensor_{sensor_key}')
            else:
                print("   Sensors are not included in Report, please specify Sensors in 'Toggle_Modules' and enable ploting")




            # save TEX files
            with open(path_report + r'\main.tex', 'w', encoding='utf-8') as file:
                file.write(TEX_Main)

            for name, tex in TEX.items():
                with open(path_report+r'\\' + name + '.tex', 'w', encoding='utf-8') as file:
                    file.write(tex)

            # compiling Latex files
            print("   compiling Latex File, please press enter a couple of times...")
            shutil.copy('./latex_templates/JBO_logo.jpg', path_report+r'\\'+'JBO_logo.jpg')

            output_pdf = ltx.compile_lualatex(path_report + r'\main.tex')

            if os.path.exists(output_pdf):
                # Open the PDF file in the default web browser (Edge)
                webbrowser.open_new(output_pdf)
            else:
                print(f"The file {output_pdf} does not exist.")


    # %% Data Out
    if INPUT["DataOut"]["CSV_out"]:
        print("saving CSV Data...")

        path_csv = os.path.join(path_out, 'csv_data')
        Coeffs_string = ("Regression Coefficients \n"
                         "order: [c_0, c_1, ..., c_deg]\n "
                         "with: y(x) = c_0 + c_1 * x + ... + c_deg * x^deg \n \n")
        try:
            # Create the new folder
            os.makedirs(path_csv, exist_ok=True)  # exist_ok=True prevents an error if the folder already exists
            print(f"Folder created successfully at: {path_csv}")
        except Exception as e:
            print(f"An error occurred: {e}")

        unpack_funcs = {"flat_data": lambda x: [seg.result for seg in x.result],
                        "flat_angles": lambda x: [seg.angles for seg in x.result],
                        "flat_exclude_omni": lambda x: [seg.result for seg in x.result if segment.angles is not None],
                        "deep_data": lambda x: [seg.result["data"] for seg in x.result],
                        "deep_coeffs": lambda x: [seg.result["coeffs"] for seg in x.result]}

        for seastate, calc in DATA_OUT["VMHS"].items():
            print(f"   VMHS {seastate}")
            data = unpack_funcs["deep_data"](calc)
            coeffs = unpack_funcs["deep_coeffs"](calc)
            table_names = unpack_funcs["flat_angles"](calc)
            table_names = ["omnidirectional" if name is None else f"{name[0]} to {name[1]}" for name in table_names]

            Coeffs = {table_name: coeff for table_name, coeff in zip(table_names, coeffs)}

            Coeffs_string += f"VMHS {seastate} \n \n"
            Coeffs_string += gl.write_dict(Coeffs) + "\n \n"

            gl.save_df_list_to_excel(path_csv + f'/VMHS_{seastate}', data, sheet_names=table_names)

        for seastate, calc in DATA_OUT["HSTP"].items():
            print(f"   HSTP {seastate}")
            data = unpack_funcs["deep_data"](calc)
            coeffs = unpack_funcs["deep_coeffs"](calc)
            table_names = unpack_funcs["flat_angles"](calc)
            table_names = ["omnidirectional" if name is None else f"{name[0]} to {name[1]}" for name in table_names]

            Coeffs = {table_name: coeff for table_name, coeff in zip(table_names, coeffs)}

            Coeffs_string += f"HSTP {seastate} \n \n"
            Coeffs_string += gl.write_dict(Coeffs) + "\n \n"

            gl.save_df_list_to_excel(path_csv + f'/HSTP_{seastate}', data, sheet_names=table_names)

        for seastate, calc in DATA_OUT["VMTP"].items():
            print(f"   VMTP {seastate}")
            data = unpack_funcs["deep_data"](calc)
            table_names = unpack_funcs["flat_angles"](calc)
            table_names = ["omnidirectional" if name is None else f"{name[0]} to {name[1]}" for name in table_names]

            gl.save_df_list_to_excel(path_csv + f'/VMTP_{seastate}', data, sheet_names=table_names)

        for seastate, calc in DATA_OUT["RWI"].items():
            print(f"   RWI {seastate}")

            df = calc.load_from_db(colnames_ini=True)
            data = unpack_funcs["flat_data"](calc)

            for i, segment in enumerate(calc.result):
                df_temp = df.loc[segment.indizes]

                data[i] = pd.concat((data[i], df_temp), axis=1)

            table_names = unpack_funcs["flat_angles"](calc)
            table_names = ["omnidirectional" if name is None else f"{name[0]} to {name[1]}" for name in table_names]

            gl.save_df_list_to_excel(path_csv + f'/RWI_{seastate}', data, sheet_names=table_names)

        for seastate, calc in DATA_OUT["WaveBreak_Steep"].items():
            print(f"   WaveBreak_Steep {seastate}")

            data = unpack_funcs["flat_data"](calc)
            table_names = unpack_funcs["flat_angles"](calc)
            table_names = ["omnidirectional" if name is None else f"{name[0]} to {name[1]}" for name in table_names]

            gl.save_df_list_to_excel(path_csv + f'/WaveBreak_Steep_{seastate}', data, sheet_names=table_names)

        for seastate, calc in DATA_OUT["table_vmhs"].items():
            print(f"   table_vmhs {seastate}")

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
            gl.save_df_list_to_excel(path_csv + f'/table_vmhs_{seastate}', data, sheet_names=table_names)

        for seastate, calc in DATA_OUT["table_vmtp"].items():
            print(f"   table_vmtp {seastate}")

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
            gl.save_df_list_to_excel(path_csv + f'/table_vmtp_{seastate}', data, sheet_names=table_names)

        if DATA_OUT["AngleDeviation"]:
            print("   AngleDeviation")

            calc = DATA_OUT["AngleDeviation"]
            data = unpack_funcs["flat_exclude_omni"](calc)
            table_names = unpack_funcs["flat_angles"](calc)
            table_names = table_names
            table_names = [f"{name[0]} to {name[1]}" for name in table_names]
            gl.save_df_list_to_excel(path_csv + r'//AngleDeviation', data, sheet_names=table_names)

        for seastate, calc in DATA_OUT["Validation"].items():
            print(f"   Validation {seastate}")

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

            gl.save_df_list_to_excel(path_csv + f'/Validation_{seastate}_hindcast_vm_vise', data_hindcast_vm_vise, sheet_names=table_names)
            gl.save_df_list_to_excel(path_csv + f'/Validation_{seastate}_condensed_vm_vise', data_condensed_vm_vise, sheet_names=table_names)
            gl.save_df_list_to_excel(path_csv + f'/Validation_{seastate}_added', table_wind_added, sheet_names=["hindcast", "condensed"])

        if len(Coeffs_string) > 0:
            with open(path_csv + "/Coeffs.txt", "w") as text_file:
                text_file.write(Coeffs_string)

    # %% MAIN - Save Infolog
    print("saving Log...")

    with open(path_out + 'Info_LOG.txt', "w") as log_text:
        log_text.write(INFO_LOG)
    del log_text

    print(f"{script_name} finished!")

    return

main()