"""
Do postprocessing of all HySEA simulations:
- postprocess_HySea is C++ code that extract info from time series
- ampfact_func.py is a script that uses the amplification factors to compute the inundation
"""

import os
import sys
import pandas as pd
from ampfact_func import ampfact_func

combo_id = ["1798777","6270947","10237389"]
combo_id_str = combo_id[0] + '_' + combo_id[1] + '_' + combo_id[2]

parent_dir = "/home/vmg/NGI/P/2022/02/20220296/Calculations/Messina_simulations_June2025//"+combo_id_str+"//"
hysea_ts=os.path.join(parent_dir,"filterkajiura_res100_ts_"+combo_id_str+".nc")
# File with amplification factors
ampfact_file="/home/vmg/NGI/P/2022/02/20220296/Calculations/testing-Messina-POIs-AmpFactor/landslide_tsunami_HySEA/factors/ampf_messina_new.txt"

# Run first step of the postprocessing that analyses the HySEA time series
command = f"./postprocess_HySea {hysea_ts}"
os.system(command)
# Run the second step of the postprocessing that computes the inundation
outfile_name = os.path.join(parent_dir, 'MIH.txt')
ampfact_func(ampfact_file, hysea_ts, outfile_name)


