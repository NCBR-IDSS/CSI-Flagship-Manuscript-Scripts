import sys
import os
import re
import datetime
import pandas as pd
import numpy as np
import argparse
from argparse import RawTextHelpFormatter
from ncbr_huse import send_update, err_out, pause_for_input
import subprocess
import json
import fnmatch
import glob
from ncbr_bsi import read_conf, send_curl, get_bsi_session, get_bsi_name, bsi_query, return_bsi_info
from generate_seqr_ped import generate_genrptlinks_ped
 
global exome_batch
global genome_batch
global include_future

    #
    # Usage statement
parseStr = 'Reads a list of available files, performs quality control on the data,\n\and outputs the files needed for GRIS.\n'

parser = argparse.ArgumentParser(description=parseStr, formatter_class=RawTextHelpFormatter)
print('here\n')
parser.add_argument('-d', '--dir', required=False, type=str, default="", help='Directory containing HGSC .cram files')
parser.add_argument('-e', '--exome', required=False, type = int, default=False, help='Specify exome batch number. Will generate processing files')
parser.add_argument('-g', '--genome', required=False, type = int, default=False, help='Specify genome batch number. Will generate processing files')
parser.add_argument('-f', '--future', required=False, action='store_true', default=False, help='Includes family members from future batches in masterkey')

args = parser.parse_args()
directory = args.dir
exome_batch = args.exome
genome_batch = args.genome
include_future = args.future

    # Mac/Linux only
if directory != "" and not directory.endswith('/'):
    directory = directory + '/'

crams = glob.glob(directory + "BCM_P*_*_*_*.cram")
    #example cram file name: BCM_P0006325_O4081393_1000086314_1.cram
    
if len(crams) == 0:
    print('Cannot locate CRAM files!')

received = pd.DataFrame({'File_Prefix': crams,
                             'Phenotips_ID': None,
                             'DLM_LIS_Number': None,
                             'Genome_ID': None})

received['File_Prefix'] = received['File_Prefix'].str.replace(directory, '', regex = False)
received['File_Prefix'] = received['File_Prefix'].str.replace('.cram', '', regex = False)
split = received['File_Prefix'].str.split("_", expand = True)

received['Phenotips_ID'] = split[1]
received['DLM_LIS_Number'] = split[2]
received['Genome_ID'] = split[3]

received.drop_duplicates(inplace = True)
received.to_csv('raw_received.csv')

    
