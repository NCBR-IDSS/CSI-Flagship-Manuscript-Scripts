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


# Python program to illustrate the intersection
# of two lists in most simple way
def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3


def query_bsi(ids, fields, query_field):
    # Set up the variables, bsi info
    cnf = os.path.expanduser('~/.my.cnf.bsi')
    url_session = 'https://rest.bsisystems.com/api/rest/EBMS/common/logon'
    url_reports = 'https://rest.bsisystems.com/api/rest/EBMS/reports/list'
    curl_get = "curl -s -X GET --header 'Accept: application/json' --header 'BSI-SESSION-ID: "

    user, pw = read_conf(cnf)
    session = get_bsi_session(url_session, user, pw)

    bsi = bsi_query(curl_get, url_reports, session, fields, ids, query_field)

    return bsi


# For Whole Genome Data:
# Reads in IDs from HGSC delivery files and gets family information
def get_all_family(received):

    print(str(received.shape[0]) + ' samples in raw directory.')

    # first get all family member IDs
    family_ids = query_bsi(received['DLM_LIS_Number'], ['Phenotips Family ID'], 'DLM LIS Number')
    family_ids = family_ids['Phenotips_Family_ID'].unique()
    # pd.DataFrame(family_ids).to_csv('fam_ids.csv')

    # Get IDs to merge from BSI
    bsi_fields = ['Phenotips ID', 'Batch Sent', 'Batch Received', 'Batch Sent Genome', 'Batch Received Genome',
                  'DLM LIS Number', 'Exome ID', 'Genome ID', 'CRIS Order Status', 'CRIS Order #', 'Archive',
                  'Active Status', 'Phenotips Family ID', 'Date of Enrollment']
    received_and_fam = query_bsi(family_ids, bsi_fields, 'Phenotips Family ID')

    # Clean up BSI data
    received_and_fam = received_and_fam[~received_and_fam['CRIS_Order_Status'].str.contains('cancel', case=False, na=False)]
    received_and_fam = received_and_fam[~received_and_fam['CRIS_Order_Status'].str.match('auto complete', case=False, na=False)]
    inactive = received_and_fam[~received_and_fam['Active status'].str.match('Active')]
    received_and_fam = received_and_fam[received_and_fam['Active status'].str.match('Active')]

    # For now, don't remove archived family and samples
    # received_and_fam = received_and_fam[~received_and_fam['Archive'].str.match('Archive')]
    # archived = received_and_fam[received_and_fam['Archive'].str.match('Archive')]
    received_and_fam.drop_duplicates(keep='first', inplace=True)
    received_and_fam.replace('', np.NaN, inplace=True)
    received_and_fam.rename(columns={'Exome ID': 'Exome_ID',
                                     'DLM LIS Number': 'DLM_LIS_Number',
                                     'Genome ID': 'Genome_ID',
                                     'Batch Sent Exome': 'Batch_Sent_Exome',
                                     'Batch Received Exome': 'Batch_Received_Exome',
                                     'Batch Sent Genome': 'Batch_Sent_Genome',
                                     'Batch Received Genome': 'Batch_Received_Genome'}, inplace=True)
    received_and_fam.to_csv('received_and_fam.csv', index=False)

    # if archived.shape[0] > 0:
    #     print('Warning! The following samples are Archived, but were delivered with the latest CIDR release: ')
    #     print(archived)
    #     print('')

    if inactive.shape[0] > 0:
        print('Warning! The following samples are Inactive, but were delivered with the latest CIDR release: ')
        print(inactive)
        print('')

    # Merge additional BSI fields and rearrange columns
    received = received.merge(received_and_fam, on='DLM_LIS_Number', how='left')
    received = received[['Exome_ID', 'Genome_ID_x', 'Genome_ID_y', 'Phenotips_ID_x', 'Phenotips_ID_y', 'DLM_LIS_Number',
                         'Batch_Sent_Exome', 'Batch_Received_Exome', 'Batch_Sent_Genome', 'Batch_Received_Genome',
                         'CRIS_Order#', 'Phenotips_Family_ID', 'Date of Enrollment', 'File_Prefix']]
    received['Notes'] = ''

    # Make sure HGSC linked correct Phenotips_ID to DLM_LIS_Number
    # Look for samples that were previously exome sequenced
    for i, row in received.iterrows():

        # flag mismatching Phenotips IDs
        if row['Phenotips_ID_x'] != row['Phenotips_ID_y']:
            print('Warning: ' + str(row['DLM_LIS_Number']) + ' associated with 2 Phenotips IDs! HGSC: ' +
                  str(row['Phenotips_ID_x']) + ' and BSI: ' + str(row['Phenotips_ID_y']) + '!')

        # check if sample was previously exome sequenced
        if not pd.isna(row['Batch_Sent_Exome']) and not pd.isna(row['Batch_Received_Exome']):
            received.at[i, 'Notes'] = received.at[i, 'Notes'] + ',WES sequenced'

        # check if sample was previously genome sequenced
        if not pd.isna(row['Genome_ID_y']) and str(row['Batch_Received_Genome']) != 'BATCH' + str(genome_batch):
            print('Warning: Sample ' + row['DLM_LIS_Number'] + ' (' + row['Phenotips_ID_x'] +
                  ') has already been assigned a Genome ID in BSI! HGSC: ' + row['Genome_ID_x'] +
                  ' and BSI: ' + str(row['Genome_ID_y']) + '!')
            received.at[i, 'Notes'] = received.at[i, 'Notes'] + ',WGS sequenced'

    received['Notes'] = received['Notes'].str.strip(' ,')
    received.drop(['Phenotips_ID_y', 'Genome_ID_y'], axis = 1, inplace = True)
    received.rename(columns = {'Phenotips_ID_x': 'Phenotips_ID',
                               'Genome_ID_x': 'Genome_ID'}, inplace = True)
    received.to_csv('received.csv')

    # Get sequenced family members
    all_family = received_and_fam[~received_and_fam['Phenotips_ID'].isin(received['Phenotips_ID'])]
    all_family.to_csv('all_family.csv')

    return received, all_family


# Makes masterkey for WES
def make_exome_masterkey(all_received, all_family):

    print('\nMaking masterkey...')

    received = all_received.copy()

    # Get sequenced and unsequenced family
    sequenced_fam = all_family[~all_family['Batch_Received_Exome'].isna()]
    sequenced_fam = sequenced_fam[~sequenced_fam['Exome_ID'].isna()]
    sequenced_fam = sequenced_fam[~sequenced_fam['DLM_LIS_Number'].isna()]
    sequenced_fam.to_csv('seq_fam_exome.csv')

    # Removes family members from future batches
    if not include_future and sequenced_fam.shape[0] > 0:
        sequenced_fam['Batch_Received_Exome'] = sequenced_fam['Batch_Received_Exome'].fillna('')
        sequenced_fam['Batch_Received_Exome'] = sequenced_fam['Batch_Received_Exome'].astype(str)
        sequenced_fam['Batch_Rec_Number'] = sequenced_fam['Batch_Received_Exome'].str.extract('(\d+)')
        sequenced_fam['Batch_Rec_Number'] = sequenced_fam['Batch_Rec_Number'].fillna('0')
        sequenced_fam['Batch_Rec_Number'] = sequenced_fam['Batch_Rec_Number'].astype(int)
        sequenced_fam = sequenced_fam[sequenced_fam['Batch_Rec_Number'] <= exome_batch]

    unsequenced_fam = all_family[~all_family['Phenotips_ID'].isin(sequenced_fam['Phenotips_ID'])]

    # if there's a previously WES version of a received sample, overwrite it
    # Remove that sample from family so it's synthetic exome reprocessed
    previous_wes = received[received['Notes'].str.contains('WES sequenced')]
    if not previous_wes.empty:
        print('Warning: ' + str(len(previous_wes)) + ' previously WES sequenced samples re-delivered')
        print('Will treat as new samples and redo synthetic exome:')
        print(previous_wes[['Phenotips_Family_ID', 'Phenotips_ID', 'Exome_ID', 'Batch_Received_Exome']])
        print('')
        sequenced_fam = sequenced_fam[~sequenced_fam['Phenotips_ID'].isin(previous_wes['Phenotips_ID'])]
        unsequenced_fam = unsequenced_fam[~unsequenced_fam['Phenotips_ID'].isin(previous_wes['Phenotips_ID'])]

    received['Exome_ID'] = received['Genome_ID']
    received['Batch_Sent_Exome'] = received['Batch_Sent_Genome']
    received['Batch_Received_Exome'] = 'BATCH' + str(exome_batch)
    sequenced_fam['File_Prefix'] = None

    masterkey = received.append(sequenced_fam, ignore_index = True, sort = False)
    print(str(received.shape[0]) + ' new samples')
    print(str(sequenced_fam.shape[0]) + ' sequenced family members')

    return masterkey, sequenced_fam, unsequenced_fam


def make_genome_masterkey(all_received, all_family):

    print('\nMaking masterkey...')
    received = all_received.copy()
    received['Batch_Received_Genome'] = 'BATCH' + str(genome_batch)

    # Get sequenced and unsequenced family
    sequenced_fam = all_family[~all_family['Batch_Received_Genome'].isna()]
    sequenced_fam = sequenced_fam[~sequenced_fam['Genome_ID'].isna()]
    sequenced_fam = sequenced_fam[~sequenced_fam['DLM_LIS_Number'].isna()]

    # Removes family members from future batches
    if not include_future and sequenced_fam.shape[0] > 0:
        sequenced_fam['Batch_Received_Genome'] = sequenced_fam['Batch_Received_Genome'].fillna('')
        sequenced_fam['Batch_Received_Genome'] = sequenced_fam['Batch_Received_Genome'].astype(str)
        sequenced_fam['Batch_Rec_Number'] = sequenced_fam['Batch_Received_Genome'].str.extract('(\d+)')
        sequenced_fam['Batch_Rec_Number'] = sequenced_fam['Batch_Rec_Number'].fillna('0')
        sequenced_fam['Batch_Rec_Number'] = sequenced_fam['Batch_Rec_Number'].astype(int)
        sequenced_fam = sequenced_fam[sequenced_fam['Batch_Rec_Number'] <= genome_batch]

    unsequenced_fam = all_family[~all_family['Phenotips_ID'].isin(sequenced_fam['Phenotips_ID'])]

    masterkey = received.append(sequenced_fam, ignore_index=True, sort=False)
    print(str(received.shape[0]) + ' new samples')
    print(str(sequenced_fam.shape[0]) + ' sequenced family members')

    return masterkey, sequenced_fam, unsequenced_fam


# writes table containing Phenotips_ID, Exome_ID, Batch_Received_Exome, CRIS_Order#
# of all samples received from HGSC for Xi to upload to BSI (one file for each split batch)
def write_bsi_info_wes(masterkey, filename):

    # Isolate only the new samples
    bsi = masterkey[masterkey['Batch_Received_Exome'].str.match('BATCH' + str(exome_batch), na = False)]
    bsi = bsi[['Phenotips_ID', 'Exome_ID', 'Batch_Received_Exome', 'CRIS_Order#']]
    bsi['Vendor'] = 'HGSC'
    bsi.rename(columns = {'CRIS_Order#': 'CRIS_Order_ID'}, inplace = True)
    bsi.to_csv(filename, index = False)


# writes table containing Phenotips_ID, Genome_ID, Batch_Received_Genome, CRIS_Order#
# of all samples received from HGSC for Xi to upload to BSI (one file for each split batch)
def write_bsi_info_wgs(masterkey, filename):

    # Isolate only the new samples
    bsi = masterkey[masterkey['Batch_Received_Genome'].str.match('BATCH' + str(genome_batch), na = False)]
    bsi = bsi[['Phenotips_ID', 'Genome_ID', 'Batch_Received_Genome', 'CRIS_Order#']]
    bsi['Vendor'] = 'HGSC'
    bsi.rename(columns = {'CRIS_Order#': 'CRIS_Order_ID'}, inplace = True)
    bsi.to_csv(filename, index = False)


# batch_num = Batch Number (i.e. 23)
# unsequenced_fam = family members of patients received in current batch (before splitting) that do not have data received/processed yet
# masterkey = masterkey file for one half of the current batch
# received = all samples pulled straight from the sample key file delivered by HGSC
# moved_samples = list of sample IDs added to current batch that were moved from previous batches into current
def make_sample_tracking_file(batch_num, unsequenced_fam, masterkey, received, moved_samples):

    masterkey = masterkey[['Phenotips_Family_ID', 'Phenotips_ID', 'Exome_ID', 'DLM_LIS_Number', 'CRIS_Order#', 'Batch_Sent', 'Batch_Received']]

    cnf, url_session, url_reports, curl_get = return_bsi_info()
    user, pw = read_conf(cnf)
    session = get_bsi_session(url_session, user, pw)

    batch_name = 'BATCH' + str(batch_num)
    batch_name_lower = 'Batch ' + str(batch_num)

    f = open('sample_tracking_summary_batch' + str(batch_num) + '.txt', 'w')
    f.write('******************************************\n*\n* Notes for samples released with Batch ' + str(batch_num) + '\n' +
            '*\n******************************************\n\n')

    fields = ['Phenotips ID', 'CRIS Order #', 'Batch Sent', 'Batch Received', 'CRIS Order Status', 'Active Status']
    sent = bsi_query(curl_get, url_reports, session, fields, [batch_name], 'Batch Sent')
    sent = sent[~sent['CRIS_Order_Status'].str.contains('cancel', case = False, na = False)]
    sent = sent[~sent['CRIS_Order_Status'].str.match('auto complete', case = False, na = False)]
    sent = sent[sent['Active status'].str.match('Active')]
    sent.drop_duplicates(keep = 'first', inplace = True)
    num_sent = len(sent['CRIS_Order#'].unique())

    f.write('Total number of samples sent to HGSC with Batch ' + str(batch_num) + ': ' + str(num_sent) + '\n')
    # f.write('** Samples sent together in ' + batch_name_lower + ' were split into two separate batches when received **\n\n')

    num_released_all = received.shape[0]
    f.write('Total number of new samples in latest HGSC release (before split): ' + str(num_released_all) + '\n\n')

    num_released = masterkey[masterkey['Batch_Received'].str.match(batch_name, na = False)].shape[0]
    f.write('Number of new samples released from HGSC in ' + batch_name_lower + ': ' + str(num_released) + '\n\n')

    num_in_masterkey = masterkey.shape[0]
    f.write('Total number of samples in the masterkey file: ' + str(num_in_masterkey) + "\n")
    f.write("Here's the breakdown: \n")
    f.write('---------------------------------------------------------------------------------------------\n\n')

    num_released_from_sent = len(intersection(sent['CRIS_Order#'].values.tolist(), masterkey['CRIS_Order#'].values.tolist()))
    f.write(str(num_released_from_sent) + ' sample(s) sent to HGSC in ' + batch_name_lower + ' were released with ' + batch_name_lower + '.\n\n')

    ##### Special Case: remove samples that are being added/moved from other batches into current batch from previous batches.
    # We'll make a separate section to outline those samples

    sent_in_other_and_released = masterkey[~masterkey['Batch_Sent'].str.match(batch_name, na = False)]
    sent_in_other_and_released = sent_in_other_and_released[sent_in_other_and_released['Batch_Received'].str.match(batch_name, na = False)]
    sent_in_other_and_released = sent_in_other_and_released[~sent_in_other_and_released['Phenotips_ID'].isin(moved_samples)]
    f.write(str(sent_in_other_and_released.shape[0]) + ' sample(s) sent to HGSC in other batches were released with ' + batch_name_lower + ':\n')
    f.write(sent_in_other_and_released.to_string(index = False) + '\n\n')

    # This section outlines samples moved to current batch:
    num_moved = len(moved_samples)
    if num_moved > 0:
        f.write(str(num_moved) + ' sample(s) released in previous batches were moved to ' + batch_name_lower + ':\n')
        f.write(masterkey[masterkey['Phenotips_ID'].isin(moved_samples)].to_string(index = False) + '\n\n')

    # all family members in the masterkey file (they're sequenced and have data on file)
    fam_members_released = masterkey[~masterkey['Batch_Received'].str.match(batch_name)]
    fam_members_released = fam_members_released[~fam_members_released['Phenotips_ID'].isin(moved_samples)]
    # fam_members_released = fam_members_released[['']]
    if fam_members_released.empty:
        f.write('0 sample(s) released in previous batches are family members of sample(s) released with ' + batch_name_lower + '.\n\n')
    else:
        f.write(str(fam_members_released.shape[0]) + ' sample(s) released in previous batches are family members of sample(s) released in ' + batch_name_lower + ':\n')
        f.write(fam_members_released.to_string(index = False) + '\n\n')

    f.write('---------------------------------------------------------------------------------------------\n\n')

    # get any unsequenced family members that are related to people in this half of the batch
    fam_members_unreleased = unsequenced_fam[unsequenced_fam['Phenotips_Family_ID'].isin(masterkey['Phenotips_Family_ID'])]
    if fam_members_unreleased.empty:
        f.write('There are no unreleased family members of sample(s) released in ' + batch_name_lower + ' on file.\n\n')
    else:
        f.write(str(fam_members_unreleased.shape[0]) + ' family members have been consented but not yet released.\n')
        fam_members_unreleased = fam_members_unreleased[['Phenotips_Family_ID', 'Phenotips_ID']]
        f.write(fam_members_unreleased.to_string(index = False))
        f.write('\n\n')

    sent_not_released = sent[~sent['CRIS_Order#'].isin(masterkey['CRIS_Order#'])]
    if sent_not_released.empty:
        f.write('All non-canceled samples sent to HGSC in ' + batch_name_lower + ' have been released in ' + batch_name_lower + '\n\n')
    else:
        f.write(str(sent_not_released.shape[0]) + ' sample(s) sent to HGSC in ' + batch_name_lower + ' were not released in ' + batch_name_lower + ':\n')
        f.write(sent_not_released.to_string(index=False) + '\n\n')

    f.close()


# Write shell script to symbolically link old bam files to new directory
# masterkey: one of the split masterkeys
def write_bam_link_script(masterkey, batch_received, field):
    batch_name = 'BATCH' + str(batch_received)
    previous_batches = masterkey[~masterkey[field].str.match(batch_name)]

    # Write script to link bams from previous batches (family/added samples) on LOCUS
    locus_dir = '/hpcdata/dir/CSI_DATA_PROCESSED'
    fname = 'link_previous_bams_' + 'batch' + str(batch_received) + '.sh'

    script = open(fname, "w")
    script.write("#!/bin/sh\nset -e\n\n")

    for i in previous_batches.index.values.tolist():
        pbatchdir = previous_batches.loc[i][field]
        pbatchdir = pbatchdir.replace('BATCH0', 'BATCH')
        phen_id = previous_batches.loc[i]['Phenotips_ID']

        bamlink = "ln -s {}/{}/BAM/{}.bam {}/{}/BAM/\n".format(locus_dir, pbatchdir, phen_id, locus_dir,
                                                               'BATCH' + str(batch_received))
        bailink = "ln -s {}/{}/BAM/{}.bam.bai {}/{}/BAM/\n".format(locus_dir, pbatchdir, phen_id, locus_dir,
                                                                   'BATCH' + str(batch_received))
        script.write(bamlink)
        script.write(bailink)
        script.write("\n")

    script.close()


def main():

    global exome_batch
    global genome_batch
    global include_future

    #
    # Usage statement
    #
    parseStr = 'Reads a list of available files, performs quality control on the data,\n\
    and outputs the files needed for GRIS.\n'

    parser = argparse.ArgumentParser(description=parseStr, formatter_class=RawTextHelpFormatter)
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

    if len(crams) == 0:
        print('Cannot locate CRAM files!')
        return

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

    # # Manual for B101
    # received = pd.read_csv('/Users/kuramvs/Documents/csi_to_gris/101/sample_mapping.csv', dtype = str)
    # received['File_Prefix'] = None

    received, all_family = get_all_family(received)

    if exome_batch:
        print('--------------------- Generating WES files ---------------------')
        wes_masterkey, wes_sequenced_fam, wes_unsequenced_fam = make_exome_masterkey(received, all_family)
        write_bsi_info_wes(wes_masterkey, 'Batch' + str(exome_batch) + '_BaylorExoIDs.csv')
        wes_masterkey = wes_masterkey[['File_Prefix', 'DLM_LIS_Number', 'Phenotips_ID', 'Exome_ID', 'Batch_Sent_Exome', 'Batch_Received_Exome']]
        wes_masterkey.to_csv('masterkey_wes_batch' + str(exome_batch) + '.txt', index = False, sep = '\t')
        write_bam_link_script(wes_masterkey, exome_batch, 'Batch_Received_Exome')

    if genome_batch:
        print('\n--------------------- Generating WGS files ---------------------')
        wgs_masterkey, wgs_sequenced_fam, wgs_unsequenced_fam = make_genome_masterkey(received, all_family)

        write_bsi_info_wgs(wgs_masterkey, 'Batch' + str(genome_batch) + '_BaylorGenoIDs.csv')
        wgs_masterkey = wgs_masterkey[['File_Prefix', 'DLM_LIS_Number', 'Phenotips_ID', 'Genome_ID', 'Batch_Sent_Genome', 'Batch_Received_Genome']]
        wgs_masterkey.to_csv('masterkey_wgs_batch' + str(genome_batch) + '.txt', index = False, sep = '\t')
        write_bam_link_script(wgs_masterkey, genome_batch, 'Batch_Received_Genome')


if __name__ == '__main__':
    main()
