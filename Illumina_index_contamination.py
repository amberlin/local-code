#!/usr/bin/env python
"""Evaluate Illumina barcodes for contamination

  NOTE:  Currently this is hardcoded in 2 places for synthetics A-F
  Adjust as necessary for future run until I figure out a dynamic way
  to know which columns are the synthetics.

"""

__author__ = "Aaron Berlin"
__copyright__ = "Copyright 2014, Enzymatics"
__credits__ = ["Aaron Berlin"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Aaron Berlin"
__email__ = "aberlin@enzymatics.com"
__status__ = "Development"

import argparse
import re
import sys
import pandas as pd
import numpy as np


def parse_cmdline_params(arg_list=None):
    """Parses commandline arguments.
    """

    description = "Look for Contamination in Illumina index panels"

    #Create instance of ArgumentParser
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument("-f",
                        "--file",
                        help="input read mapping file",
                        required=True)

    parser.add_argument("-t",
                        "--target_min",
                        help="Minimum fraction on target to look for contamination",
                        default=0.1,
                        type=float,
                        required=False)

    parser.add_argument("-m",
                        "--syn_min",
                        help="Minimum percent of total reads needed to count a synthetic",
                        type=float,
                        default=0.10,
                        required=False)

    # Parse options
    opts = parser.parse_args(args=arg_list)

    return opts


def load_read_data(in_file):

    # Load data from tab-delmin file
    df = pd.read_csv(in_file, sep='\t', dtype={'designed_target': np.str})

    # Strip whitespace from Synthetics
    df["Synthetics"] = df["designed_target"]  # .map(str.strip)

    # Calculate Total Reads for each index
    df["Total_Reads"] = df.loc[:, 'A':'F'].sum(axis=1, numeric_only=True)

    # Create empty columns for future calculations - controls order
    df["Fraction_On_Target"] = None
    df["Pct_of_Real_Mean"] = None
    df["Contamination"] = None
    df["Obs_Synthetics"] = None
    df["Expected_Synthetics_Present"] = None
    df["Contaminating_Index"] = None

    return df


### Return the fraction of reads that align to the expected synthetics
def get_fraction_on_target(df_index):
    on_reads = 0
    for target in df_index["Synthetics"]:
        if target in df_index:
            on_reads += df_index[target]

    return on_reads/float(df_index["Total_Reads"])


### Check that the expected synthetics are present
def expected_present(df_index, syn_type):

    if str(df_index[syn_type]) == 'nan':
        return None
    for expected in df_index[syn_type]:
        if not re.search(expected, df_index["Obs_" + syn_type]):
            return False
    else:
        return True


### Find all observed synthetics
def find_observed(df_index, min_to_count):
    obs_syn = []

    # index is hardcoded for current experiment - could be smarter about this
    for col in df_index.loc['A':'F'].index.values:
        if df_index[col] > df_index["Total_Reads"] * min_to_count:
            obs_syn.append(col)

    return ','.join(obs_syn)


### Compare observed synthetics to all expected synthetics to find all contaminating libraries
def find_contaminating_index(df_index, index_def):
    index_list = []
    syn_list = []
    for index in index_def.index:

        # Skip itself
        if index == df_index.name:
            continue

        # If all of the index synthetics are observed then if could be a contaminate
        if set(index_def.ix[index, "Synthetics"]).issubset(set(df_index["Obs_Synthetics"].split(','))):

            index_list.append(str(index_def.ix[index, "Library_Name"].split("_")[0]))
            syn_list.append(str(index_def.ix[index, "Synthetics"]))

    return index_list, syn_list


def main(args):
    opts = parse_cmdline_params(args[1:])

    data = load_read_data(opts.file)

    # Create a new dataframe with just the real MBC libraries
    real_lib_data = data[data['Library_Name'].str.contains("MBC")]

    # Calc the percent of the mean reads in the real MBC libraries
    data["Pct_of_Real_Mean"] = (data["Total_Reads"] / real_lib_data['Total_Reads'].mean()) * 100

    for index in data.index:
        fract_on_target = None
        syn_expected = set()

        # Lookup all the synthetics with read depth above the provided threshold
        obs_syn = find_observed(data.ix[index], opts.syn_min)
        data.ix[index, "Obs_Synthetics"] = obs_syn

        # If synthetics are present - check they match what is expected and calc the on target reads
        if not str(data.ix[index, "Synthetics"]) == 'nan':
            data.ix[index, "Expected_Synthetics_Present"] = expected_present(data.ix[index], "Synthetics")
            fract_on_target = get_fraction_on_target(data.ix[index])
            syn_expected = set(str(data.ix[index, "Synthetics"]))

        data.ix[index, "Fraction_On_Target"] = fract_on_target

        # If we have enough reads and there are observed synthetics that are not expected - contamination
        if data.ix[index, "Pct_of_Real_Mean"] > opts.target_min and not set(obs_syn.split(',')) == syn_expected:

            # Look up all possible contaminating indexes
            contam_index, contam_syn = find_contaminating_index(data.ix[index], real_lib_data)

            # Not contamination if there are no valid indexes with the observed synthetics
            if not contam_index:
                data.ix[index, "Contamination"] = "PASS"
            else:
                data.ix[index, "Contaminating_Index"] = ','.join(contam_index)
                data.ix[index, "Contaminating_Syn"] = ','.join(contam_syn)
                data.ix[index, "Contamination"] = "Contaminated"

            # Compare observed synthetics to the expected synthetics for the real library with the same P5
            real_expected_syn = real_lib_data.query('P5 == ' + str(data.ix[index, "P5"]))["Synthetics"]

            # If all the observed synthetics match expected synthetics from the real library mark them
            if not real_expected_syn.empty and set(obs_syn.split(',')) == set(real_expected_syn.iloc[0]):
                data.ix[index, "Expected_Synthetics_Present"] = "MATCH"

        else:
            data.ix[index, "Contamination"] = "PASS"

    data.to_csv("output.tab", sep='\t')
    print data

if __name__ == "__main__":
    main(sys.argv)