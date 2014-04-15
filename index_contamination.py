#!/usr/bin/env python
__author__ = 'aberlin'

import argparse
import re
import sys
import pandas as pd


def parse_cmdline_params(arg_list=None):
    """Parses commandline arguments.
    """

    description = "Look for Contamination in index panels"

    #Create instance of ArgumentParser
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument("-f",
                        "--file",
                        help="input read mapping file",
                        required=True)

    parser.add_argument("-t",
                        "--target_min",
                        help="Minimum fraction on target to look for contamination",
                        default=0.999,
                        required=False)

    parser.add_argument("-m",
                        "--syn_min",
                        help="Minimum number of reads needed to count a synthetic",
                        default=10,
                        required=False)

    # Parse options
    opts = parser.parse_args(args=arg_list)

    return opts


def load_read_data(in_file):

    # Load data from tab-delmin file
    df = pd.DataFrame.from_csv(in_file, sep='\t')

    # Strip whitespace from Synthetics
    df["Synthetics"] = df["Synthetics"].map(str.strip)

    # Calculate Total Reads for each index
    df["Total_Reads"] = df.sum(axis=1)

    # Create empty columns for future calculations - controls order
    df["Fraction_On_Target"] = None
    df["Contamination"] = None
    df["Obs_Synthetics"] = None
    df["Obs_Primers"] = None
    df["Expected_Synthetics_Present"] = None
    df["Expected_Primers_Present"] = None
    df["Contaminating_Index"] = None

    return df


def get_fraction_on_target(df_index):
    on_reads = 0
    for target in df_index["Synthetics"]:
        if target in df_index:
            on_reads += df_index[target]

    on_reads += df_index[df_index["Primers"]]

    return on_reads/float(df_index["Total_Reads"])


def expected_present(df_index, type):

    for expected in df_index[type]:
        if not re.search(expected, df_index["Obs_" + type]):
            return False
    else:
        return True


def find_observed(df_index, min_to_count):
    obs_syn = []
    obs_primers = []
    for col in df_index[2:8].index.values:
        if df_index[col] >= min_to_count:
            obs_syn.append(col)

    for col in df_index[8:11].index.values:
        if df_index[col] >= min_to_count:
            obs_primers.append(col)

    return ','.join(obs_syn), ','.join(obs_primers)


def find_contaminating_index(df_index, index_def):
    index_list = []
    syn_list = []
    for index in index_def.index:
        all_match = True

        if index == df_index.name:
            continue

        if re.search(index_def.ix[index, "Primers"], df_index["Obs_Primers"]):
           # print "Primer Match -  " + str(index) + " : " + index_def.ix[index, "Primers"] + " and " + df_index["Obs_Primers"]

            for syn in index_def.ix[index, "Synthetics"]:
                if not re.search(syn, df_index["Obs_Synthetics"]):
                    all_match = False
                    continue

            if all_match:
               # print "Syn Match -  " + str(index) + " : " + index_def.ix[index, "Synthetics"] + " and " + df_index["Obs_Synthetics"]
                index_list.append(str(index))
                syn_list.append(str(index_def.ix[index, "Synthetics"]))

        else:
            continue

    return index_list, syn_list


def main(args):
    opts = parse_cmdline_params(args[1:])

    data = load_read_data(opts.file)

    for index in data.index:

        obs_syn, obs_primers = find_observed(data.ix[index], opts.syn_min)
        data.ix[index, "Obs_Synthetics"] = obs_syn
        data.ix[index, "Obs_Primers"] = obs_primers

        data.ix[index, "Expected_Synthetics_Present"] = expected_present(data.ix[index], "Synthetics")
        data.ix[index, "Expected_Primers_Present"] = expected_present(data.ix[index], "Primers")

        fract_on_target = get_fraction_on_target(data.ix[index])
        data.ix[index, "Fraction_On_Target"] = fract_on_target

        if fract_on_target < opts.target_min:
            data.ix[index, "Contamination"] = "Contaminated"
            contam_index, contam_syn = find_contaminating_index(data.ix[index], data.iloc[:, :2])
            data.ix[index, "Contaminating_Index"] = ','.join(contam_index)
            data.ix[index, "Contaminating_Syn"] = ','.join(contam_syn)

        else:
            data.ix[index, "Contamination"] = "PASS"
        #print data.ix[index]

    data.to_csv("output.tsv", sep='\t')
    print data

if __name__ == "__main__":
    main(sys.argv)