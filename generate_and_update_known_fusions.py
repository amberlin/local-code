#!/usr/bin/env python
"""Compare different fusion DB and lists"""

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
import os
from collections import defaultdict
from Bio import Entrez
import urllib2
import pickle
import sqlite3 as lite

Entrez.email = "aberlin@enzymatics.com"


# Parse command line arguments
def parse_cmdline_params(arg_list=None):
    """Parses commandline arguments.
    """

    description = "Compare different fusion DBs."

    #Create instance of ArgumentParser
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument("-c",
                        "--cosmic_file",
                        type=argparse.FileType('r'),
                        help="COSMIC fusion file",
                        required=False)

    parser.add_argument("-d",
                        "--database_name",
                        help="Database Name",
                        required=False)

    parser.add_argument("-i",
                        "--chitars_file",
                        type=argparse.FileType('r'),
                        help="ChiTars breakpoint file",
                        required=False)

    parser.add_argument("-e",
                        "--ensembl_file",
                        type=argparse.FileType('r'),
                        help="Ensembl transcripts gtf file",
                        required=False)

    parser.add_argument("-p",
                        "--pickled_data",
                        help="Data from a previous run",
                        required=False)

    parser.add_argument("-t",
                        "--test_data",
                        type=argparse.FileType('r'),
                        help="Test data",
                        required=False)

    # Parse options
    opts = parser.parse_args(args=arg_list)

    return opts


# Allows for automatically nested defaultdicts
def my_dd():
    return defaultdict(my_dd)


# Parse the COSMIC fusion from hgvs format
def __parse_fusion(fusion, return_type="All"):

    # Super complicated regex to deal with the hgvs format
    regex = re.compile("(?P<gene_names>[a-zA-Z0-9]+)\{(?P<transcripts>\w+\.?\d?)\}:r.(?P<coords>.*?)(?=_[A-Z]|_o|$)")

    parsed_results = defaultdict(list)
    for result in regex.finditer(fusion):
        parsed_results["pairs"].append([result.group("gene_names"),
                                        result.group("transcripts")])

        parsed_results["full_data"].append([result.group("gene_names"),
                                            result.group("transcripts"),
                                            result.group("coords")])
        for group in result.groupdict():
            parsed_results[group].append(result.group(group))

    if parsed_results:

        if return_type == "genes":
            return parsed_results['gene_names']
        elif return_type == "transcripts":
            return parsed_results['transcripts']
        elif return_type == "pairs":
            return parsed_results['pairs']
        else:
            return parsed_results
    else:
        print "WARN: Fusion not able to be parsed - " + fusion
        return None


# Return the all information available from the COSMIC fusion id
def get_all_from_fusion(fusion):
    # Returns a list of [gene_name, transcript, coordinates in gene] for each transcript in the fusion
    return __parse_fusion(fusion, "All")


# Return the transcripts from the COSMIC fusion id
def get_transcripts(fusion):
    return __parse_fusion(fusion, "transcripts")


# Return the gene name and transcript pairs from the COSMIC fusion id
def get_gene_transcript_pairs(fusion):
    return __parse_fusion(fusion, "pairs")


# Return the gene names from the COSMIC fusion id
def get_genes(fusion):
    return __parse_fusion(fusion, "genes")


# Parse genes from FusionExport File
def parse_cosmic(cosmic_lines):
    gene_list = []
    seen_list = {}
    for line in cosmic_lines:
        fields = re.split("\t", line)
        if len(fields) >= 8:
            fusion = fields[7]
            if not fusion in seen_list:
                genes = get_genes(fusion)
               # print genes
                gene_list.append(genes)
                seen_list[fusion] = 1

    return gene_list


def get_exons_involved(fusion_coordinates, exon_coordinates):

    if not len(fusion_coordinates) == len(exon_coordinates):
        print "WARN: Fusion and Exons coordinates are different lengths ... something is crossed up"

    exon_list = []

    #print fusion_coordinates
    #print exon_coordinates

    for i, transcript_coords in enumerate(fusion_coordinates):
        breakpoint_edge_coords = []
        comments = ""

        if transcript_coords == "?":
            exon_list.append("Unknown")
            continue

        if re.search("\(", transcript_coords):
            exon_list.append("Unknown")
            print "Weird coord w/ (): " + transcript_coords
            print fusion_coordinates
            continue

        # Parse the fusion coordinates into breakpoints that need to be looked up
        coords = transcript_coords.split("_")

        if len(coords) == 3:
            comments = coords[2]
            if not re.search("ins|del", coords[2]):
                print "Weird coord: " + transcript_coords
                print fusion_coordinates

        if i == 0:
            breakpoint_edge_coords.append(coords[1])
        elif i == len(fusion_coordinates) - 1:
            breakpoint_edge_coords.append(coords[0])
        else:
            breakpoint_edge_coords.append(coords[0])
            breakpoint_edge_coords.append(coords[1])

        # Lookup each edge from the corresponding list of exon coordinates
        for breakpoint_edge in breakpoint_edge_coords:
            #print "Coord: " + str(breakpoint_edge)

            # If the edge coordinate has a + or - then the break point is actually X bases in to the intron
            if re.search("\+", breakpoint_edge):
                breakpoint_edge, bp_into_intron = breakpoint_edge.split("+")
                comments += "(intron +%sbp)" % bp_into_intron

            if re.search("\-", breakpoint_edge):
                breakpoint_edge, bp_into_intron = breakpoint_edge.split("-")
                comments += "(intron -%sbp)" % bp_into_intron

            # If breakpoint edge coordinate has some weird text lets grab it and skip it for now
            try:
                int(breakpoint_edge)
            except ValueError:
                print "Weird breakpoint: " + breakpoint_edge
                exon_list.append("Unknown")
                continue

            exon_coords = exon_coordinates[i]
            # If we don't have an exon coordinate list for the transcript
            if not exon_coords:
                exon_list.append("Unknown")
                continue

            for j in range(len(exon_coords) - 1):
                if int(breakpoint_edge) > int(exon_coords[len(exon_coords) - 1]) and re.search("ins", comments):
                    exon_list.append(str(len(exon_coords)-1) + comments)
                    break

                if int(breakpoint_edge) > int(exon_coords[len(exon_coords) - 1]):
                    exon_list.append(">"+str(len(exon_coords)-1))
                    break

                if int(exon_coords[j]) < int(breakpoint_edge) <= int(exon_coords[j+1]):
                    exon_list.append(str(j+1) + comments)
                    #print "Exon : " + str(j+1)
                    break

    return exon_list


# Parse the gene record returned by Entrez
def parse_gene_record(record):
    omim_id = None
    locus_id = None
    ensembl_id = None

    try:
        gene_handle = Entrez.efetch(db="gene", rettype="xml", id=record["IdList"][0])
        gene_record = Entrez.read(gene_handle)[0]
    except urllib2.HTTPError:
        print "WARN: " + str(record["IdList"][0]) + " is not valid"
        return None
    except IOError:
        print "Problem connecting to NCBI"
        return None

    gene_record_locus = gene_record["Entrezgene_locus"][0]["Gene-commentary_seqs"][0]["Seq-loc_int"]["Seq-interval"]

    # Pull start, stop and strand
    gene_start_coord = gene_record_locus["Seq-interval_from"]
    gene_end_coord = gene_record_locus["Seq-interval_to"]
    gene_strand = gene_record_locus["Seq-interval_strand"]["Na-strand"].attributes["value"]

    # Pull reference information
    gene_ref_accession = gene_record["Entrezgene_locus"][0]["Gene-commentary_accession"]
    gene_chromosome = gene_record["Entrezgene_source"]["BioSource"]["BioSource_subtype"][0]["SubSource_name"]
    cytoband_location = gene_record["Entrezgene_gene"]["Gene-ref"].get("Gene-ref_maploc", None)

    # Pull MIM and Ensembl link out ids
    for db in gene_record["Entrezgene_gene"]["Gene-ref"]["Gene-ref_db"]:
        if db['Dbtag_db'] == "MIM":
            omim_id = db['Dbtag_tag']['Object-id'].get('Object-id_id', None)
        if db['Dbtag_db'] == "Ensembl":
            ensembl_id = db['Dbtag_tag']['Object-id'].get('Object-id_str', None)

    # Pull NCBI LocusID
    for db in gene_record['Entrezgene_unique-keys']:
        if db['Dbtag_db'] == "LocusID":
            locus_id = db['Dbtag_tag']['Object-id'].get('Object-id_id', None)

    #print omim_id, locus_id, str(ensembl_id)
    #print cytoband_location, gene_ref_accession, gene_chromosome, gene_start_coord, gene_end_coord, gene_strand
    return [omim_id, locus_id, cytoband_location, gene_ref_accession, str(ensembl_id),
            gene_chromosome, gene_start_coord, gene_end_coord, str(gene_strand)]


# Run the search on Entrez - NCBI
def run_search(query):
    search_handle = Entrez.esearch(db="gene", term=query)
    record = Entrez.read(search_handle)

    if int(record["Count"]) == 1:
        return parse_gene_record(record)

    # Warn if no records are found
    elif int(record["Count"]) == 0:
        print "WARN: No genes found: " + str(record["WarningList"].get("OutputMessage")) + \
              ": Search Attempted: " + str(record["QueryTranslation"])
        #print record
        return None

    # Warn if more than one record is found
    else:
        #print record
        print "WARN: More than one gene ID found: " + record["Count"] + " records found" + \
              ": Search Attempted: " + str(record["QueryTranslation"])
        return None


# Lookup information by gene name
def lookup_by_gene(gene_name):

    return run_search(gene_name + "[gene name] AND Homo[organism] AND alive[prop]")


# Lookup information by transcript id from Ensembl gtf and NCBI queries
def lookup_transcript(transcript, gene_name, ensembl_transcript_details):

    # remove the version from accession numbers for better search results
    if re.search('\.\d+$', transcript):
        transcript = re.sub('\.\d+$', "", transcript)

    # Search for transcript and ignore discontinued records
    ncbi_info = run_search(transcript + " AND alive[prop]")

    # If the transcript is not found look the gene name
    if not ncbi_info:
        print "Trying to lookup by gene name"
        ncbi_info = lookup_by_gene(gene_name)

    ensembl_info = None
    exon_info = []

    if transcript in ensembl_transcript_details:
        ensembl_info = ensembl_transcript_details[transcript]["details"]
        exon_info = ensembl_transcript_details[transcript]["exons"]

    if not ncbi_info and not ensembl_info:
        return None, None

    if not ensembl_info:
        return ncbi_info, None

    if not ncbi_info:
        ncbi_info = ["", "", "", ""]

    return ncbi_info[:4] + ensembl_info, exon_info


# Parse the Ensembl gtf annotations for transcript coordinates
def parse_transcript_gtf(transcript_file):
    gene_id = None
    trans_id = None
    transcript_details = {}
    running_total = 0

    for line in transcript_file:
        if re.match("#", line):
            continue

        fields = line.split('\t')
        chromosome = fields[0]
        start_coord = fields[3]
        end_coord = fields[4]
        strand = fields[6]
        comments = fields[8].split(';')

        if fields[2] == "transcript":

            for entry in comments:
                if re.search("gene_id", entry):
                    name, gene_id = re.sub("\"", "", entry).split()
                if re.search("transcript_id", entry):
                    name, trans_id = re.sub("\"", "", entry).split()
            #transcript_details[trans_id] = [gene_id, chromosome, start_coord, end_coord, strand]
            transcript_details[trans_id] = {}
            transcript_details[trans_id]["exons"] = [0]
            transcript_details[trans_id]["details"] = [gene_id, chromosome, start_coord, end_coord, strand]
            running_total = 0

        if fields[2] == "exon":

            for entry in comments:
                if re.search("transcript_id", entry):
                    name, trans_id = re.sub("\"", "", entry).split()

            running_total += int(end_coord) - int(start_coord) + 1

            # Store the coordinate of the the exon end on the transcript
            transcript_details[trans_id]["exons"].append(running_total)

    return transcript_details


#Parse genes and PubMed reference from CosmicFusionExport_vXX.tsv file
def parse_cosmic_basic(cosmic_lines, ensembl_transcript_details):
    gene_list = my_dd()  # nested defaultdict
    transcript_details = {}
    seen_fusion_list = {}
    seen_cosmic_ids = {}

    for line in cosmic_lines:
        fields = re.split("\t", line)

        if len(fields) >= 8:
            cosmic_id = fields[7]
            fusion_id = fields[8]

            # Store the Pubmed reference id if present
            if len(fields) >= 12:
                reference = fields[12]
            else:
                reference = "NA"

            # Skip identical fusions with the same pubmed reference
            if not fusion_id+reference in seen_fusion_list:

                seen_fusion_list[fusion_id+reference] = 1
                fusion_details = get_all_from_fusion(fusion_id)

                if fusion_details:
                    #print "\n" + fusion_id + " : " + reference
                    gene_pair = ":".join(fusion_details['gene_names'])

                    #Just pubmed id if we have already seen an identical fusion
                    if gene_pair in gene_list and fusion_id in gene_list[gene_pair]:
                        gene_list[gene_pair][fusion_id]["references"]["pubmed"].append(reference)
                        continue

                    # Store fusion information
                    gene_list[gene_pair][fusion_id]["references"]["pubmed"] = [reference]
                    gene_list[gene_pair][fusion_id]["references"]["cosmic_id"] = [cosmic_id]
                    gene_list[gene_pair][fusion_id]["transcript_details"] = []
                    gene_list[gene_pair][fusion_id]["exon_details"] = []
                    gene_list[gene_pair][fusion_id]["exons-involved"] = []
                    gene_list[gene_pair][fusion_id]["sample_count"] = 1

                    seen_cosmic_ids[cosmic_id] = 1

                    # Store information for each transcript in the fusion
                    for gene_name, transcript in fusion_details['pairs']:
                        if not transcript in transcript_details:

                            # Query NCBI for transcripts details
                            details, exon_details = lookup_transcript(transcript, gene_name, ensembl_transcript_details)

                            # Save the transcript details so we don't have to query NCBI again
                            transcript_details[transcript] = {}
                            transcript_details[transcript]["details"] = details
                            transcript_details[transcript]["exons"] = exon_details

                        #print "\t" + transcript + ": " + str(transcript_details[transcript])
                        gene_list[gene_pair][fusion_id]["transcript_details"].append(transcript_details[transcript]["details"])
                        gene_list[gene_pair][fusion_id]["exon_details"].append(transcript_details[transcript]["exons"])

                    gene_list[gene_pair][fusion_id]["exons-involved"] = \
                        get_exons_involved(fusion_details['coords'], gene_list[gene_pair][fusion_id]["exon_details"])

                    #print gene_list[gene_pair][fusion_id]
            else:

                gene_list[":".join(get_genes(fusion_id))][fusion_id]["sample_count"] += 1

            # Store new COSMIC id for a previously seen fusion
            if not cosmic_id in seen_cosmic_ids:
                gene_pair = ":".join(get_genes(fusion_id))
                gene_list[gene_pair][fusion_id]["references"]["cosmic_id"].append(cosmic_id)
                seen_cosmic_ids[cosmic_id] = 1

    return gene_list


# Merge COSMIC and ChiTars data structures
def merge_data(cosmic_data, chitars_data):
    chitars_only = 0
    shared = 0
    total_cosmic = len(cosmic_data)

    print "Total gene pairs in COSMIC: %d" % total_cosmic
    print "Total gene pairs in ChiTars: %d" % len(chitars_data)

    for gene_pair in chitars_data:
        if not gene_pair in cosmic_data:
            cosmic_data[gene_pair] = chitars_data[gene_pair]
            chitars_only += 1
        else:
            for breakpoint in chitars_data[gene_pair]:
                cosmic_data[gene_pair][breakpoint] = chitars_data[gene_pair][breakpoint]
            shared += 1

    print "gene pairs in ChiTars Only: %d" % chitars_only
    print "gene pairs in COSMIC Only: %d" % (total_cosmic - shared)
    print "Shared gene pairs: %d" % shared

    return cosmic_data


# Parse GenBank record and returns all pubmed ids
def parse_pubmed_id(genbank_record):
    new_ids = []
    for line in genbank_record:
        if re.search("PUBMED", line):
            new_ids.append(line.split()[1])
    return new_ids


# Query NCBI to get pubmed id from accession number
def lookup_accession(number):
    pubmed_ids = []

    try:
        handle = Entrez.efetch(db="nucleotide", rettype="gb", id=number)
        pubmed_ids = parse_pubmed_id(handle.read().split("\n"))

    except urllib2.HTTPError:
        print str(number) + " is not valid"
    except IOError:
        print "Problem connecting to NCBI"

    return pubmed_ids


# Check that id listed is a Pubmed id if not look up accession id and swap fields
def fix_pubmed_ids(ids):
    pubmed_ids = []
    seq_ids = []
    for pubmed_id in ids.split("_"):
        if re.search("[A-Z][0-9]", pubmed_id):
            seq_ids.append(pubmed_id)
            new_pubmed_id = lookup_accession(pubmed_id)
            if new_pubmed_id:
                pubmed_ids += new_pubmed_id
        else:
            pubmed_ids.append(pubmed_id)

    return set(list(pubmed_ids)), set(list(seq_ids))


# Parse all data exported from ChiTars
def parse_chitars(lines):
    breakpoints = my_dd()  # nested defaultdict
    for line in lines:
        if re.match("#", line):
            continue

        fields = line.split("\t")
        pubmed_ids = fields[0]
        breakpoint = fields[1]
        sequence_ids = fields[2].split("_")
        gene_pair = fields[3]+"-"+fields[4]

        # If a pubmed id looks like an accession number then fix the ids
        if re.search("[A-Z]", pubmed_ids):
            pubmed_ids, new_sequence_ids = fix_pubmed_ids(pubmed_ids)
            sequence_ids += new_sequence_ids
        else:
            pubmed_ids = pubmed_ids.split("_")

        breakpoints[gene_pair][breakpoint]["references"]["pubmed"] = list(set(pubmed_ids))
        breakpoints[gene_pair][breakpoint]["references"]["sequence"] = list(set(sequence_ids))
        breakpoints[gene_pair][breakpoint]["disease"] = fields[-1]

    return breakpoints


# Parse the Cancer breakpoints file downloaded from ChiTars
def parse_chitars_breakpoint(breakpoints_lines, database):

    breakpoints = my_dd()  # nested defaultdict
    breakpoint_entry = None
    breakpoint = None

    conn = lite.connect(database)

    conn.row_factory = lite.Row
    cur = conn.cursor()
    cur.execute("SELECT NAME, RESOURCE_ID FROM Resources")

    resource_id_lookup = {}
    for row in cur.fetchall():
        resource_id_lookup[row['NAME']] = row['RESOURCE_ID']


    for line in breakpoints_lines:
        lines = line.split('\r')  # Handle conversion from excel tab-del format

        for split_line in lines[1:]:  # Skip header line
            fields = split_line.split('\t')
            # 0 : Pubmed Reference id
            # 1 : Database providing pubmed id
            # 2 : GenBank Sequence id
            # 3 : Database providing sequence id
            # 4 : Breakpoint
            # 5 : Gene 1 id
            # 6 : Gene 2 id
            # 7 : Diseases

            # Check if line has a new breakpoint
            if fields[4]:

                breakpoint = fields[4]
                gene_pair = fields[5]+":"+fields[6]
                databases = set()
                for db in (fields[1], fields[3]):
                    if db:
                        databases.add(db)
                        conn.execute('INSERT INTO Evidence VALUES (?,?,NULL,NULL)',
                                    (resource_id_lookup[db], str(breakpoint)))

                conn.execute('INSERT INTO Known_Fusions VALUES (?,?,?)', (breakpoint, gene_pair, fields[7]))
                if fields[0]:
                    conn.execute('INSERT INTO Evidence VALUES (?,?,?,NULL)',
                                (resource_id_lookup["PubMed"], str(breakpoint), str(fields[0])))
                if fields[2]:
                    conn.execute('INSERT INTO Evidence VALUES (?,?,?,NULL)',
                                (resource_id_lookup["GenBank"], str(breakpoint), str(fields[2])))

                breakpoint_entry = breakpoints[gene_pair][breakpoint]
                breakpoint_entry["references"]["pubmed"] = {fields[0]}
                breakpoint_entry["references"]["sequence"] = {fields[2]}
                breakpoint_entry["disease"] = fields[7]
                breakpoint_entry["databases"] = databases

            # If not - add information to current breakpoint
            else:
                # Store pubmed id
                if fields[0]:
                    pubmed_ids = {fields[0]}
                    sequence_ids = None

                    #print breakpoint
                    #print resource_id_lookup[fields[1]]
                    if not fields[1] in resource_id_lookup:
                        print str(fields[1]) + " not in resource_id_lookup"

                    # If a pubmed id looks like an accession number then fix the ids
                    if re.search("[A-Z]", str(pubmed_ids)):
                        pubmed_ids, sequence_ids = fix_pubmed_ids(str(pubmed_ids.pop()))
                    if pubmed_ids:
                        for pubmed_id in pubmed_ids:
                            conn.execute('INSERT INTO Evidence VALUES (?,?,?,NULL)',
                                        (resource_id_lookup["PubMed"], str(breakpoint), str(pubmed_id)))
                        breakpoint_entry["references"]["pubmed"] = breakpoint_entry["references"]["pubmed"].union(pubmed_ids)
                    if sequence_ids:
                        for sequence_id in sequence_ids:
                            conn.execute('INSERT INTO Evidence VALUES (?,?,?,NULL)',
                                        (resource_id_lookup["GenBank"], str(breakpoint), str(sequence_id)))
                        breakpoint_entry["references"]["sequence"] = breakpoint_entry["references"]["sequence"].union(sequence_ids)

                    breakpoint_entry["databases"].add(fields[1])


                    cur.execute("SELECT * FROM Evidence WHERE RESOURCE_ID = ? AND FUSION_ID = ?",
                                (resource_id_lookup[fields[1]], str(breakpoint)))

                    if len(cur.fetchall()) == 0:
                        conn.execute('INSERT INTO Evidence VALUES (?,?,NULL,NULL)',
                                    (resource_id_lookup[fields[1]], str(breakpoint)))

                # Store sequence id
                if fields[2]:
                    breakpoint_entry["references"]["sequence"].add(fields[2])
                    breakpoint_entry["databases"].add(fields[3])

                    conn.execute('INSERT INTO Evidence VALUES (?,?,?,NULL)',
                                (resource_id_lookup["GenBank"], str(breakpoint), str(fields[2])))

                    cur.execute("SELECT * FROM Evidence WHERE RESOURCE_ID = ? AND FUSION_ID = ?",
                                (resource_id_lookup[fields[3]], str(breakpoint)))
                    if len(cur.fetchall()) == 0:
                        conn.execute('INSERT INTO Evidence VALUES (?,?,NULL,NULL)',
                                    (resource_id_lookup[fields[3]], str(breakpoint)))

    conn.commit()
    conn.close()
    return breakpoints


def parse_gene_info(data_lines):
    genes = []
    gene_list = []
    for line in data_lines:
        if re.search("GENES", line):
            if gene_list:
                genes.append(gene_list)
            gene_list = [line.split('\t')[1]]

        if re.search("ANNOTATION", line):
            gene_list.append(line.split('\t')[1])

    genes.append(gene_list)

   # for gene in genes:
   #     print gene

    return genes


def lookup_fusion(fusion_info, fusion_database):
    gene_pair = fusion_info[0]
    annotations = fusion_info[1:]

    orig_gene_pair = gene_pair

    if gene_pair in fusion_database:
        print "%s: Found" % gene_pair
        for annotation in annotations:
            print annotation
        print_gene_pair_entry(fusion_database[gene_pair])

    else:
        tmp = gene_pair.split(":")
        gene_pair = tmp[1] + ":" + tmp[0]

        if gene_pair in fusion_database:
            print "%s: Found as swap" % orig_gene_pair
            for annotation in annotations:
                print annotation
            print_gene_pair_entry(fusion_database[gene_pair])
        else:
            print "%s: Not Found" % orig_gene_pair


def print_gene_pair_entry(entry):
    ref_list = None
    for fusion_id in entry:
        print "\t" + fusion_id
        for key in entry[fusion_id]:
            if key == "references":
                print "\t\tReferences"
                for ref_type in entry[fusion_id]["references"]:

                    ref_list = ','.join(entry[fusion_id][key][ref_type])
                    print "\t\t\t%s: %s" % (ref_type, ref_list)
                    if ref_type == "pubmed":
                        print "\t\t\tLink: http://www.ncbi.nlm.nih.gov/pubmed/%s?" % ref_list
                    if ref_type == "sequence":
                        print "\t\t\tLink: http://www.ncbi.nlm.nih.gov/nuccore/%s?" % ref_list
                    if ref_type == "cosmic_id":
                        for cosmic_id in entry[fusion_id][key][ref_type]:
                            print "\t\t\tLink: http://cancer.sanger.ac.uk/cosmic/fusion/summary?id=%s" % cosmic_id
                    if ref_type == "cosmic_gene_ids":
                        print "\t\t\tLink: http://cancer.sanger.ac.uk/cosmic/fusion/overview?fid=%s&gid=%s" % \
                              (entry[fusion_id][key][ref_type][0], entry[fusion_id][key][ref_type][1])

            elif key == "transcript_details":
                print "\t\tTranscript Details"
                for transcript in entry[fusion_id]["transcript_details"]:
                    print "\t\t\t%s" % str(transcript)
            elif key == "databases":
                print "\t\t%s: %s" % (key, ','.join(entry[fusion_id][key]))

            else:
                print "\t\t%s: %s" % (key, entry[fusion_id][key])

# Link to NCBI queries
#http://www.ncbi.nlm.nih.gov/nuccore/<ids>?
#http://www.ncbi.nlm.nih.gov/pubmed/<ids>?

# Link to exact fusion
#http://cancer.sanger.ac.uk/cosmic/fusion/summary?id=
# Link to all fusion btw gene pair
#http://cancer.sanger.ac.uk/cosmic/fusion/overview?fid=<gene1id>&gid=<gene2_id>

# Link to all breakpoints btw gene pair
#http://chitars.bioinfo.cnio.es/cgi-bin/breakpoints.pl?refdis=3&searchstr=<gene1_name>%20%26%20<gene2_name>


def create_database(db_name):

    conn = None

    try:
        conn = lite.connect(db_name)
        #c = conn.cursor()

        # Create Resource Table
        conn.execute('''CREATE TABLE Resources
           (RESOURCE_ID INTEGER PRIMARY KEY,
           NAME         TEXT    NOT NULL,
           RANK         INT    NOT NULL,
           URL          TEXT,
           CAN_QUERY    INT);''')

        # Create Known Fusion Evidence table
        conn.execute('''CREATE TABLE Evidence
           (RESOURCE_ID    INT    NOT NULL,
           FUSION_ID       TEXT   NOT NULL,
           ACCESSION       TEXT,
           SAMPLE_COUNT    INT);''')

        # Create Known Fusions table
        conn.execute('''CREATE TABLE Known_Fusions
           (FUSION_ID     TEXT    NOT NULL,
           GENE_ORDER     TEXT,
           DISEASE        TEXT);''')

        # Create Known Fusion Details
        conn.execute('''CREATE TABLE Known_Fusions_Details
           (FUSION_ID             TEXT  PRIMARY KEY NOT NULL,
           FUSION_ORDER           INT,
           GENE_NAME              TEXT,
           GENOME_POS             INT,
           EXON                   TEXT,
           INTRON                 INT,
           CYTOBAND               TEXT,
           OMIM_ID                INT,
           LOCUS_ID               INT,
           ENSEMBL_TRANSCRIPT_ID  TEXT);''')

        # Load Resources
        resources = [
            ('COSMIC', 1, 'http://cancer.sanger.ac.uk/cosmic/fusion/summary?id=', 1),
            ('ChiTaRS', 2, 'http://chitars.bioinfo.cnio.es/cgi-bin/breakpoints.pl?'
                              'refdis=3&searchstr=<gene1_name>%20%26%20<gene2_name>', 1),
            ('TICdb', 3, 'http://www.unav.es/genetica/TICdb/', 1),
            ('ChimerDB', 4, 'http://biome.ewha.ac.kr:8080/FusionGene/', 1),
            ('dbCRID', 5, 'http://dbcrid.biolead.org/records.php', 1),
            ('Mitelman', 6, 'http://cgap.nci.nih.gov/Chromosomes/Mitelman', 0),
            ('PubMed', 7, 'http://www.ncbi.nlm.nih.gov/pubmed/<ids>?', 1),
            ('GenBank', 8, 'http://www.ncbi.nlm.nih.gov/nuccore/<ids>?', 1),
            ('OMIM', 9, 'http://www.omim.org/', 1),
            ('NCBI_Genes', 10, 'http://www.ncbi.nlm.nih.gov/gene/<ids>?', 1),
        ]

        conn.executemany('INSERT INTO Resources VALUES (NULL,?,?,?,?)', resources)

        conn.commit()

    except lite.Error, e:

        if conn:
            conn.rollback()

        print "Error %s:" % e.args[0]
        sys.exit(1)

    finally:

        if conn:
            conn.close()


def main(args):
    opts = parse_cmdline_params(args[1:])
    merged_data = None

    if opts.database_name:
        if not os.path.isfile(opts.database_name):
            create_database(opts.database_name)

    conn = lite.connect(opts.database_name)

    conn.row_factory = lite.Row
    cur = conn.cursor()

    cur.execute("SELECT NAME, RESOURCE_ID FROM Resources ORDER BY RANK")

    rows = cur.fetchall()

    print rows

    cur.execute("SELECT * FROM Resources ORDER BY RANK")

    for row in cur.fetchall():

        print row['RANK'], row['RESOURCE_ID'], row['NAME'], row['URL']

    cur.execute("SELECT * FROM Evidence")

    for row in cur.fetchall():

        print row

    if opts.cosmic_file:

        transcript_coords = defaultdict

        if opts.ensembl_file:
            transcript_coords = parse_transcript_gtf(opts.ensembl_file)
            print "Loaded GTF file"

        cosmic_lines = [f.rstrip() for f in opts.cosmic_file][1:]
        cosmic_gene_fusions = parse_cosmic_basic(cosmic_lines, transcript_coords)
        print "COSMIC processed"

    if opts.chitars_file:
        chitars_lines = [f.rstrip() for f in opts.chitars_file]
        chitars_breakpoints = parse_chitars_breakpoint(chitars_lines, opts.database_name)
        #chitars_breakpoints = parse_chitars(chitars_lines)
        print "CHiTars processed"

    if opts.cosmic_file and opts.chitars_file:
        merged_data = merge_data(cosmic_gene_fusions, chitars_breakpoints)
        print "Merged"

        pickle.dump(merged_data, open("merged_fusion_data.pickle", "wb"))
        print "Pickled"

    if opts.pickled_data:
        merged_data = pickle.load(open(opts.pickled_data, "rb"))

    if opts.test_data:
        test_lines = [f.rstrip() for f in opts.test_data]
        gene_info = parse_gene_info(test_lines)
        for fusion_info in gene_info:
            lookup_fusion(fusion_info, merged_data)

    print_data = False
    if merged_data and print_data:
        for gene_pair in merged_data:
            print gene_pair
            print_gene_pair_entry(merged_data[gene_pair])

    return 1


if __name__ == "__main__":
    main(sys.argv)
