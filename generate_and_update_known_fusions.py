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
import sqlite3 as lite

Entrez.email = "aberlin@enzymatics.com"


# Parse command line arguments
def parse_cmdline_params(arg_list=None):
    """Parses commandline arguments.
    """

    description = "Compare different fusion DBs."

    #Create instance of ArgumentParser
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument("-d",
                        "--database_name",
                        help="Database Name",
                        required=True)

    parser.add_argument("-c",
                        "--cosmic_file",
                        type=argparse.FileType('r'),
                        help="COSMIC fusion file",
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

    parser.add_argument("-a",
                        "--gene_aliases",
                        type=argparse.FileType('r'),
                        help="HGNC dataset file",
                        required=False)

    parser.add_argument("-t",
                        "--test_data",
                        type=argparse.FileType('r'),
                        help="Test data",
                        required=False)

    # Parse options
    opts = parser.parse_args(args=arg_list)

    return opts


### Functions to parse COSMIC fusion ID
def __parse_fusion(fusion, return_type="All"):
    """Private function to parse the COSMIC fusion in the hgvs format"""

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


def load_gene_aliases(gene_alias_lines, database):

    print "Loading Gene Alias"
    conn = lite.connect(database)
    for line in gene_alias_lines:
        # ignore withdrawn symbols
        if re.search("withdrawn", line):
            continue
        fields = line.split("\t")
        # ignore symbols without aliases
        if not fields[6] and not fields[8]:
            continue
        approved_name = fields[1]

        for alias in ", ".join([fields[6], fields[8]]).split(", "):
            if alias:  # skip blank entries
                conn.execute('INSERT INTO Gene_Name_Alias VALUES (?,?)',
                            (approved_name, alias))

    conn.commit()
    conn.close()


def get_all_from_fusion(fusion):
    """Return the all information available from the COSMIC fusion id"""
    # Returns a list of [gene_name, transcript, coordinates in gene] for each transcript in the fusion
    return __parse_fusion(fusion, "All")


def get_transcripts(fusion):
    """Return the transcripts from the COSMIC fusion id"""
    return __parse_fusion(fusion, "transcripts")


def get_gene_transcript_pairs(fusion):
    """Return the gene name and transcript pairs from the COSMIC fusion id"""
    return __parse_fusion(fusion, "pairs")


def get_genes(fusion):
    """Return the gene names from the COSMIC fusion id"""
    return __parse_fusion(fusion, "genes")


### Maintain data integrity for parsing external data
def get_exons_involved(breakpoints, transcript_details):
    """Looks up the exon closest to the breakpoint.  Can handle
     insertions, deletions and breakpoints in introns.
    """
    exon_list = []
    empty_list = ["Unknown", None, None, None]

    exon_coords = transcript_details["exons"]
    exon_positions = transcript_details["exon_positions"]

    for breakpoint_edge in breakpoints:
            intron_bases = 0
            indel = None

            # If we don't have an exon coordinate list for the transcript
            # All the ? should have already been converted into Null but just in case
            # Ignore Null coordinates
            if not exon_coords or breakpoint_edge == "?" or breakpoint_edge == "Null":
                exon_list.append(empty_list)
                continue

            # Not sure how to handle the () yet - Skip them for now
            if re.search("\(", breakpoint_edge):
                exon_list.append(empty_list)
                print "WARN: Coordinate w/ () is being skipped: " + breakpoint_edge
                continue

            # If the edge coordinate has a + or - then the break point is actually X bases in to the intron
            if re.search("\+", breakpoint_edge):
                breakpoint_edge, bp_into_intron = breakpoint_edge.split("+")
                intron_bases = bp_into_intron

            if re.search("\-", breakpoint_edge):
                breakpoint_edge, bp_into_intron = breakpoint_edge.split("-")
                intron_bases = "-" + bp_into_intron

            # Handle if the coordinate has an insertion or deletion present
            if re.search("ins|del", breakpoint_edge):
                breakpoint_edge, indel = breakpoint_edge.split('_')

            # If breakpoint edge coordinate has some weird text lets grab it and skip it for now
            try:
                breakpoint_edge = int(breakpoint_edge)
            except ValueError:
                print "WARN: Weird breakpoint is being skipped: " + breakpoint_edge + " " + str(breakpoints)
                exon_list.append(empty_list)
                continue

            for j in range(len(exon_coords) - 1):

                #Handle insertions
                if breakpoint_edge > int(exon_coords[len(exon_coords) - 1]) and indel:
                    exon_list.append([str(len(exon_coords)-1), 1, intron_bases, indel])
                    break

                #Handle offset
                if breakpoint_edge > int(exon_coords[len(exon_coords) - 1]):
                    exon_list.append([">"+str(len(exon_coords)-1), 0, intron_bases, indel])
                    break

                #find pair of exon coordinates that contain breakpoint
                if int(exon_coords[j]) < breakpoint_edge <= int(exon_coords[j+1]):
                    # if breakpoint matches an exon edge then assume it is not an exact breakpoint
                    if int(exon_coords[j]) + 1 == breakpoint_edge or breakpoint_edge == int(exon_coords[j+1]):
                        exact_breakpoint = 0
                    else:
                        exact_breakpoint = exon_positions[j] - 1 + (breakpoint_edge - exon_coords[j])
                    exon_list.append([str(j+1), exact_breakpoint, intron_bases, indel])
                    break

    return exon_list


def parse_pubmed_id(genbank_record):
    """Parse GenBank record and return all pubmed ids"""

    new_ids = []
    for line in genbank_record:
        if re.search("PUBMED", line):
            new_ids.append(line.split()[1])
    return new_ids


def lookup_accession(number):
    """Query NCBI to get pubmed id from accession number"""

    pubmed_ids = []
    try:
        handle = Entrez.efetch(db="nucleotide", rettype="gb", id=number)
        pubmed_ids = parse_pubmed_id(handle.read().split("\n"))

    except urllib2.HTTPError:
        print "WARN: " + str(number) + " is not valid"
    except IOError:
        print "Problem connecting to NCBI"

    return pubmed_ids


def fix_pubmed_ids(ids):
    """Check that id listed is a Pubmed id if not look up accession id and swap fields"""
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

    return pubmed_ids, seq_ids


def load_pubmed_sequence_evidence(breakpoint, pubmed_id, sequence_id, database_handle):

    resource_id_lookup = get_resources_lookup(database_handle)
    #pubmed_ids = pubmed_id
    sequence_ids = [sequence_id]
    if re.search("[A-Z]", str(pubmed_id)):
        pubmed_ids, new_sequence_ids = fix_pubmed_ids(pubmed_id)
        sequence_ids += new_sequence_ids
    else:
        pubmed_ids = [pubmed_id]

    if pubmed_ids:
        for pubmed_id in pubmed_ids:
            database_handle.execute('INSERT INTO Evidence VALUES (?,?,?,NULL)',
                                   (resource_id_lookup["PubMed"], str(breakpoint), str(pubmed_id)))
    if sequence_ids:
        for sequence_id in sequence_ids:
            if sequence_id:
                database_handle.execute('INSERT INTO Evidence VALUES (?,?,?,NULL)',
                                       (resource_id_lookup["GenBank"], str(breakpoint), str(sequence_id)))


def split_cytoband(cytoband):
    #Example format: t(11;22)(p13;q22)

    gene_cytobands = []
    cytoband_parts = cytoband.split("(")

    if len(cytoband_parts) < 3:
        print "WARN Cytoband does not look valid: " + cytoband
        return []

    prefix = cytoband_parts[0]
    chromosomes = cytoband_parts[1].rstrip(")").split(";")
    arms = cytoband_parts[2].rstrip(")").split(";")

    if not len(chromosomes) == len(arms):
        print "WARN: Number of Chromosomes does not match Number of Arms: " + cytoband
        return []

    for i, chromosome in enumerate(chromosomes):
        gene_cytobands.append(chromosome + arms[i])

    return gene_cytobands


def valid_fusion(fusion_id):

    fusion_details = get_all_from_fusion(fusion_id)

    for gene_name, transcript, transcript_range in fusion_details['full_data']:
        coords = transcript_range.split("_")

        # Check for _[ins|del]* or other wierdness
        if len(coords) >= 3:
            if not re.search("ins|del", coords[2]):
                #print "WARN: New coordinate format being skipped: " + transcript_range
                return False

        # Not sure how to handle the () yet - Skip them for now
        if re.search("\(", transcript_range):
            return False

        # Check if range values are all ints
        for breakpoint_edge in coords:
            # Fusion coordinates that are "?" or ins/del are ok and handled properly
            if breakpoint_edge == "?" or re.search("ins|del", breakpoint_edge):
                continue

            # If the edge coordinate has a + or - then the break point is actually X bases in to the intron
            if re.search("\+", breakpoint_edge):
                breakpoint_edge, bp_into_intron = breakpoint_edge.split("+")
            if re.search("\-", breakpoint_edge):
                breakpoint_edge, bp_into_intron = breakpoint_edge.split("-")

            try:
                int(breakpoint_edge)
            except ValueError:
                #print "WARN: Weird breakpoint is being skipped: " + breakpoint_edge + " " + str(coords)
                return False

    else:
        return True


### Functions to lookup transcript and gene information from gtf and Entrez queries
def get_gene_record(record_id):
    """Queries NCBI by gene ID for a specific gene record
    """
    gene_record = None
    try:
        gene_handle = Entrez.efetch(db="gene", rettype="xml", id=record_id)
        gene_record = Entrez.read(gene_handle)[0]
    except urllib2.HTTPError:
        print "WARN: " + str(record_id) + " is not valid"
        return None
    except IOError:
        print "Problem connecting to NCBI"
        return None
    finally:
        return gene_record


def parse_gene_record(gene_record):
    """Parses the gene record from Entrez for specific fields we want
    """

    omim_id = None
    locus_id = None
    ensembl_id = None

    if not gene_record or not gene_record["Entrezgene_locus"][0]:
        return None
    if not "Gene-commentary_seqs" in gene_record["Entrezgene_locus"][0]:
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
    if "Gene-ref_db" in gene_record["Entrezgene_gene"]["Gene-ref"]:
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


def run_search(query):
    """Run the search on Entrez - NCBI
    """
    search_handle = Entrez.esearch(db="gene", term=query)
    record = Entrez.read(search_handle)

    if int(record["Count"]) == 1:

        return parse_gene_record(get_gene_record(record['IdList'][0]))

    # Warn if no records are found
    elif int(record["Count"]) == 0:
        #print record
        print "WARN: No genes found, Search Attempted: " + str(record["QueryTranslation"])
        return None

    # Warn if more than one record is found
    else:
        # Search the multiple records returned for one with the same gene name not gene alias
        target_gene = query.split("[")[0]
        for gene_id in record["IdList"]:
            gene_record = get_gene_record(gene_id)
            if target_gene == gene_record["Entrezgene_gene"]["Gene-ref"]["Gene-ref_locus"]:
                return parse_gene_record(gene_record)
        print "WARN: More than one gene ID found: " + record["Count"] + " records found" + \
              ", Search Attempted: " + str(record["QueryTranslation"])
        return None


def query_by_gene(gene_name):
    """Query Entrez by Gene Name"""

    return run_search(gene_name + "[gene name] AND Homo[organism] AND alive[prop]")


def query_by_transcript(transcript, gene_name, ensembl_transcript_details):
    """Lookup transcript information from Ensembl gtf and NCBI queries"""

    # remove the version from accession numbers for better search results
    if re.search('\.\d+$', transcript):
        transcript = re.sub('\.\d+$', "", transcript)

    # Search for transcript and ignore discontinued records
    ncbi_info = run_search(transcript + " AND alive[prop]")

    # If the transcript is not found lookup the gene name
    if not ncbi_info:
        print "\tQuerying by gene name: " + gene_name,
        ncbi_info = query_by_gene(gene_name)
        if not ncbi_info:
            print "- No records found"
        else:
            print " - Found"

    ensembl_info = None
    exon_info = []

    if transcript in ensembl_transcript_details:
        ensembl_info = ensembl_transcript_details[transcript]["details"]
        exon_info = ensembl_transcript_details[transcript]["exons"]

    if not ncbi_info and not ensembl_info:
        return None, None, None

    if not ensembl_info:
        return ncbi_info, None, None

    if not ncbi_info:
        ncbi_info = ["", "", "", ""]

    return ncbi_info[:4] + ensembl_info, exon_info, ensembl_transcript_details[transcript]["exon_positions"]


def parse_transcript_gtf(transcript_file):
    """Parse the Ensembl gtf annotations for transcript coordinates"""
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

            transcript_details[trans_id] = {}
            # exon position on the transcript
            transcript_details[trans_id]["exons"] = [0]
            # exon position on the genome
            transcript_details[trans_id]["exon_positions"] = []
            # details about the transcript
            transcript_details[trans_id]["details"] = [gene_id, chromosome, start_coord, end_coord, strand]
            running_total = 0

        if fields[2] == "exon":

            for entry in comments:
                if re.search("transcript_id", entry):
                    name, trans_id = re.sub("\"", "", entry).split()

            running_total += int(end_coord) - int(start_coord) + 1

            # Store the coordinate of the the exon end on the transcript
            try:
                transcript_details[trans_id]["exons"].append(running_total)
                transcript_details[trans_id]["exon_positions"].append(int(start_coord))
            except KeyError:
                print "ERROR: transcript file isn't valid"
                sys.exit(-1)

    return transcript_details


### Parse COSMIC database and load into internal DB

def parse_cosmic_basic(cosmic_lines, ensembl_transcript_details, database):
    """Parse genes and PubMed reference from CosmicFusionExport_vXX.tsv file"""
    transcript_details = {}
    seen_fusions = set()

    conn = lite.connect(database)
    resource_id_lookup = get_resources_lookup(conn)

    for line in cosmic_lines:
        #Skip lines that are not fusions
        if not re.search("{", line):
                continue
        fields = re.split("\t", line)

        if len(fields) >= 8:
            cosmic_id = fields[7]
            fusion_id = fields[8]
            disease = "{}-{}".format(fields[3], fields[4])

            # Store the Pubmed reference id if present
            if len(fields) >= 12:
                reference = fields[12]
            else:
                reference = None

            ## Catch and skip invalid fusion IDs
            if not valid_fusion(fusion_id):
                ## Try to correct fusion containing ()
                edited_id = re.sub("\(|\)", "", fusion_id)
                if valid_fusion(edited_id):
                    fusion_id = edited_id
                else:
                    print "WARN: Fusions is not valid and will be skipped: {}, {}".format(fusion_id, edited_id)
                    continue

            fusion_details = get_all_from_fusion(fusion_id)
            gene_names = fusion_details['gene_names']
            gene_pair = ":".join(gene_names)

            # Store Fusion and supporting evidence
            conn.execute('INSERT INTO Known_Fusions VALUES (?,?,?)', (fusion_id, gene_pair, disease))
            if reference:
                conn.execute('INSERT INTO Evidence VALUES (?,?,?,NULL)',
                            (resource_id_lookup["PubMed"], str(fusion_id), str(reference)))
            conn.execute('INSERT INTO Evidence VALUES (?,?,?,?)',
                        (resource_id_lookup["COSMIC"], str(fusion_id), str(cosmic_id), 0))
            conn.execute('UPDATE Evidence '
                         'SET SampleCount = SampleCount + 1 '
                         'WHERE ResourceID = ? AND FusionID = ? AND Accession = ?',
                         (resource_id_lookup["COSMIC"], str(fusion_id), str(cosmic_id)))

            # Add the gene and exon information for new fusions
            if not fusion_id in seen_fusions:

                counter = 0
                genes_involved = len(gene_names)

                # Load all the gene and exon data
                for gene_name, transcript, transcript_range in fusion_details['full_data']:
                    if not transcript in transcript_details:

                        # Query NCBI for transcripts details
                        details, exon_details, exon_positions = query_by_transcript(transcript, gene_name, ensembl_transcript_details)

                        # Save the transcript details so we don't have to query NCBI again
                        transcript_details[transcript] = {}
                        transcript_details[transcript]["details"] = details
                        # 0: omim_id,
                        # 1: locus_id,
                        # 2: cytoband_location,
                        # 3: gene_ref_accession,
                        # 4: ensembl_id,
                        # 5: gene_chromosome,
                        # 6: gene_start_coord,
                        # 7: gene_end_coord,
                        # 8: gene_strand

                        transcript_details[transcript]["exons"] = exon_details
                        transcript_details[transcript]["exon_positions"] = exon_positions

                    conn.execute('INSERT INTO Gene_Evidence VALUES (?,?,?)',
                                (resource_id_lookup["OMIM"], str(gene_name),
                                 str(transcript_details[transcript]["details"][0])))
                    conn.execute('INSERT INTO Gene_Evidence VALUES (?,?,?)',
                                (resource_id_lookup["NCBI_Genes"], str(gene_name),
                                 str(transcript_details[transcript]["details"][1])))

                    details = transcript_details[transcript]["details"]

                    if details[8] == "+" or details[8] == "plus":
                        strand = 1
                    else:
                        strand = -1

                    conn.execute('INSERT INTO Gene_Details VALUES (?,?,?,?,?,?)',
                                (gene_name, transcript_details[transcript]["details"][5],
                                 transcript_details[transcript]["details"][6],
                                 transcript_details[transcript]["details"][7],
                                 strand,
                                 str(transcript_details[transcript]["details"][2])
                                 ))

                    coords = transcript_range.split("_")

                    if coords[0] == "?":
                        coords = ["Null", "Null"]

                    # Check for _[ins|del]* or other weirdness
                    if len(coords) == 3:
                        coords[1] += "_" + coords[2]
                        if not re.search("ins|del", coords[2]):
                            print "WARN: New coordinate format being skipped: " + transcript_range
                            continue

                    # Only load the 5' side of the fusion
                    if counter == 0:
                        exons = get_exons_involved([coords[1]], transcript_details[transcript])[0]

                        conn.execute('INSERT INTO Fusion_Part_Details VALUES '
                                     '(?,?,?,?,?,?,?,?,?,?,?,?,?,NULL,NULL,NULL,NULL)',
                                    (fusion_id, counter, gene_name, details[5], details[6],
                                     details[7], strand, transcript,
                                     coords[0], coords[1], exons[0], exons[1], exons[2]))

                    # Only load the 3' side of the fusion
                    elif counter == genes_involved - 1:
                        exons = get_exons_involved([coords[0]], transcript_details[transcript])[0]

                        conn.execute('INSERT INTO Fusion_Part_Details VALUES '
                                     '(?,?,?,?,?,?,?,?,?,?,NULL,NULL,NULL,?,?,?,NULL)',
                                    (fusion_id, counter, gene_name, details[5], details[6],
                                     details[7], strand, transcript,
                                     coords[0], coords[1], exons[0], exons[1], exons[2]))
                    # Load both sides of the fusion
                    else:
                        exons5, exons3 = get_exons_involved([coords[0], coords[1]],
                                                            transcript_details[transcript])

                        conn.execute('INSERT INTO Fusion_Part_Details VALUES '
                                     '(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,NULL)',
                                    (fusion_id, counter, gene_name, details[5], details[6],
                                     details[7], strand, transcript, coords[0], coords[1],
                                     exons5[0], exons5[1], exons5[2], exons3[0], exons3[1], exons3[2]))
                    counter += 1

            seen_fusions.add(fusion_id)

    conn.commit()
    conn.close()


def parse_chitars_breakpoint(breakpoints_lines, database):
    """Parse the Cancer breakpoints file downloaded from ChiTars"""
    breakpoint = None
    breakpoint_valid = True
    failed_gene_search = set()
    conn = lite.connect(database)
    resource_id_lookup = get_resources_lookup(conn)

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

                pubmed_ids = fields[0]
                sequence_ids = fields[2]
                breakpoint = fields[4]
                gene_pair = fields[5]+":"+fields[6]
                gene_list = [fields[5], fields[6]]
                breakpoint_valid = True

                # Split cytoband and load into Fusion Parts
                cytoband_parts = split_cytoband(breakpoint)
                if not cytoband_parts:
                    breakpoint_valid = False
                    print "WARN: ChiTars breakpoint is invalid skipping : " + breakpoint
                    continue

                for i, loci in enumerate(cytoband_parts):
                    chromosome = loci.split("(")[0]
                    if i >= 2:
                        gene_name = "Unknown"
                    else:
                        gene_name = gene_list[i]
                    conn.execute('INSERT INTO Fusion_Part_Details VALUES '
                                 '(?,?,?,?,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,?)',
                                 (breakpoint, i, gene_name, chromosome, loci))

                # Load into database tables
                for db in (fields[1], fields[3]):
                    if db:
                        conn.execute('INSERT INTO Evidence VALUES (?,?,?,NULL)',
                                    (resource_id_lookup[db], str(breakpoint), "db"))

                conn.execute('INSERT INTO Known_Fusions VALUES (?,?,?)', (breakpoint, gene_pair, fields[7]))

                # Load Pubmed and sequence evidence - check validity
                load_pubmed_sequence_evidence(breakpoint, pubmed_ids, sequence_ids, conn)

                # If Gene is not in the Gene_Details table lookup info from NCBI and populate
                for gene_name in gene_list:
                    cursor = conn.cursor()
                    cursor.execute('SELECT rowid FROM Gene_Details WHERE GeneName = ?', (gene_name,))
                    data = cursor.fetchone()
                    if not data and not gene_name in failed_gene_search:
                        gene_details = query_by_gene(gene_name)
                            # 0: omim_id,
                            # 1: locus_id,
                            # 2: cytoband_location,
                            # 3: gene_ref_accession,
                            # 4: ensembl_id,
                            # 5: gene_chromosome,
                            # 6: gene_start_coord,
                            # 7: gene_end_coord,
                            # 8: gene_strand
                        if not gene_details:
                            failed_gene_search.add(gene_name)
                            continue

                        if gene_details[8] == "+" or gene_details[8] == "plus":
                            strand = 1
                        else:
                            strand = -1

                        conn.execute('INSERT INTO Gene_Details VALUES (?,?,?,?,?,?)',
                                    (gene_name, gene_details[5], gene_details[6],
                                     gene_details[7], strand, str(gene_details[2])))
                        conn.execute('INSERT INTO Gene_Evidence VALUES (?,?,?)',
                                    (resource_id_lookup["OMIM"], str(gene_name), str(gene_details[0])))
                        conn.execute('INSERT INTO Gene_Evidence VALUES (?,?,?)',
                                    (resource_id_lookup["NCBI_Genes"], str(gene_name), str(gene_details[1])))

            # If not - add information to current breakpoint
            elif breakpoint_valid:
                # Store pubmed id and external database supporting breakpoint
                if fields[0]:
                    pubmed_ids = fields[0]
                    sequence_ids = None

                    if not fields[1] in resource_id_lookup:
                        print "WARN: " + str(fields[1]) + " not in the Resources Table"

                    load_pubmed_sequence_evidence(breakpoint, pubmed_ids, sequence_ids, conn)


                    conn.execute('INSERT INTO Evidence VALUES (?,?,?,NULL)',
                                (resource_id_lookup[fields[1]], str(breakpoint), "db"))

                # Store sequence id and external database supporting breakpoint
                if fields[2]:
                    if not fields[3] in resource_id_lookup:
                        print "WARN: " + str(fields[3]) + " not in the Resources Table"
                    conn.execute('INSERT INTO Evidence VALUES (?,?,?,NULL)',
                                (resource_id_lookup["GenBank"], str(breakpoint), str(fields[2])))
                    conn.execute('INSERT INTO Evidence VALUES (?,?,?,NULL)',
                                (resource_id_lookup[fields[3]], str(breakpoint), "db"))

    conn.commit()
    conn.close()
    #return breakpoints


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
        #print_gene_pair_entry(fusion_database[gene_pair])

    else:
        tmp = gene_pair.split(":")
        gene_pair = tmp[1] + ":" + tmp[0]

        if gene_pair in fusion_database:
            print "%s: Found as swap" % orig_gene_pair
            for annotation in annotations:
                print annotation
           # print_gene_pair_entry(fusion_database[gene_pair])
        else:
            print "%s: Not Found" % orig_gene_pair


# Link to NCBI queries
#http://www.ncbi.nlm.nih.gov/nuccore/<ids>?
#http://www.ncbi.nlm.nih.gov/pubmed/<ids>?

# Link to exact fusion
#http://cancer.sanger.ac.uk/cosmic/fusion/summary?id=
# Link to all fusion btw gene pair
#http://cancer.sanger.ac.uk/cosmic/fusion/overview?fid=<gene1id>&gid=<gene2_id>

# Link to all breakpoints btw gene pair
#http://chitars.bioinfo.cnio.es/cgi-bin/breakpoints.pl?refdis=3&searchstr=<gene1_name>%20%26%20<gene2_name>

### Functions for Creating and Querying the internal DB
def create_database(db_name):
    """ Create and setup a sqlite3 database for known fusion information
    """
    conn = None

    try:
        conn = lite.connect(db_name)

        # Create Resource Table
        conn.execute('''CREATE TABLE Resources
           (ResourceID   INTEGER   PRIMARY KEY,
           Name         TEXT      NOT NULL,
           Rank         INT       NOT NULL,
           Url          TEXT,
           CanQuery     INT
           );''')

        # Create Known Fusion Evidence table
        conn.execute('''CREATE TABLE Evidence
           (ResourceID     INT   NOT NULL,
           FusionID       TEXT   NOT NULL,
           Accession      TEXT,
           SampleCount    INT,
           PRIMARY KEY (ResourceID, FusionID, Accession) ON CONFLICT IGNORE
           );''')

        # Create Known Fusions table
        conn.execute('''CREATE TABLE Known_Fusions
           (FusionID       TEXT  PRIMARY KEY  ON CONFLICT IGNORE  NOT NULL,
           GeneOrder      TEXT,
           Disease        TEXT
           );''')

        # Create Fusion Part Details table
        conn.execute('''CREATE TABLE Fusion_Part_Details
           (FusionID             TEXT  NOT NULL,
           FusionOrder           INT   NOT NULL,
           GeneName              TEXT,
           Chromosome            TEXT, --Don't forget about X and Y
           GenomePosStart        INT,
           GenomePosEnd          INT,
           Strand                INT, -- +1 or -1
           TranscriptID          TEXT,
           TranscriptPosStart    INT,
           TranscriptPosEnd      INT,
           ExonFusion5prime      TEXT, --Closest exon to the 5' side of the fusion
           ExactBreakpoint5prime INT,  --0 if exact breakpoint unknown; otherwise genome position is listed
           IntronBases5prime     INT,  --Bases into the intron if breakpoint is not exonic
           ExonFusion3prime      TEXT, --Closest exon to the 3' side of the fusion
           ExactBreakpoint3prime INT,  --0 if exact breakpoint unknown; otherwise genome position is listed
           IntronBases3prime     INT,  --Bases into the intron if breakpoint is not exonic
           Cytoband              TEXT,
           PRIMARY KEY (FusionID, FusionOrder) ON CONFLICT IGNORE
           );''')

        # Create Gene Details table
        conn.execute('''CREATE TABLE Gene_Details
            (GeneName            TEXT   PRIMARY KEY  ON CONFLICT IGNORE  NOT NULL,
            Chromosome           TEXT, -- Don't forget about X and Y
            GenomePosStart       INT,
            GenomePosEnd         INT,
            Strand               INT,  -- +1 or -1
            Cytoband             TEXT
            );''')

        # Create Gene Evidence table
        conn.execute('''CREATE TABLE Gene_Evidence
           (ResourceID     INT   NOT NULL,
           GeneName       TEXT   NOT NULL,
           Accession      TEXT,
           PRIMARY KEY (ResourceID, GeneName, Accession) ON CONFLICT IGNORE
           );''')

        # Create Gene name alias table
        conn.execute('''CREATE TABLE Gene_Name_Alias
            (ApprovedGeneName            TEXT NOT NULL,
            GeneNameAlias                TEXT NOT NULL,
            PRIMARY KEY (ApprovedGeneName, GeneNameAlias) ON CONFLICT IGNORE
            );''')

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
            ('OMIM', 9, 'http://www.omim.org/entry/<id>#cytogenetics', 1),
            ('NCBI_Genes', 10, 'http://www.ncbi.nlm.nih.gov/gene/<ids>?', 1),
            ('Enzymatics', 11, 'http://archer.dev.enzymatics.com/', 0)
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


def get_resources_lookup(database_handle):
    """Query database and create a lookup dict for resource name and ID"""
    #conn = lite.connect(database)

    # Store resources by Name to lookup the ID
    database_handle.row_factory = lite.Row
    cur = database_handle.cursor()
    cur.execute("SELECT Name, ResourceID FROM Resources")
    resource_id_lookup = {}
    for row in cur.fetchall():
        resource_id_lookup[row['Name']] = row['ResourceID']

    return resource_id_lookup


def main(args):
    opts = parse_cmdline_params(args[1:])

    if opts.database_name:
        if not os.path.isfile(opts.database_name):
            create_database(opts.database_name)
        ###Confirm that Database is valid

    if opts.gene_aliases:
        alias_lines = [f.rstrip() for f in opts.gene_aliases][1:]
        load_gene_aliases(alias_lines, opts.database_name)

    if opts.cosmic_file:

        transcript_coords = defaultdict

        if opts.ensembl_file:
            transcript_coords = parse_transcript_gtf(opts.ensembl_file)
            print "Loaded GTF file"
            sys.stdout.flush()

        cosmic_lines = [f.rstrip() for f in opts.cosmic_file][1:]
        parse_cosmic_basic(cosmic_lines, transcript_coords, opts.database_name)
        print "COSMIC processed"
        sys.stdout.flush()

    if opts.chitars_file:
        chitars_lines = [f.rstrip() for f in opts.chitars_file]
        parse_chitars_breakpoint(chitars_lines, opts.database_name)
        print "CHiTars processed"
        sys.stdout.flush()

    if opts.test_data:
        test_lines = [f.rstrip() for f in opts.test_data]

        print test_lines
    else:
        test_lines = None

    conn = lite.connect(opts.database_name)

    conn.row_factory = lite.Row
    cur = conn.cursor()

    cur.execute("SELECT * FROM Evidence WHERE ResourceID= 7 ORDER BY Accession")
    #for row in cur.fetchall():
     #   print row

    cur.execute("SELECT * FROM Resources ORDER BY Rank")
    #for row in cur.fetchall():
        #print row['Rank'], row['ResourceID'], row['Name'], row['Url']

    cur.execute("SELECT * FROM Known_Fusions")
    #for row in cur.fetchall():
        #print row
    if test_lines:
        for gene in test_lines:
            print gene
            cur.execute("SELECT GeneOrder, Disease,  Evidence.Accession FROM Known_Fusions, Evidence \
                         WHERE Evidence.FusionID = Known_Fusions.FusionID AND GeneOrder LIKE '%"+gene+"%' \
                         AND Evidence.ResourceID='7' AND length(Evidence.Accession) >2")

            pairs = set()
            pubmed = set()
            disease = set()
            for row in cur.fetchall():
                pairs.add(row[0])
                pubmed.add(row[2])
                disease.add(row[1])

            for item in pairs:
                print item
            print ", ".join(pubmed)
            for item in disease:
                print item

    cur.execute("SELECT * FROM Fusion_Part_Details ORDER BY FusionID, FusionOrder")
    #for row in cur.fetchall():
        #print row

    cur.execute("SELECT * FROM Evidence ORDER BY ResourceID, FusionID, Accession")
    #for row in cur.fetchall():
        #print row['ResourceID'], row['FusionID'], row['Accession'], row['SampleCount']

    return 1

if __name__ == "__main__":
    main(sys.argv)
