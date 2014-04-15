#!/usr/bin/env python
"""Lookup Accession number and replace it with the corresponding PubMed id"""

__author__ = 'aberlin'

import argparse
import re
import sys
import urllib2
from Bio import Entrez


def parse_cmdline_params(arg_list=None):
    """Parses commandline arguments.
    """

    description = "Compare different fusion DBs."

    #Create instance of ArgumentParser
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument("-f",
                        "--file",
                        help="fusion file",
                        required=True)

    parser.add_argument("-o",
                        "--output_file",
                        help="output file",
                        required=True)

    # Parse options
    opts = parser.parse_args(args=arg_list)

    return opts


def parse_pubmed_id(genbank_record):
    for line in genbank_record:
        if re.search("PUBMED", line):
            return line.split()[1]
    else:
        return None


def lookup_accession(number):
    Entrez.email = "aberlin@enzymatics.com"
    handle = Entrez.efetch(db="nucleotide", rettype="gb", id=number)
    pubmed_id = parse_pubmed_id(handle.read().split("\n"))
    return pubmed_id


def lookup_transcript(transcript):
    Entrez.email = "aberlin@enzymatics.com"
    omim_id = None
    locus_id = None
    ensembl_id = None

    search_handle = Entrez.esearch(db="gene", term=transcript)
    record = Entrez.read(search_handle)

    if int(record["Count"]) == 1:
        try:
            gene_handle = Entrez.efetch(db="gene", rettype="xml", id=record["IdList"][0])
            gene_record = Entrez.read(gene_handle)[0]
        except urllib2.HTTPError:
            print str(record["IdList"][0]) + " is not valid"
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

        print omim_id, locus_id, str(ensembl_id)
        print cytoband_location, gene_ref_accession, gene_chromosome, gene_start_coord, gene_end_coord, gene_strand

    else:
        print "WARN: More than on gene ID found: "
        print record
        return None

def main(args):
    #opts = parse_cmdline_params(args[1:])

    lookup_transcript("NM_005656")

    exit()

    with open(opts.output_file, 'w') as output:
        with open(opts.file, 'r') as input_file:
            fusions = input_file.read().split("\n")
            for fusion in fusions:
                elements = fusion.rstrip().split("\t")
                print elements
                if re.search("[A-Z]", elements[2]):
                    pubmed_id = lookup_accession(elements[2])
                    if pubmed_id:
                        #print "Replacing %s with %s" % (elements[2], pubmed_id)
                        fusion = "%s\t%s\t%s\t%s" % (elements[0], elements[1], pubmed_id, elements[3])

                output.write(fusion.rstrip() + "\n")

    return 1


if __name__ == "__main__":
    main(sys.argv)
