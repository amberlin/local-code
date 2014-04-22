__author__ = 'aberlin'


def sort_sam_by_name(samfile):
    unsorted_bam = re.sub("\.sam", ".bam", samfile)
    sorted_bam = re.sub("\.sam", ".sorted.bam", samfile)
    sorted_sam = re.sub("\.bam", ".sam", sorted_bam)

    #convert SAM -> BAM
    command("samtools view -bS " + samfile + " > " + unsorted_bam)
    #Sort BAM by query name
    pysam.sort("-nf", unsorted_bam, sorted_bam)
    #convert BAM -> SAM
    command("samtools view -h " + sorted_bam + " > " + sorted_sam)

    return sorted_sam

def bed2fasta(bed_file, fasta_file):
    command("fastaFromBed -fi $ARCHER_GENOME -s -name -fo " + fasta_file + " -bed " + bed_file)


def create_bowtie_index(reference_fasta):
    command("bowtie2-build " + reference_fasta + " " + reference_fasta)


def align_with_bowtie(reads_fastq, reference_fasta, output):

    # Run Alignment
    command("bowtie2 -U " + reads_fastq + " -x " + reference_fasta + " --local --very-sensitive-local -S " + output)


def search_by_direct_alignment(opts, tmp_dir, header):
    GSP2_reference = tmp_dir + "GSP2_seqs.fasta"

    # Convert Bed to Fasta
    bed2fasta(opts.GSP2_bed, GSP2_reference)

    # Align reads to GSP2 sequences
    align_with_bowtie(opts.fastq, GSP2_reference, tmp_dir + header + ".sam", "-a")

    # Sort the aligned SAM by query name
    sorted_sam = sort_sam_by_name(tmp_dir + header + ".sam")
    #sorted_sam = (tmp_dir + header + ".sorted.sam")

    samfile = pysam.Samfile(sorted_sam, "r")

    GSP2_lengths = samfile.lengths
    GSP2_start = defaultdict(int)
    GSP2_read_though = defaultdict(int)
    GSP2_references = samfile.references

    read_name = None
    found_at_start = False
    potential_read_through = []

    for aligned_read in samfile.fetch():
        if not aligned_read.is_unmapped:
            full_length = False
            if not read_name or not read_name == aligned_read.qname:
                read_name = aligned_read.qname
                found_at_start = False
                potential_read_through = []

            GSP2_name = samfile.getrname(aligned_read.tid)

            # First element in the cigar string in a M and the length is the full length of the primer
            if aligned_read.cigar[0][0] == 0 and aligned_read.cigar[0][1] == GSP2_lengths[aligned_read.tid]:
                found_at_start = True
                full_length = True
                GSP2_start[GSP2_name] += 1
                #print "Start - " + GSP2_name
                if potential_read_through:
                    for GSP2 in potential_read_through:
                        GSP2_read_though[GSP2] += 1

            #check if the GSP2 is completely covered bu not at the start
            elif aligned_read.alen == GSP2_lengths[aligned_read.tid]:
                full_length = True
                if found_at_start:
                    GSP2_read_though[GSP2_name] += 1
                    #print "Read Through - " + GSP2_name
                #Save for later in case the next alignment is at the start
                else:
                    potential_read_through.append(GSP2_name)

            print "%s\t%s\t%s\t%s\t%s" % (aligned_read.qname, GSP2_name,
                                          aligned_read.cigarstring, found_at_start, full_length)

    for GSP2 in GSP2_references:
        print "%s\t%d\t%d" % (GSP2, GSP2_start[GSP2], GSP2_read_though[GSP2])