#!/usr/bin/env python
import sys
import os
import logging
from optparse import OptionParser
import glob
import shutil

from utils.GenomeLoader import GenomeLoader
from utils import utils_logging
import pysam
import command_runner


def read_white_list(white_list_file):
    all_consensus = []
    with open(white_list_file) as open_file:
        for line in open_file:
            sp_line = line.strip().split()
            if sp_line and len(sp_line) > 0:
                all_consensus.append(sp_line[0])
    return all_consensus


_all_open_bam = {}


def _get_open_bam_file(bam_file):
    if _all_open_bam.has_key(bam_file):
        return _all_open_bam.get(bam_file)

    logging.info('open %s' % bam_file)
    open_bam = pysam.Samfile(bam_file, "rb")
    _all_open_bam[bam_file] = open_bam
    return open_bam


def close_bam_file(bam_file):
    if _all_open_bam.has_key(bam_file):
        open_bam = _all_open_bam.pop(bam_file)
        logging.info('close %s' % bam_file)
        open_bam.close()


def get_pysam_iterator(bam_file, options, consensus):
    sam_file = _get_open_bam_file(bam_file)
    return sam_file.fetch(consensus)


def get_read_group_from_bam_files(bam_files):
    all_read_groups = []
    for bam_file in bam_files:
        open_bam_file = _get_open_bam_file(bam_file)
        read_groups = open_bam_file.header.get('RG')
        out = ['@RG']
        for read_group in read_groups:
            out.append('ID:%s' % read_group.pop('ID'))
            for key in read_group.keys():
                out.append('%s:%s' % (key, read_group.get(key)))
        all_read_groups.append('\t'.join(out))
    return all_read_groups


def load_from_sites_generator2(bam_file, options='', consensus=''):
    """This function return a generator that iterates over read pair where both pairs are mapping.
    @return a tuple containing read1 read2 and the restriction site position."""
    iter_sam = get_pysam_iterator(bam_file, options=options, consensus=consensus)
    all_unmatched_read1 = {}
    all_unmatched_read2 = {}
    count_line = 0
    for align_read in iter_sam:
        count_line += 1
        if align_read.is_read1:
            align_read_r2 = all_unmatched_read2.pop(align_read.qname, None)
            if align_read_r2:
                yield ((align_read, align_read_r2))
            else:
                all_unmatched_read1[align_read.qname] = align_read
        else:
            sam_record_r1 = all_unmatched_read1.pop(align_read.qname, None)
            if sam_record_r1:
                yield ((sam_record_r1, align_read))
            else:
                all_unmatched_read2[align_read.qname] = align_read


def process_one_bam_file_one_consensus(bam_file, consensus_name):
    sam_pair_generator = load_from_sites_generator2(bam_file, consensus=consensus_name)
    all_first_reads_for_consensus = []
    all_second_reads_for_consensus = []
    for aligned_read_r1, aligned_read_r2, in sam_pair_generator:
        rgid = aligned_read_r1.opt('RG')
        all_first_reads_for_consensus.append(aligned_read_to_fastq(aligned_read_r1, rgid=rgid))
        all_second_reads_for_consensus.append(aligned_read_to_fastq(aligned_read_r2, rgid=rgid))
        aligned_read_r1 = None
        aligned_read_r2 = None
    return all_first_reads_for_consensus, all_second_reads_for_consensus


def extract_reads_from_one_consensus(bam_files, output_dir, consensus_name, consensus_sequence):
    all_read1_for_that_consensus = []
    all_read2_for_that_consensus = []
    for bam_file in bam_files:
        all_first_reads_for_consensus, all_second_reads_for_consensus = process_one_bam_file_one_consensus(bam_file,
                                                                                                           consensus_name)
        all_read1_for_that_consensus.extend(all_first_reads_for_consensus)
        all_read2_for_that_consensus.extend(all_second_reads_for_consensus)

    consensus_directory = os.path.join(output_dir, consensus_name + '_dir')
    if not os.path.exists(consensus_directory):
        os.mkdir(consensus_directory)
    read1_file = os.path.join(consensus_directory, consensus_name + "_1.fastq")
    read2_file = os.path.join(consensus_directory, consensus_name + "_2.fastq")
    read1_consensus = os.path.join(consensus_directory, consensus_name + "_1.fa")
    open_file = open(read1_consensus, 'w')
    open_file.write('>%s\n%s\n' % (consensus_name, consensus_sequence))
    open_file.close()

    open_file = open(read1_file, 'w')
    open_file.write('\n'.join(all_read1_for_that_consensus))
    open_file.close()

    open_file = open(read2_file, 'w')
    open_file.write('\n'.join(all_read2_for_that_consensus))
    open_file.close()

    return read1_consensus, read1_file, read2_file


open_fastq_files = {}


def _get_open_fastq_files(fastq_file):
    if open_fastq_files.has_key(fastq_file):
        return open_fastq_files.get(fastq_file)
    open_file = open(fastq_file, 'w')
    open_fastq_files[fastq_file] = open_file
    return open_file


def close_fastq_files():
    for filename in open_fastq_files.keys():
        open_fastq_files.pop(filename).close()


def extract_reads_from_one_bam_file(bam_file, output_dir, list_consensus, genome_loader):
    for consensus_name in list_consensus:
        consensus_name, consensus_sequence = genome_loader.get_chr(consensus_name)
        all_first_reads_for_consensus, all_second_reads_for_consensus = process_one_bam_file_one_consensus(bam_file,
                                                                                                           consensus_name)
        consensus_directory = os.path.join(output_dir, consensus_name + '_dir')
        read1_file = os.path.join(consensus_directory, consensus_name + "_1.fastq")
        read2_file = os.path.join(consensus_directory, consensus_name + "_2.fastq")
        if not os.path.exists(consensus_directory):
            os.mkdir(consensus_directory)
            read1_consensus = os.path.join(consensus_directory, consensus_name + "_1.fa")
            open_file = open(read1_consensus, 'w')
            open_file.write('>%s\n%s\n' % (consensus_name, consensus_sequence))
            open_file.close()
            # open_file1 = open(read1_file,'w')
            #      open_file2 = open(read2_file,'w')
            #  else:
            #      open_file1 = open(read1_file,'a')
            #      open_file2 = open(read2_file,'a')

        if all_first_reads_for_consensus:
            open_file1 = _get_open_fastq_files(read1_file)
            open_file1.write('\n'.join(all_first_reads_for_consensus) + '\n')
        # open_file1.close()
        if all_second_reads_for_consensus:
            open_file2 = _get_open_fastq_files(read2_file)
            open_file2.write('\n'.join(all_second_reads_for_consensus) + '\n')
        logging.info("Extract %s reads from %s" % (len(all_first_reads_for_consensus), consensus_name))
        #open_file2.close()
    close_bam_file(bam_file)


def split_whitelist(white_list_file, nb_consensus_per_dir):
    list_consensus = read_white_list(white_list_file)
    list_of_list_consensus = []
    temp_list = []
    i = 0
    for consensus in list_consensus:
        i += 1
        temp_list.append(consensus)
        if i % nb_consensus_per_dir == 0:
            list_of_list_consensus.append(temp_list)
            temp_list = []
    if len(temp_list) > 0:
        list_of_list_consensus.append(temp_list)
    # print "split into %s jobs of %s consensus"%(len(list_of_list_consensus),nb_consensus_per_dir)
    for i, list_of_consensus in enumerate(list_of_list_consensus):
        output_dir = os.path.join(os.path.curdir, '%s_dir' % (i + 1))
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
        output_file = os.path.join(output_dir, 'whitelist.txt')
        with open(output_file, 'w') as open_file:
            open_file.write('\n'.join(list_of_consensus))
        yield (output_dir, output_file)


def extract_reads_from_all_bam_files_set_of_consensus_old(bam_files, list_consensus, output_dir, genome_loader=None,
                                                          all_read1_consensus_file=None):
    if genome_loader is None:
        genome_loader = GenomeLoader(all_read1_consensus_file, keep_until_done=True)
    for consensus_name in list_consensus:
        logging.info("Extract reads from %s " % consensus_name)
        consensus_name, consensus_sequence = genome_loader.get_chr(consensus_name)
        extract_reads_from_one_consensus(bam_files, output_dir, consensus_name, consensus_sequence)


def extract_reads_from_all_bam_files_set_of_consensus(bam_files, list_consensus, output_dir, genome_loader=None,
                                                      all_read1_consensus_file=None):
    all_previous_dir = glob.glob(os.path.join(output_dir, '*_dir'))
    if len(all_previous_dir):
        logging.info("cleanup previous run in %s" % output_dir)
        for dir in all_previous_dir:
            shutil.rmtree(dir)
    if genome_loader is None:
        genome_loader = GenomeLoader(all_read1_consensus_file, keep_until_done=True)
    for bam_file in bam_files:
        extract_reads_from_one_bam_file(bam_file, output_dir, list_consensus, genome_loader)
    # All the read have been extract now close the fastq files
    close_fastq_files()


def aligned_read_to_fastq(aligned_read, rgid=None):
    out = []
    if aligned_read.is_read1:
        suffix = '/1'
    elif aligned_read.is_read2:
        suffix = '/2'
    else:
        suffix = ''
    if rgid:
        rgid_str = 'RGID:%s' % (rgid)
    else:
        rgid_str = ''
    out.append('@%s%s%s' % (aligned_read.qname, rgid_str, suffix))
    if aligned_read.is_read2:
        # Do not reverse complement the fastq file
        #out.append(DNA_tools.rev_complements(aligned_read.seq))
        #out.append('+')
        #out.append(aligned_read.qual[::-1])
        out.append(aligned_read.seq)
        out.append('+')
        out.append(aligned_read.qual)
    else:
        out.append(aligned_read.seq)
        out.append('+')
        out.append(aligned_read.qual)
    return '\n'.join(out)


def main():
    # initialize the logging
    utils_logging.init_logging(logging.INFO)
    #Setup options
    optparser = _prepare_optparser()
    (options, args) = optparser.parse_args()
    #verify options
    arg_pass = _verifyOption(options)
    if not arg_pass:
        logging.warning(optparser.get_usage())
        logging.critical("Non valid arguments: exit")
        sys.exit(1)
    if options.debug:
        utils_logging.init_logging(logging.DEBUG)
    bam_files = [options.bam_file]
    if len(args) > 0:
        bam_files.extend(args)
    command_runner.set_command_to_run_localy()
    if options.output_dir:
        list_consensus = read_white_list(options.white_list_file)
        extract_reads_from_all_bam_files_set_of_consensus(bam_files, list_consensus, options.output_dir,
                                                          genome_loader=None,
                                                          all_read1_consensus_file=options.read1_consensus_file)
    else:
        iter_dir_and_file = split_whitelist(white_list_file=options.white_list_file,
                                            nb_consensus_per_dir=options.nb_consensus_per_dir)
        for output_dir, sub_whitelist in iter_dir_and_file:
            command = 'python %s -g %s -w %s -o %s -b %s' % (
            sys.argv[0], options.read1_consensus_file, sub_whitelist, output_dir, ' '.join(bam_files))
            print command


def _prepare_optparser():
    """Prepare optparser object. New options will be added in this
    function first.
    """
    usage = """usage: %prog <-b bam_file> <-g genome_file> <-k known_sites>"""
    description = """This script extract reads from an aligned bam file and create the corresponding fastq files."""

    optparser = OptionParser(description=description, usage=usage, add_help_option=False)
    optparser.add_option("-h", "--help", action="help", help="show this help message and exit.")
    optparser.add_option("-b", "--bam_file", dest="bam_file", type="string",
                         help="The bam file from which the reads should be extracted.")
    optparser.add_option("-g", "--read1_consensus_file", dest="read1_consensus_file", type="string",
                         help="The fasta file containing the read1 consensus file corresponding to the bam file.")
    optparser.add_option("-w", "--white_list_file", dest="white_list_file", type="string",
                         help="The file containing the name of the consensus to assemble")
    optparser.add_option("-o", "--output_dir", dest="output_dir", type="string",
                         help="This act as a flag to bypass the nb_consensus_per_dir. This means that all consensus in the whitelist will be extracted in the output_dir.")
    optparser.add_option("-p", "--number_of_process", dest="number_of_process", type="int", default=10,
                         help="The number of process used to extract the information from the bam files.")
    optparser.add_option("-n", "--nb_consensus_per_dir", dest="nb_consensus_per_dir", type="int", default=50,
                         help="The number of that should be held in each directory.")
    optparser.add_option("--debug", dest="debug", action='store_true', default=False,
                         help="Output debug statment. Default: %default")

    return optparser


def _verifyOption(options):
    """Check if the mandatory option are present in the options objects.
    @return False if any argument is wrong."""
    arg_pass = True

    if not options.bam_file:
        logging.error("You must specify a bam file -b.")
        arg_pass = False
    if not options.read1_consensus_file:
        logging.error("You must specify a fasta file with the read1 consensus with -g.")
        arg_pass = False
    elif not os.path.exists(options.read1_consensus_file):
        logging.error("You must specify a an existing fasta file with -g.")
        arg_pass = False
    if not options.white_list_file:
        logging.error("You must specify a list of consensus to use with -w.")
        arg_pass = False
    elif not os.path.exists(options.white_list_file):
        logging.error("You must specify a an existing list of consensus to use with -w.")
        arg_pass = False
    return arg_pass


if __name__ == "__main__":
    main()

if __name__ == "1__main__":
    bam_files = sys.argv[1:]
    for bam_file in bam_files:
        open_bam = _get_open_bam_file(bam_file)
    