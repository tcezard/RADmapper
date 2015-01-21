#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 29 August 2013
@author: tcezard
'''
import sys
from optparse import OptionParser


def generate_ustack_command(fastq_file, output_dir, min_depth, max_dist_1, max_dist_2, counter):
    """ustacks -t file_type -f file_path [-d] [-r] [-o path] [-i id] [-e errfreq] [-m min_cov] [-M max_dist] [-p num_threads] [-R] [-H] [-h]
  p: enable parallel execution with num_threads threads.
  t: input file Type. Supported types: fasta, fastq, bowtie, sam, tsv.
  f: input file path.
  o: output path to write results.
  i: SQL ID to insert into the output to identify this sample.
  m: Minimum depth of coverage required to create a stack (default 2).
  M: Maximum distance (in nucleotides) allowed between stacks (default 2).
  N: Maximum distance allowed to align secondary reads to primary stacks (default: M + 2).
  d: enable the Deleveraging algorithm, used for resolving over merged tags.
  r: enable the Removal algorithm, to drop highly-repetitive stacks (and nearby errors) from the algorithm.
  e: specify the barcode error frequency (0 < e < 1) if using the 'fixed' model.
  R: retain unused reads.
  H: disable calling haplotypes from secondary reads.
  h: display this help messsage.
"""
    ustacks_bin = 'ustacks'
    command = "%s -t fastq -p 2 -m %s -M %s -N %s -H -o %s -f %s -i %s" % (
    ustacks_bin, min_depth, max_dist_1, max_dist_2, output_dir, fastq_file, counter)
    return command


def run_on_cluster(commands):
    pass


def cluster_read1_one_sample(fastq_file):
    command = generate_ustack_command(fastq_file)


def cluster_read1_accross_all_samples(fastq_files):
    for fastq_file in fastq_files:
        cluster_read1_one_sample(fastq_file)


def main():
    # initialize the logging
    utils_logging.init_logging(logging.INFO)
    #utils_logging.init_logging(logging.CRITICAL)
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
    if not options.print_command:
        command_runner.set_command_to_run_localy()
    run_all_fastq_files(options.consensus_dir)


def _prepare_optparser():
    """Prepare optparser object. New options will be added in this
    function first.
    """
    usage = """usage: %prog <-b bam_file> [ -o output_file]"""
    description = """This script will take aligned RAD read to the consensuses and calculate per consensus coverage."""

    optparser = OptionParser(version="None", description=description, usage=usage, add_help_option=False)
    optparser.add_option("-h", "--help", action="help", help="show this help message and exit.")
    optparser.add_option("-d", "--consensus_dir", dest="consensus_dir", type="string",
                         help="Path to a directory containing fastq file (only extension .fastq will be processed). Default: %default")
    optparser.add_option("--print", dest="print_command", action='store_true', default=False,
                         help="print the commands instead of running them. Default: %default")
    optparser.add_option("--debug", dest="debug", action='store_true', default=False,
                         help="Output debug statment. Default: %default")
    return optparser


def _verifyOption(options):
    """Check if the mandatory option are present in the options objects.
    @return False if any argument is wrong."""
    arg_pass = True

    return arg_pass


if __name__ == "__main__":
    main()
