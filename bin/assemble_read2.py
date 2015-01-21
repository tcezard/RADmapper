#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 14 November 2013
@author: tcezard
'''
import sys
import os
import logging
from optparse import OptionParser

import command_runner
import RAD_assemble_read2
import RAD_smalt_align_reads_back_to_consensus
import RAD_merge_results
import RAD_bam_to_fastq
import RAD_summarize_assembly
from utils import utils_logging
from utils import utils_commands


def get_readgroup_from_bam(bam_files):
    all_read_groups=[]
    for bam_file in bam_files:
        command = "samtools view -H %s | grep '^@RG' " % bam_file
        stdout, process = utils_commands.get_output_stream_from_command(command)
        for line in stdout:
            all_read_groups.append(line.strip())
    return all_read_groups


def extract_and_assemble(bam_files, genome_file, white_list_file, output_dir, assembly_function_list,force_merge=False):
    logging.info('processing %s'%output_dir)
    logging.info("--------------------------")
    logging.info("Extract the reads from bam")
    logging.info("--------------------------")
    list_consensus=RAD_bam_to_fastq.read_white_list(white_list_file)
    RAD_bam_to_fastq.extract_reads_from_all_bam_files_set_of_consensus(bam_files, list_consensus=list_consensus,
                                                                       output_dir=output_dir, all_read1_consensus_file=genome_file )
    logging.info("--------------")
    logging.info("Assemble read2")
    logging.info("--------------")
    RAD_assemble_read2.run_all_fastq_files(output_dir, assembly_function_list, 600, force_merge=force_merge)
    logging.info("----------------------------")
    logging.info("Extract read groups from bam")
    logging.info("----------------------------")
    all_read_groups=get_readgroup_from_bam(bam_files)
    logging.info("-------------------------------------")
    logging.info("Align both reads back to the assembly")
    logging.info("-------------------------------------")
    RAD_smalt_align_reads_back_to_consensus.run_all_fastq_files(output_dir, all_read_groups=all_read_groups, snp_call=True)
    logging.info("---------------------")
    logging.info("Merge all the results")
    logging.info("---------------------")
    RAD_merge_results.merge_results(output_dir)
    logging.info("------------------------------")
    logging.info("Get the summary of the results")
    logging.info("------------------------------")
    RAD_summarize_assembly.run_all_summaries(output_dir)


def main():
    #initialize the logging
    utils_logging.init_logging(logging.INFO)
    #Setup options
    optparser=_prepare_optparser()
    (options,args) = optparser.parse_args()
    #verify options
    arg_pass=_verifyOption(options)
    if not arg_pass:
        logging.warning(optparser.get_usage())
        logging.critical("Non valid arguments: exit")
        sys.exit(1)
    if options.debug:
        utils_logging.init_logging(logging.DEBUG)
    utils_logging.init_logging(output_level=None,
                               log_file_name=os.path.join(options.output_dir,"extract_and_assemble.log"))

    bam_files=[options.bam_files]
    if len(args)>0:
        bam_files.extend(args)

    all_assembler_to_try = options.assembler_name.split(',')
    assembly_function_list=[]
    for assembler_name in all_assembler_to_try:
        assembly_function=RAD_assemble_read2.get_assembly_function(assembler_name)
        assembly_function_list.append(assembly_function)
    command_runner.set_command_to_run_localy()
    extract_and_assemble(bam_files, options.read1_consensus_file, options.white_list_file, options.output_dir, assembly_function_list)

def _prepare_optparser():
    """Prepare optparser object. New options will be added in this
    function first.
    """
    usage = """usage: %prog <-b bam_file> [ -o output_file]"""
    description = """This script will take aligned RAD read to the consensuses and calculate per consensus coverage."""
    
    optparser = OptionParser(version="None",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="show this help message and exit.")
    optparser.add_option("-b","--bam_files",dest="bam_files",type="string",default=None,
                         help="The path to bam files that will be used to extrat the reads. Default: %default")
    optparser.add_option("-a","--assembler_name",dest="assembler_name",type="string",default=None,
                         help="The name of the assembler that will be used on the fastq files. Default: %default")
    optparser.add_option("-g","--read1_consensus_file",dest="read1_consensus_file",type="string",
                         help="The fasta file containing the read1 consensus file corresponding to the bam file.")
    optparser.add_option("-w","--white_list_file",dest="white_list_file",type="string",
                         help="The file containing the name of the consensus to assemble")
    optparser.add_option("-o","--output_dir",dest="output_dir",type="string",
                         help="This act as a flag to bypass the nb_consensus_per_dir. This means that all consensus in the whitelist will be extracted in the output_dir.")
    optparser.add_option("--force_merge",dest="force_merge",action='store_true',default=False,
                         help="Force merged the consensus: add a run of 100 N in between each sequence. Default: %default")
    optparser.add_option("--print",dest="print_command",action='store_true',default=False,
                         help="print the commands instead of running them. Default: %default")
    optparser.add_option("--debug",dest="debug",action='store_true',default=False,
                         help="Output debug statment. Default: %default")
    return optparser


def _verifyOption(options):
    """Check if the mandatory option are present in the options objects.
    @return False if any argument is wrong."""
    arg_pass=True
    
    return arg_pass



if __name__=="__main__":
    main()
    
