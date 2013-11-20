#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 19 March  2012
@author: tcezard
'''
import sys, os
from utils import utils_logging, utils_commands, longest_common_substr_from_start, utils_param, sort_bam_file_per_coordinate
import logging, threading, re
from optparse import OptionParser
from glob import glob
import command_runner
from utils.FastaFormat import FastaReader
import time
from RAD_merge_read1_and_read2 import merge_2_contigs
from utils.parameters import Config_file_error
from RAD_assemble_read2 import run_all_fastq_files
from utils.utils_commands import get_output_stream_from_command
from collections import Counter
import multiprocessing

max_nb_file=100

def submit_command_to_sge(command):
    print "submit %s to cluster"%command
    
def merge_all_bam_files(all_bam_files, current_directory):
    if len(all_bam_files) == 1:
        return all_bam_files[0]
    elif len(all_bam_files)>max_nb_file:
        result_from_merge=[]
        for i in range(0,len(all_bam_files),max_nb_file):
            output_file=''
            result_from_merge.append(output_file)
            #commands.append(generate_merge_command(all_bam_files[i:i+max_nb_file], output_file))
            
        return merge_all_bam_files(result_from_merge)
    else:
        #Merge max_nb_file or less bam files
        output_file=''
        generate_merge_command(all_bam_files, output_file)
        
def generate_merge_command(all_bam_files, output_file):
    command = 'java -jar -Xmx2G  ~/workspace/Picard/dist/MergeSamFilesWithCat.jar VALIDATION_STRINGENCY=SILENT CAT_SEQUENCE_DICTIONARIES=True USE_THREADING=True O=%s '%(output_file)
    inputs=['I=%s'%file for file in all_bam_files]
    command += ' '.join(inputs)
    
        

def merge_all_bam_files_from_directory(directory):
    directory=os.path.abspath(directory)
    all_bam_files = glob(os.path.join(directory,'*_dir','*_dir','*_corrected_sorted_mrk_dup_fixed.bam'))
    bam_file = merge_all_bam_files(all_bam_files)
    
def merge_bam_files(directory):
    directory=os.path.abspath(directory)
    all_bam_files = glob(os.path.join(directory,'*_dir','*_corrected_sorted_mrk_dup_fixed.bam'))
    output_file = os.path.join(directory,'%s_files.bam'%len(all_bam_files))
    command = 'java -jar -Xmx2G  ~/workspace/Picard/dist/MergeSamFilesWithCat.jar VALIDATION_STRINGENCY=SILENT CAT_SEQUENCE_DICTIONARIES=True USE_THREADING=True O=%s '%(output_file)
    inputs=['I=%s'%file for file in all_bam_files]
    command += ' '.join(inputs)
    return command_runner.run_command(command)


def merge_contigs_files(directory):
    all_fasta_files = glob(os.path.join(directory,'*_dir','best_assembly.fa'))
    output_file = os.path.join(directory,'%s_best_assembly.fa'%len(all_fasta_files))
    command = 'cat %s > %s '%(' '.join(all_fasta_files), output_file)
    return command_runner.run_command(command)

def merge_snps_files(directory):
    return_code=0
    all_vcf_files = glob(os.path.join(directory,'*_dir','*samtools.vcf'))
    output_file_body = os.path.join(directory,'%s_snps_files.vcf.body'%len(all_vcf_files))
    command = 'cat %s  | egrep -v "^#" > %s '%(' '.join(all_vcf_files), output_file_body)
    if return_code==0:
        return_code = command_runner.run_command(command)
    output_file_header = os.path.join(directory,'%s_snps_files.vcf.header'%len(all_vcf_files))
    command = 'grep  "^#" %s > %s '%(all_vcf_files[0], output_file_header)
    if return_code==0:
        return_code = command_runner.run_command(command)
    output_file = os.path.join(directory,'%s_snps_files.vcf'%len(all_vcf_files))
    command = 'cat %s %s > %s '%(output_file_header, output_file_body, output_file)
    if return_code==0:
        return_code = command_runner.run_command(command)
    command = 'rm %s %s'%(output_file_header, output_file_body)
    if return_code==0:
        return_code = command_runner.run_command(command)
    return return_code


def merge_all_bam_files_from_directories(directory):
    directory=os.path.abspath(directory)
    all_bam_files = glob(os.path.join(directory,'*_dir','*_files.bam'))
    output_file = os.path.join(directory,'all_consensus_merged.bam')
    command = 'java -jar -Xmx2G  ~/workspace/Picard/dist/MergeSamFilesWithCat.jar VALIDATION_STRINGENCY=SILENT CAT_SEQUENCE_DICTIONARIES=True USE_THREADING=True O=%s '%(output_file)
    inputs=['I=%s'%file for file in all_bam_files]
    command += ' '.join(inputs)
    return command_runner.run_command(command)

def merge_all_contigs_files_from_directories(directory):
    all_fasta_files = glob(os.path.join(directory,'*_dir','*_best_assembly.fa'))
    output_file = os.path.join(directory,'all_consensus_assembly.fa')
    command = 'cat %s > %s '%(' '.join(all_fasta_files), output_file)
    return command_runner.run_command(command)

def merge_all_snps_files_from_directories(directory):
    return_code=0
    all_vcf_files = glob(os.path.join(directory,'*_dir','*_snps_files.vcf'))
    output_file_body = os.path.join(directory,'all_consensus_snps_files.vcf.body')
    command = 'cat %s  | egrep -v "^#" > %s '%(' '.join(all_vcf_files), output_file_body)
    if return_code==0:
        return_code = command_runner.run_command(command)
    output_file_header = os.path.join(directory,'all_consensus_snps_files.vcf.header')
    command = 'grep  "^#" %s > %s '%(all_vcf_files[0], output_file_header)
    if return_code==0:
        return_code = command_runner.run_command(command)
    output_file = os.path.join(directory,'all_consensus_snps_files.vcf')
    command = 'cat %s %s > %s '%(output_file_header, output_file_body, output_file)
    if return_code==0:
        return_code = command_runner.run_command(command)
    command = 'rm %s %s'%(output_file_header, output_file_body)
    if return_code==0:
        return_code = command_runner.run_command(command)
    return return_code



def merge_results(directory):
    return_code=0
    return_code = merge_contigs_files(directory)
    #if return_code==0:
    #    return_code = merge_bam_files(directory)
    if return_code==0:
        return_code = merge_snps_files(directory)
    return return_code

def merge_all_results(directory):
    return_code=0
    #return_code = merge_all_bam_files_from_directories(directory)
    if return_code==0:
        return_code = merge_all_contigs_files_from_directories(directory)
    if return_code==0:
        return_code = merge_all_snps_files_from_directories(directory)
    return return_code

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
    if not options.print_commands:
        command_runner.set_command_to_run_localy()
    if options.debug:
        utils_logging.init_logging(logging.DEBUG)
    if options.final_merge:
        return merge_all_results(options.consensus_dir)
    else:
        return merge_results(options.consensus_dir)
    

def _prepare_optparser():
    """Prepare optparser object. New options will be added in this
    function first.
    """
    usage = """usage: %prog <-b bam_file> [ -o output_file]"""
    description = """This script will take aligned RAD read to the consensuses and calculate per consensus coverage."""
    
    optparser = OptionParser(version="None",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="show this help message and exit.")
    optparser.add_option("-d","--consensus_dir",dest="consensus_dir",type="string",
                         help="Path to a directory containing fastq file (only extension .fastq will be processed). Default: %default")
    optparser.add_option("--final_merge",dest="final_merge",action='store_true',default=False,
                         help="Merge the already merged file. Default: %default")
    optparser.add_option("--print",dest="print_commands",action='store_true',default=False,
                         help="print commands instead of running them. Default: %default")
    optparser.add_option("--debug",dest="debug",action='store_true',default=False,
                         help="Output debug statements. Default: %default")
    return optparser


def _verifyOption(options):
    """Check if the mandatory option are present in the options objects.
    @return False if any argument is wrong."""
    arg_pass=True
    
    return arg_pass


if __name__=="__main__":
    main()
    
    