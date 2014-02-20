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
from RAD_merge_bam_files import merge_bam_files
from utils.parameters import Config_file_error
from utils.utils_commands import get_output_stream_from_command
from collections import Counter
import multiprocessing


#get the path to the current script to infer the path to RAD_set_read1_consensus_to_read2.py
RADmapper_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
path_to_picard = os.path.join(RADmapper_dir,"picard")
mergeSamFilesWithCat_jar=os.path.join(path_to_picard,'MergeSamFilesWithCat.jar')

##### Generic merging functions
def merge_bam_files_with_picard(list_of_file, output_file=None, **kwargs):
    """This is a generic merging function for bam files.
    It assumes that all the bam file comes from mapping to independent contigs"""

    if not output_file:
        #Create a generic name and put it in the current working directory
        working_directory=os.getcwd()
        i=1
        output_file_template=os.path.join(working_directory,'tmp_merge_bam_%s.bam')
        output_file=output_file_template%i
        while os.path.exists(output_file):
            i+=1
            output_file=output_file_template%i
    command = 'java -jar -Xmx2G %s VALIDATION_STRINGENCY=SILENT CAT_SEQUENCE_DICTIONARIES=True USE_THREADING=True O=%s '%(mergeSamFilesWithCat_jar,output_file)
    inputs=['I=%s'%file for file in list_of_file]
    command += ' '.join(inputs)
    return_code=command_runner.run_command(command)
    if return_code==0:
        return output_file
    else:
        return None


def concatenate_file(list_of_file,output_file=None, **kwargs):
    """This is a generic merging function for concatenating text files.
    It can take a filter keyword argument to grep out using the provided value"""
    if not output_file:
        #Create a generic name and put it in the current working directory
        working_directory=os.getcwd()
        i=1
        output_file_template=os.path.join(working_directory,'tmp_concatenate_%s')
        output_file=output_file_template%i
        while os.path.exists(output_file):
            i+=1
            output_file=output_file_template%i

    if kwargs.has_key('filter'):
        filter_on = kwargs.get('filter')
        command = 'cat %s | egrep -v %s > %s '%(' '.join(list_of_file), filter_on, output_file)
    else:
        command = 'cat %s > %s '%(' '.join(list_of_file), output_file)
    return_code=command_runner.run_command(command)
    if return_code==0:
        return output_file
    else:
        return None


def merge_by_chunck(file_to_merge, function_to_merge, output_file=None, max_nb_file=100, **kwargs):
    """This function merge file using a generic merge function. It merges chunk of max_nb_file at one time"""
    if len(file_to_merge) == 1:
        if output_file:
            os.rename(file_to_merge[0], output_file)
        else:
            output_file=file_to_merge[0]
    elif len(file_to_merge)>max_nb_file:
        new_file_to_merge=[]
        for i in range(0,len(file_to_merge),max_nb_file):
            tmp_merge_file = function_to_merge(file_to_merge[i:i+max_nb_file], **kwargs)
            new_file_to_merge.append(tmp_merge_file)
        output_file = merge_by_chunck(new_file_to_merge, function_to_merge, output_file, **kwargs)
        for tmp_file in new_file_to_merge:
            logging.info('Remove %s'%tmp_file)
            os.remove(tmp_file)
    else:
        output_file= function_to_merge(file_to_merge, output_file, **kwargs)
    return output_file


    
def merge_all_bam_files_from_directory(directory):
    """This function merge bam file from a single directory"""
    directory=os.path.abspath(directory)
    all_bam_files = glob(os.path.join(directory,'*_dir','*_corrected_sorted_mrk_dup_fixed.bam'))
    output_file = os.path.join(directory,'%s_files.bam'%len(all_bam_files))
    output_file = merge_bam_files_with_picard(all_bam_files, output_file)
    if not output_file:
        logging.error("Merging bam files in %s failed"%(directory))
        #TODO do something about it
    return output_file


def merge_contigs_files(directory):
    """This function merge bam file from a single directory"""
    all_fasta_files = glob(os.path.join(directory,'*_dir','best_assembly.fa'))
    output_file = os.path.join(directory,'%s_best_assembly.fa'%len(all_fasta_files))
    output_file = concatenate_file(all_fasta_files,output_file)
    if not output_file:
        logging.error("Merging assemblies in %s failed"%(directory))
        #TODO do something about it
    return output_file


def merge_snps_files(directory):
    """This function merge snps files from a single directory"""
    return_code=0
    all_vcf_files = glob(os.path.join(directory,'*_dir','*samtools.vcf'))
    output_file_body = os.path.join(directory,'%s_snps_files.vcf.body'%len(all_vcf_files))
    output_file_body = concatenate_file(all_vcf_files, output_file_body, filter="^#")
    if output_file_body:
        return_code=0
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
    """This function will merge the bam files across all the directories"""
    directory=os.path.abspath(directory)
    all_bam_files = glob(os.path.join(directory,'*_dir','*_files.bam'))
    #Need to sort as glob retuns the file in random order
    all_bam_files.sort()
    output_file = os.path.join(directory,'all_consensus_merged.bam')
    merge_bam_files(all_bam_files,output_file)

    #output_file = merge_by_chunck(all_bam_files, merge_bam_files_with_picard, output_file)
    if output_file:
        return 0

def merge_all_contigs_files_from_directories(directory):
    """This function will merge the contigs files across all the directories"""
    all_fasta_files = glob(os.path.join(directory,'*_dir','*_best_assembly.fa'))
    #Need to sort as glob retuns the file in random order
    all_fasta_files.sort()

    output_file = os.path.join(directory,'all_consensus_assembly.fa')
    output_file = merge_by_chunck(all_fasta_files, concatenate_file, output_file)
    if output_file:
        return 0


def merge_all_snps_files_from_directories(directory):
    """This function will merge the snps files across all the directories"""
    return_code=0
    all_vcf_files = glob(os.path.join(directory,'*_dir','*_snps_files.vcf'))
    output_file_body = os.path.join(directory,'all_consensus_snps_files.vcf.body')
    output_file_body = merge_by_chunck(all_vcf_files, concatenate_file, output_file_body, filter="^#")
    if output_file_body:
        return_code=0
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

def merge_all_summary_files_from_directories(directory):
    """This function will merge the summary files across all the directories"""
    return_code=0
    all_summary_files = glob(os.path.join(directory,'*_dir','*summary_stat.txt'))
    output_file_body = os.path.join(directory,'all_summary_stat.txt.body')
    output_file_body = merge_by_chunck(all_summary_files, concatenate_file, output_file_body, filter="^name")
    if output_file_body:
        return_code=0
    output_file_header = os.path.join(directory,'all_summary_stat.txt.header')
    command = 'head -n 1 %s > %s '%(all_summary_files[0], output_file_header)
    if return_code==0:
        return_code = command_runner.run_command(command)
    output_file = os.path.join(directory,'all_summary_stat.txt')
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
    return_code = merge_all_bam_files_from_directory(directory)
    return_code = merge_snps_files(directory)

    return return_code

def merge_all_results(directory):
    return_code=0
    if return_code==0:
        return_code = merge_all_contigs_files_from_directories(directory)
    if return_code==0:
        return_code = merge_all_snps_files_from_directories(directory)
    if return_code==0:
        return_code = merge_all_summary_files_from_directories(directory)
    if return_code==0:
        return_code = merge_all_bam_files_from_directories(directory)
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
        code = merge_all_results(options.consensus_dir)
    else:
        code = merge_results(options.consensus_dir)
    sys.exit(code)

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
    
    
