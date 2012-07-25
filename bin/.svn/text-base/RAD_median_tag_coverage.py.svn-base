#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 20 October 2011
@author: tcezard
'''
import sys, os
from utils.Binning import Distribution_holder
from utils import get_mpileup_from_bam, utils_param, utils_logging
import logging
from utils.parameters import Config_file_error
from optparse import OptionParser


def RAD_median_coverage(bam_files,output_file):
    try:
        pipeline_param=utils_param.get_pipeline_parameters()
        samtools_dir=pipeline_param.get_samtools_dir()
    except Config_file_error, e:
        #logging.exception('Config_file_error:')
        logging.warning("You'll need to have samtools in your path")
        samtools_dir=''
    samtools_bin=os.path.join(samtools_dir,"samtools")
    bam_file_str=' '.join(bam_files)
    all_dists=[]
    pileup_stream = get_mpileup_from_bam(bam_file_str, genome_file=None, samtools_bin=samtools_bin, options="-d 100000 -A")
    if output_file:
        open_output=utils_logging.open_output_file(output_file)
    else:
        open_output=sys.stdout
    bam_file_names=[]
    for file in bam_files:
        bam_file_names.append(os.path.basename(file))
    open_output.write("Consensus\t%s\n"%("\t".join(bam_file_names))) 
    line = pileup_stream.readline()
    sp_line=line.strip().split()
    curr_contig=sp_line[0]
    for i in range(len(sp_line)/3-1):
        all_dists.append(Distribution_holder());
        all_dists[i].add_value(sp_line[(i+1)*3])
    
    for line in pileup_stream:
        sp_line=line.strip().split()
        if curr_contig!=sp_line[0]:
            open_output.write(curr_contig)
            for i in range(len(sp_line)/3-1):
                open_output.write('\t%s'%all_dists[i].get_percentiles(50)[0])
                all_dists[i]=Distribution_holder();
                all_dists[i].add_value(sp_line[(i+1)*3])
            curr_contig=sp_line[0]
            open_output.write("\n") 
        for i in range(len(sp_line)/3-1):
            all_dists[i].add_value(sp_line[(i+1)*3])
    pileup_stream.close()
    open_output.close()
    
def main():
    #initialize the logging
    utils_logging.init_logging()
    #Setup options
    optparser=_prepare_optparser()
    (options,args) = optparser.parse_args()
    #verify options
    arg_pass=_verifyOption(options)
    if not arg_pass:
        logging.warning(optparser.get_usage())
        logging.critical("Non valid arguments: exit")
        sys.exit(1)
    utils_logging.change_log_stdout_to_log_stderr()
    bam_files=[options.bam_files]
    for file in args:
        if os.path.exists(file):
            bam_files.append(file)
    RAD_median_coverage(bam_files,options.output_file)
    

def _prepare_optparser():
    """Prepare optparser object. New options will be added in this
    function first.
    """
    usage = """usage: %prog <-b bam_file> [ -o output_file]"""
    description = """This script will take aligned RAD read to the consensuses and calculate per consensus coverage."""
    
    optparser = OptionParser(version="None",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="show this help message and exit.")
    optparser.add_option("-b","--bam_files",dest="bam_files",type="string",
                         help="Path to one or several bam files from which the coverage should be extracted. Default: %default")
    optparser.add_option("-o","--output_file",dest="output_file",type="string",
                         help="The path to the file where the data shoule be output. If not set, the results will be print to stdout. Default: %default")
    
    return optparser


def _verifyOption(options):
    """Check if the mandatory option are present in the options objects.
    @return False if any argument is wrong."""
    arg_pass=True
    
    if not options.bam_files :
        logging.error("You must specify a bam file.")
        arg_pass=False
    return arg_pass



if __name__=="__main__":
    main()
    
