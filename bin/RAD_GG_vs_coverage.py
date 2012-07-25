#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 20 October 2011
@author: tcezard
'''
import sys, os
from utils import  utils_param, utils_logging
import logging, threading
from utils.parameters import Config_file_error
from optparse import OptionParser
from utils.utils_commands import get_output_stream_from_command

class FuncThread(threading.Thread):
    def __init__(self, target, *args):
        self._target = target
        self._args = args
        threading.Thread.__init__(self)
 
    def run(self):
        self._target(*self._args)
 

class AllContigsInfo():
    def __init__(self):
        self.per_contigs_info={}
        self.per_contig_per_sample_info={}
        self.lock = threading.RLock()
        
    def add_values(self,contig, nb_GG, coverage, duplicate, sample=None):
        with self.lock:
            info = self.per_contigs_info.get(contig)
            if not info:
                info={'nb_GG':0,'coverage':0, 'duplicate':0}
                self.per_contigs_info[contig]=info
            
            info['nb_GG']+=nb_GG
            info['coverage']+=coverage
            info['duplicate']+=duplicate
            if sample:
                info[sample]={'nb_GG':nb_GG,'coverage':coverage, 'duplicate':duplicate}
            
    def get_all_contigs(self):
        return self.per_contigs_info.keys()
    
    def get_contig_GG_value(self, contig, sample=None):
        info = self.per_contigs_info.get(contig)
        if info:
            if sample:
                tmp = info.get(sample)
                if tmp: 
                    return tmp.get('nb_GG')
            else:
                return info.get('nb_GG')
        
    def get_contig_coverage_value(self, contig,sample=None):
        info = self.per_contigs_info.get(contig)
        if info:
            if sample:
                tmp = info.get(sample)
                if tmp: 
                    return tmp.get('coverage')
            else:
                return info.get('coverage')
            
    def get_contig_duplicate_value(self, contig,sample=None):
        info = self.per_contigs_info.get(contig)
        if info:
            if sample:
                tmp = info.get(sample)
                if tmp: 
                    return tmp.get('duplicate')
            else:
                return info.get('duplicate')

def process_single_samtools_run(bam_file,all_contigs_info,samtools_bin):
    command="%s view -F 132 %s"%(samtools_bin, bam_file)
    open_stream, process=get_output_stream_from_command(command)
    current_contig=None
    nb_GG=0
    coverage=0
    duplicate=0
    sample_name, ext = os.path.splitext(bam_file)
    for line in open_stream:
        sp_line=line.strip().split()
        if current_contig!=sp_line[2] and current_contig != None:
            all_contigs_info.add_values(current_contig,nb_GG,coverage,duplicate,sample = sample_name)
            nb_GG=0
            coverage=0
            duplicate=0
        current_contig=sp_line[2]
        if int(sp_line[3])==1:
            if sp_line[9][4:6]=='GG':
                nb_GG+=1
            if int(sp_line[1]) & 1024 == 1024:
                duplicate+=1
            coverage+=1
            
    open_stream.close()


def RAD_median_coverage(bam_files,output_file):
    try:
        pipeline_param=utils_param.get_pipeline_parameters()
        samtools_dir=pipeline_param.get_samtools_dir()
    except Config_file_error, e:
        #logging.exception('Config_file_error:')
        logging.warning("You'll need to have samtools in your path")
        samtools_dir=''
    samtools_bin=os.path.join(samtools_dir,"samtools")
    all_threads=[]
    all_contigs_info=AllContigsInfo()
    sample_names=[]
    for bam_file in bam_files:
        sample_name, ext = os.path.splitext(bam_file)
        sample_names.append(sample_name)
        t = FuncThread(process_single_samtools_run, bam_file, all_contigs_info, samtools_bin)
        t.start()
        all_threads.append(t)
    for t in all_threads:
        t.join()
    all_contigs = all_contigs_info.get_all_contigs()
    all_contigs.sort()
    open_output=utils_logging.open_output_file(output_file)
    open_output.write("#contig\tnb_GG\tcoverage\tduplicate")
    for sample in sample_names:
        open_output.write("\t%s\t\t"%sample)
    open_output.write("\n")
    for contig in all_contigs:
        
        nb_GG=all_contigs_info.get_contig_GG_value(contig)
        coverage=all_contigs_info.get_contig_coverage_value(contig)
        duplicate=all_contigs_info.get_contig_duplicate_value(contig)
        open_output.write("%s\t%s\t%s\t%s"%(contig,nb_GG,coverage,duplicate))
        for sample in sample_names:
            nb_GG = all_contigs_info.get_contig_GG_value(contig,sample)
            coverage = all_contigs_info.get_contig_coverage_value(contig,sample)
            duplicate=all_contigs_info.get_contig_duplicate_value(contig,sample)
            if not nb_GG:
                nb_GG=0
            if not coverage:
                coverage=0
            if not coverage:
                duplicate=0
            open_output.write("\t%s\t%s\t%s"%(nb_GG,coverage,duplicate))
        open_output.write("\n")
        
    
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
    
