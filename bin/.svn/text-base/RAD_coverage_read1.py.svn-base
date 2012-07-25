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
from IO_interface.samIterator import Sam_record

class FuncThread(threading.Thread):
    def __init__(self, target, *args):
        self._target = target
        self._args = args
        threading.Thread.__init__(self)
        self.returned_value=None
 
    def run(self):
        self.returned_value = self._target(*self._args)
        
    def get_returned_value(self):
        return self.returned_value
 

class AllContigsInfo():
    def __init__(self):
        self.per_contigs_info={}
        self.per_contig_per_sample_info={}
        self.sample_list=set()
        self.lock = threading.RLock()
        
    def add_values(self,contig, coverage, duplicate, sample=None):
        with self.lock:
            info = self.per_contigs_info.get(contig)
            if not info:
                info={'coverage':0, 'coverage_mrk_dup':0}
                self.per_contigs_info[contig]=info
            info['coverage']+=coverage
            info['coverage_mrk_dup']+=coverage-duplicate
            if sample:
                self.sample_list.add(sample)
                if info.has_key(sample):
                    info[sample]['coverage']+=coverage
                    info[sample]['coverage_mrk_dup']+=coverage-duplicate
                else:
                    info[sample]={'coverage':coverage, 'coverage_mrk_dup':coverage-duplicate}
            
    def get_all_contigs(self):
        return self.per_contigs_info.keys()
    
    def get_all_samples(self):
        return list(self.sample_list)
    
        
    def get_contig_coverage_value(self, contig,sample=None):
        info = self.per_contigs_info.get(contig)
        if info:
            if sample:
                tmp = info.get(sample)
                if tmp: 
                    return tmp.get('coverage')
            else:
                return info.get('coverage')
            
    def get_contig_coverage_mrk_dup_value(self, contig,sample=None):
        info = self.per_contigs_info.get(contig)
        if info:
            if sample:
                tmp = info.get(sample)
                if tmp: 
                    return tmp.get('coverage_mrk_dup')
            else:
                return info.get('coverage_mrk_dup')

def process_single_samtools_run(bam_file,all_contigs_info,samtools_bin):
    command="%s view -F 132 %s"%(samtools_bin, bam_file)
    open_stream, process=get_output_stream_from_command(command)
    current_contig=None
    coverage=0
    duplicate=0
    sample_name, ext = os.path.splitext(bam_file)
    for line in open_stream:
        sp_line=line.strip().split()
        if current_contig!=sp_line[2] and current_contig != None:
            all_contigs_info.add_values(current_contig,coverage,duplicate,sample = sample_name)
            coverage=0
            duplicate=0
        current_contig=sp_line[2]
        if int(sp_line[3])==1:
            if int(sp_line[1]) & 1024 == 1024:
                duplicate+=1
            coverage+=1
            
    open_stream.close()

def process_single_samtools_run_with_read_group(bam_file,all_contigs_info,samtools_bin):
    command="%s view -h -F 132 %s"%(samtools_bin, bam_file)
    open_stream, process=get_output_stream_from_command(command)
    current_contig=None
    sample_name, ext = os.path.splitext(bam_file)
    read_groups={}
    try:
        for line in open_stream:
            if not line.startswith("@"):
                break
            if line.startswith("@RG"):
                sp_line = line.strip().split()
                rg_id=rg_sample=rg_library=None
                for value in sp_line:
                    if value.startswith("ID"):
                        rg_id=value[3:]
                    elif value.startswith("SM"):
                        rg_sample=value[3:]
                    elif value.startswith("LB"):
                        rg_library=value[3:]
                if rg_id:
                    if rg_sample:
                        read_groups[rg_id]=rg_sample
                    elif rg_library:
                        read_groups[rg_id]=rg_library
                    else:
                        read_groups[rg_id]=rg_id
        all_sample_coverage={}
        all_sample_duplicate={}
        for sample in read_groups.values():
            all_sample_coverage[sample]=0
            all_sample_duplicate[sample]=0
        #process the first read
        sam_record = Sam_record(line.strip())
        current_contig = sam_record.get_reference_name()
        if not sam_record.is_unmapped():
            rg_id = sam_record.get_tag("RG")
            if sam_record.is_duplicate_read():
                all_sample_duplicate[read_groups.get(rg_id)]+=1
            all_sample_coverage[read_groups.get(rg_id)]+=1
        i=1
        #process all the others
        for line in open_stream:
            i+=1
            if i%1000000==0:
                print i
            sam_record = Sam_record(line.strip())
            if current_contig != sam_record.get_reference_name() and current_contig != None:
                for sample in read_groups.values():
                    all_contigs_info.add_values(current_contig, all_sample_coverage.get(sample),
                                                all_sample_duplicate.get(sample), sample = sample)
                    all_sample_coverage[sample]=0
                    all_sample_duplicate[sample]=0
            current_contig = sam_record.get_reference_name()
            
            if not sam_record.is_unmapped():
                rg_id = sam_record.get_tag("RG")
                if sam_record.is_duplicate_read():
                    all_sample_duplicate[read_groups.get(rg_id)]+=1
                all_sample_coverage[read_groups.get(rg_id)]+=1
    finally:
        open_stream.close()


def RAD_median_coverage(bam_files,output_file, with_rg=False):
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
    
    if with_rg:
        function = process_single_samtools_run_with_read_group
    else:
        function = process_single_samtools_run
        
    for bam_file in bam_files:
        sample_name, ext = os.path.splitext(bam_file)
        t = FuncThread(function, bam_file, all_contigs_info, samtools_bin)
        t.start()
        all_threads.append(t)
    for t in all_threads:
        t.join()
    all_contigs = all_contigs_info.get_all_contigs()
    all_contigs.sort()
    sample_names=all_contigs_info.get_all_samples()
    sample_names.sort()
    open_output=utils_logging.open_output_file(output_file)
    open_output.write("#contig\tcoverage\tcoverage_mrk_dup\tnb_sample")
    for sample in sample_names:
        open_output.write("\t%s\t%s_mrk_dup"%(sample,sample))
    open_output.write("\n")
    for contig in all_contigs:
        out=[]
        nb_sample=0
        coverage=all_contigs_info.get_contig_coverage_value(contig)
        coverage_mrk_dup=all_contigs_info.get_contig_coverage_mrk_dup_value(contig)
        out.append("%s\t%s\t%s"%(contig,coverage,coverage_mrk_dup))
        for sample in sample_names:
            coverage_mrk_dup=all_contigs_info.get_contig_coverage_mrk_dup_value(contig,sample)
            if coverage_mrk_dup and coverage_mrk_dup>2:
                nb_sample+=1
        out.append("%s"%(nb_sample))
        for sample in sample_names:
            coverage = all_contigs_info.get_contig_coverage_value(contig,sample)
            coverage_mrk_dup=all_contigs_info.get_contig_coverage_mrk_dup_value(contig,sample)
            if not coverage:
                coverage=0
            if not coverage_mrk_dup:
                coverage_mrk_dup=0
            out.append("%s\t%s"%(coverage,coverage_mrk_dup))
            
        open_output.write("%s\n"%('\t'.join(out)))
        
    
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
    RAD_median_coverage(bam_files,options.output_file, with_rg=options.with_read_group)
    

def _prepare_optparser():
    """Prepare optparser object. New options will be added in this
    function first.
    """
    usage = """usage: %prog <-b bam_file> [ -o output_file]"""
    description = """This script will take aligned RAD read to the consensus and calculate per consensus coverage."""
    
    optparser = OptionParser(version="None",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="show this help message and exit.")
    optparser.add_option("-b","--bam_files",dest="bam_files",type="string",
                         help="Path to one or several bam files from which the coverage should be extracted. Default: %default")
    optparser.add_option("-o","--output_file",dest="output_file",type="string",
                         help="The path to the file where the data should be output. If not set, the results will be print to stdout. Default: %default")
    optparser.add_option("-r","--with_read_group",dest="with_read_group",action="store_true",default=False,
                         help="Make the script use the read group of the bam files to determine samples instead of the file name. Default: %default")
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
    
