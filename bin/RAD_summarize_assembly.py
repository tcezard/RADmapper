#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 19 March  2012
@author: tcezard
'''
import sys, os
import copy
from utils import utils_logging, utils_commands,\
    longest_common_substr_from_start, utils_param, sort_bam_file_per_coordinate
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
from RAD_coverage_read1 import RAD_median_coverage
from collections import Counter

INVARIABLE_HEADERS=['name',"nb_contigs","max_contig_length","merge_status","nb_reads","nb_reads_mapped",
                    "nb_reads_duplicates","nb_SNPs"]

def summarize_stats(name, output_file, consensus_file, read1_fastq, read2_fastq, bam_file, vcf_file):
    if not os.path.exists(consensus_file):
        return ''
    output_dict={}
    output_dict['name']=name
    number_read2_contig, longuest_contigs, is_merged = analyse_consensus_file(consensus_file)
    output_dict['nb_contigs']=str(number_read2_contig)
    output_dict['max_contig_length']=str(longuest_contigs)
    if is_merged:output_dict['merge_status']='merged'
    else:output_dict['merge_status']='concatenated'
    
    total_nb_reads, all_read_groups=count_reads_in_fastq(read1_fastq)
    if os.path.exists(bam_file):
        dummy,duplicates,mapped,properly_paired = run_stats_from_bam(bam_file)
        output_dict['nb_reads']=str(total_nb_reads*2)
        output_dict['nb_reads_mapped']=str(mapped)
        output_dict['nb_reads_duplicates']=str(duplicates)


    if os.path.exists(vcf_file):
        nb_snps = count_SNPs(vcf_file)
        output_dict['nb_SNPs']=str(nb_snps)
    headers=copy.copy(INVARIABLE_HEADERS)

    coverage_info = get_coverage_from_bam(bam_file)
    if coverage_info:
        output_dict.update(coverage_info)
        sample_names=sorted(coverage_info.keys())
        headers.extend(sample_names)
    else:
        sample_names=[]
    with open(output_file,'w') as open_output:
        open_output.write( '\t'.join(headers))
        open_output.write('\n')
        open_output.write( print_output_dict(output_dict, headers))
        open_output.write('\n')
    return output_dict,sample_names

def print_output_dict(output_dict, header_order):
    out = [output_dict.get(header,"0") for header in header_order]
    return '\t'.join(out)


def count_reads_in_fastq(fastq_file):
    command = '''awk '{if (NR%%4==1){split($1,array,"RGID:"); print array[2]}}' %s| uniq -c'''%(fastq_file)
    logging.info(command)
    stream, process = get_output_stream_from_command(command)
    total=0
    all_read_groups=Counter()
    for line in stream:
        if len(line.strip().split())==2
            count, rgid = line.strip().split()
            count=int(count)
            total+=count
            all_read_groups[rgid]
    return total, all_read_groups

        
def analyse_consensus_file(consensus_file):
    with open(consensus_file) as open_file:
        fasta_reader = FastaReader(open_file)
        number_read2_contig=0
        longuest_contigs=0
        is_merged=False
        for header, sequence in fasta_reader:
            #_merged_os0_71D21M193I or _pair_1_length_452
            sp_header = header.split('_')
            if len(sp_header)>2 and sp_header[-3] == 'merged':
                type='merged'
                number_read2_contig+=1
                if len(sequence) > longuest_contigs:
                    longuest_contigs=len(sequence)
                is_merged=True
            elif len(sp_header)>3 and sp_header[-4] == 'pair':
                type='read2_contig'
                number_read2_contig+=1
                if len(sequence) > longuest_contigs:
                    longuest_contigs=len(sequence)
            else:
                type='read1_contig'
    return number_read2_contig, longuest_contigs, is_merged
    
def run_stats_from_bam(bam_file):
    try:
        pipeline_param=utils_param.get_pipeline_parameters()
        samtools_dir=pipeline_param.get_samtools_dir()
    except Config_file_error, e:
        logging.warning("You'll need to have samtools in your path")
        samtools_dir=''
    samtools_bin=os.path.join(samtools_dir,'samtools ')
    if os.path.exists(bam_file+".stat"):
        bam_stat_file=bam_file+".stat"
    else:
        command='%s flagstat %s > %s'%(samtools_bin,bam_file,bam_file+".stat")
        return_code = command_runner.run_command(command)
        bam_stat_file=bam_file+".stat"
    return parse_stat_file(bam_stat_file)
    

def parse_stat_file(stat_file):
    open_file = utils_logging.open_input_file(stat_file,pipe=False)
    line_number=0
    total=0
    duplicates=0
    mapped=0
    properly_paired=0

    for line in open_file:
        line_number+=1
        if line_number==1:
            total = int(line.split()[0])
        elif line_number==2:
            duplicates = int(line.split()[0])
        elif line_number==3:
            mapped = int(line.split()[0])
        elif line_number==4:
            dummy = int(line.split()[0])
        elif line_number==5:
            dummy = int(line.split()[0])
        elif line_number==6:
            dummy = int(line.split()[0])
        elif line_number==7:
            properly_paired = int(line.split()[0])
        elif line_number==8:
            dummy = int(line.split()[0])
        elif line_number==9:
            dummy = int(line.split()[0])
        elif line_number==10:
            dummy = int(line.split()[0])
        elif line_number==11:
            dummy = int(line.split()[0])
    return (total,duplicates,mapped,properly_paired)


    open_file.close()


def get_coverage_from_bam(bam_file):
    name, ext=os.path.splitext(bam_file)
    output_file=name+'_coverage.txt'
    coverage_info=RAD_median_coverage([bam_file],output_file, with_rg=True)
    all_contigs = coverage_info.get_all_contigs()
    if len(all_contigs)!=1:
        logging.error("can't calculate coverage for %s --> %s read1 contigs found"%(bam_file, len(all_contigs)))
        return None
    contig=all_contigs[0]

    sample_names=coverage_info.get_all_samples()
    coverage_info_to_return={}
    for sample in sample_names:
        coverage = coverage_info.get_contig_coverage_value(contig,sample)
        coverage_mrk_dup=coverage_info.get_contig_coverage_mrk_dup_value(contig,sample)
        coverage_info_to_return[sample]=str(coverage)
        coverage_info_to_return[sample+"_mrk_dup"]=str(coverage_mrk_dup)
    return coverage_info_to_return


def count_SNPs(vcf_file):
    count_SNPs=0
    if os.path.exists(vcf_file):
        with open(vcf_file) as open_vcf:
            for line in open_vcf:
                if line.startswith('#'):
                    pass
                else:
                    count_SNPs+=1
    return count_SNPs
        

       
def run_all_summaries(directory):
    directory=os.path.abspath(directory)
    all_dirs = glob(os.path.join(directory,'*_dir'))
    summary_stats_file=os.path.join(directory,'summary_stat.txt')
    all_samples=set()
    all_output_dicts=[]
    for sub_dir in all_dirs:
        name=os.path.basename(sub_dir)[:-len("_dir")]
        consensus_file=os.path.join(sub_dir,'best_assembly.fa')
        read1_fastq=os.path.join(sub_dir,name+"_1.fastq")
        read2_fastq=os.path.join(sub_dir,name+"_2.fastq")
        bam_file = os.path.join(sub_dir,name+"_corrected_sorted_mrk_dup_fixed.bam")
        vcf_file = os.path.join(sub_dir,name+"_corrected_sorted_mrk_dup_fixed_samtools_filterd20q60.vcf")
        output_file=os.path.join(sub_dir,name+'_summary_stat.txt')

        output_dict, sample_names = summarize_stats(name,output_file,consensus_file, read1_fastq, read2_fastq, bam_file, vcf_file)
        all_output_dicts.append(output_dict)
        all_samples.update(set(sample_names))
    headers=copy.copy(INVARIABLE_HEADERS)
    headers.extend(sorted(all_samples))

    with open(summary_stats_file,'w') as open_output:
        open_output.write( '\t'.join(headers))
        open_output.write('\n')
        for output_dict in all_output_dicts:
            open_output.write( print_output_dict(output_dict, headers))
            open_output.write('\n')


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
    command_runner.set_command_to_run_localy()
    if options.debug:
        utils_logging.init_logging(logging.DEBUG)
    run_all_summaries(options.consensus_dir)

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
    
    