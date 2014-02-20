#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 19 March  2012
@author: tcezard
'''
import sys, os
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
from utils.utils_commands import get_output_stream_from_command
from collections import Counter, defaultdict
from IO_interface.vcfIO import VcfReader
import shutil
import copy

iupac_alphabet_list=['A', 'C', 'B', 'D', 'G', 'H', 'K', 'M', 'N', 'S', 'R', 'T', 'W', 'V', 'Y', 'a', 'c', 'b', 'd', 'g', 'h', 'k', 'm', 'n', 's', 'r', 't', 'w', 'v', 'y']

def parse_RG_line(line):
    dict = {}
    sp_line = line.split('\t')
    for element in sp_line[1:]:
        tmp = element.split(':')
        dict[tmp[0]]=':'.join(tmp[1:])
    return dict


def generate_readgroup_exclusion_file_per_samples(bam_file):
    directory = os.path.dirname(os.path.abspath(bam_file))
    command = 'samtools view -H %s | grep @RG'%(bam_file)
    stream, process = get_output_stream_from_command(command)
    all_samples=set()
    all_samples2id=defaultdict(list)
    for line in stream:
        RG_dict = parse_RG_line(line)
        all_samples.add(RG_dict.get('SM'))
        all_samples2id[RG_dict.get('SM')].append(RG_dict.get('ID'))
    
    all_samples2exclusion_id_file={}
    for sample in all_samples:
        exclusion_id = []
        exclusion_samples = all_samples.difference(set([sample]))
        for exclusion_sample in exclusion_samples:
            exclusion_id.extend(all_samples2id.get(exclusion_sample))
        sample_exclusion_file=os.path.join(directory,'exclusion_id_for_%s.txt'%sample)
        
        with open(sample_exclusion_file,'w') as open_file: open_file.write('\n'.join(exclusion_id))
        all_samples2exclusion_id_file[sample]= sample_exclusion_file  
    return all_samples2exclusion_id_file
        

def calculate_base_frequency_for_snps(directory, name):
    bam_file = os.path.join(directory, name+'_corrected_sorted_mrk_dup_fixed.bam')
    vcf_file = os.path.join(directory, name+'_corrected_sorted_mrk_dup_fixed_samtools_filterd20q60.vcf')
    samples=[]
    if os.path.exists(vcf_file):
        list_snp_position=[]
        with open(vcf_file) as open_vcf:
            reader=VcfReader(open_vcf)
            list_snp_position = ['%s\t%s'%(rec.get_reference(),rec.get_position()) for rec in reader]
        if len(list_snp_position)>0:
            snp_positions_file=os.path.join(directory, 'snps_positions.txt')
            with open(snp_positions_file,'w') as open_file:open_file.write('\n'.join(list_snp_position))
            all_samples2exclusion_id_file = generate_readgroup_exclusion_file_per_samples(bam_file)
            samples=all_samples2exclusion_id_file.keys()
            for sample in all_samples2exclusion_id_file.keys():
                exclusion_id_file = all_samples2exclusion_id_file.get(sample)
                output_file=os.path.join(directory, name+'_corrected_sorted_mrk_dup_fixed_samtools_%s.allelefreq'%sample)
                allele_freq_from_bam_and_list_pos(output_file, bam_file, snp_positions_file, list_snp_position, exclusion_id_file)
    return samples

def replace_number_base(match_object):
    s=match_object.group()
    m=re.search('[0-9]+',s)
    if m:
        return s[m.end()+int(m.group()):]
    else: return ''


def get_base_frequency_from_line(line, bas_qual_threshold, map_qual_threshold, test_mapping_qual=True):
    ATCG={'A':0,'T':0,'C':0,'G':0}
    ATCG_filtered={'A':0,'T':0,'C':0,'G':0}
    
    sp_line = line.strip().split()
    ## do not process line specifying the deletion
    if sp_line[2]=='*':
        return 
    ##remove the insertion
    iupac_alphabet=''.join(iupac_alphabet_list)
    bases=re.sub('\+[0-9]+['+iupac_alphabet+']+',replace_number_base,sp_line[4])
    ##remove the deletion 
    bases=re.sub('\-[0-9]+['+iupac_alphabet+']+',replace_number_base,bases)
    ##remove weird character for start and end of the read
    bases=re.sub('\^.','',bases)
    bases=re.sub('\$','',bases)
    sp_bases=re.findall('['+iupac_alphabet+'.,*><]', bases) # sp_line[8] is the bases
    if len(sp_line)<=6:
        test_mapping_qual=False
    if len(sp_bases)==len(sp_line[5]) and (not test_mapping_qual or len(sp_bases)==len(sp_line[6])):
        
        for i, base in enumerate(sp_bases):
            bas_qual=ord(sp_line[5][i])-33
            if test_mapping_qual:
                map_qual=ord(sp_line[6][i])-33
            else:
                map_qual=40
            if bas_qual>=bas_qual_threshold and map_qual>=map_qual_threshold:
                if base=='.' or base==',':
                    base=sp_line[2]
                if base.upper() in ATCG.keys():
                    ATCG[base.upper()]+=1
            else:
                if base=='.' or base==',':
                    base=sp_line[2]
                if base.upper() in ATCG.keys():
                    ATCG_filtered[base.upper()]+=1
    else:
        print 'problem in line %s'%line
        print '%s (%s) and %s (%s) and %s (%s) have different length'%(''.join(sp_bases), len(sp_bases) ,sp_line[-2],
                                                                       len(sp_line[-2]), sp_line[-1], len(sp_line[-1]))
        return None
    return ATCG,ATCG_filtered



def format_base_freq_line(chr, position, ref_base, consensus_base,  ATCG, ATCG_filtered, coverage):
    out=[]
    out.append(chr)
    out.append('%s'%position)
    out.append(ref_base)
    out.append(consensus_base)
    out.append('%s'%(coverage))
    out.append('A:%s:%s'%(ATCG['A'],ATCG_filtered['A']))
    out.append('T:%s:%s'%(ATCG['T'],ATCG_filtered['T']))
    out.append('C:%s:%s'%(ATCG['C'],ATCG_filtered['C']))
    out.append('G:%s:%s'%(ATCG['G'],ATCG_filtered['G']))
    return '\t'.join(out)

def get_mpileup_from_bam(bam_file, options=''):
    try:
        pipeline_parm=utils_param.get_pipeline_parameters()
        samtools_bin=os.path.join(pipeline_parm.get_samtools_dir(),'samtools')
    except Config_file_error, e:
        logging.warning("Can't find the configuration file you'll need to have samtools in you path.")
        samtools_bin='samtools'
    if bam_file=='PIPE':
        bam_file='-'
    else:
        command = '%s mpileup -A %s %s'%(samtools_bin, bam_file, options)
    stream, process = utils_commands.get_output_stream_from_command(command, logger_name=None)
    return stream


def allele_freq_from_bam_and_list_pos(output_file, input_file, list_position_file, all_positions_loaded, exclusion_id_file, bas_qual_threshold=20,
                                   map_qual_threshold=10, coverage_threshold=6):
    input_stream = get_mpileup_from_bam(input_file, options='-s -l %s -G %s'%(list_position_file,exclusion_id_file))
    all_positions_loaded=copy.copy(all_positions_loaded)
    if input_stream is not None:
        open_output=open(output_file,'w')
        for line in input_stream:
            sp_line = line.strip().split()
            position = '%s\t%s'%(sp_line[0],sp_line[1])
            
            if position in all_positions_loaded :
                all_positions_loaded.remove(position)
            else:
                continue
            
            info=get_base_frequency_from_line(line, bas_qual_threshold, map_qual_threshold)
            if info:
                ATCG,ATCG_filtered=info
            else:
                ATCG={'A':0,'T':0,'C':0,'G':0}
                ATCG_filtered={'A':0,'T':0,'C':0,'G':0}
            
            ##Calculate overall coverage
            coverage=ATCG['A']+ATCG['C']+ATCG['G']+ATCG['T']
            #if coverage<coverage_threshold:
            #    continue
            
            ## Get the most frequent base
            sorted_list=sorted(ATCG, key=lambda x: ATCG[x], reverse=True)
            
            if sorted_list[0]!=sp_line[2].upper():
                snp_base=sorted_list[0]
            else:
                snp_base=sorted_list[1]
            
            output_line=format_base_freq_line(chr=sp_line[0], position=sp_line[1], ref_base=sp_line[2].upper(),
                                  consensus_base=snp_base, ATCG=ATCG, ATCG_filtered=ATCG_filtered,
                                  coverage=coverage)
            open_output.write('%s\n'%(output_line))
        input_stream.close()
        for position in all_positions_loaded:
            reference, coordinate = position.split('\t')
            ATCG={'A':0,'T':0,'C':0,'G':0}
            ATCG_filtered={'A':0,'T':0,'C':0,'G':0}
            output_line=format_base_freq_line(chr=reference, position=coordinate, ref_base='N',
                                  consensus_base='N', ATCG=ATCG, ATCG_filtered=ATCG_filtered,
                                  coverage=0)
            open_output.write('%s\n'%(output_line))
        open_output.close()


def run_all_fastq_files(directory):
    directory=os.path.abspath(directory)
    all_dirs = glob(os.path.join(directory,'*_dir'))
    all_samples=set()
    for sub_dir in all_dirs:
        print sub_dir
        name=os.path.basename(sub_dir)[:-len("_dir")]
        samples=calculate_base_frequency_for_snps(sub_dir, name)
        all_samples.update(set(samples))
        
    for sample in all_samples:
        #concatenate the allele frequency file per samples
        merged_file = os.path.join(directory,'samtools_snps_%s.allelefreq'%sample)
        command = 'cat %s/*_dir/*_%s.allelefreq > %s'%(directory, sample, merged_file)
        command_runner.run_command(command)
    return

def main():
    #initialize the logging
    utils_logging.init_logging(logging.INFO)
    #utils_logging.init_logging(logging.CRITICAL)
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
    if not options.print_command:
        command_runner.set_command_to_run_localy()
    run_all_fastq_files(options.consensus_dir)
    
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
