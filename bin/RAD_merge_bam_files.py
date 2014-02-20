#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 19 March  2012
@author: tcezard
'''
import sys, os
import pysam
from utils import utils_logging, utils_commands, longest_common_substr_from_start, utils_param, sort_bam_file_per_coordinate
import logging, threading, re
from optparse import OptionParser
from glob import glob
import command_runner
from utils.FastaFormat import FastaReader
import time
from RAD_merge_read1_and_read2 import merge_2_contigs
from utils.utils_commands import get_output_stream_from_command
from collections import defaultdict, OrderedDict, Counter


#get the path to the current script to infer the path to RAD_set_read1_consensus_to_read2.py
RADmapper_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
path_to_picard = os.path.join(RADmapper_dir,"picard")
mergeSamFilesWithCat_jar=os.path.join(path_to_picard,'MergeSamFilesWithCat.jar')

def get_readgroup_and_seq_dict_from_bam(bam_files, allow_collision=False):
    print "Gather seq dict and read groups from %s bam files"%len(bam_files)
    all_read_groups={}
    all_seq_dict=OrderedDict()
    for bam_file in bam_files:
        command = "samtools view -H %s | egrep '@RG|@SQ' "%bam_file
        stdout, process = utils_commands.get_output_stream_from_command(command)
        for line in stdout:
            if line.startswith('@RG'):
                read_group_dict={}
                for element in line.strip().split('\t'):
                    if element != '@RG':
                        key,value = element.split(':')
                        read_group_dict[key]=value
                if read_group_dict.has_key('ID') and read_group_dict.get('ID') not in all_read_groups:
                    all_read_groups[read_group_dict.get('ID')]=read_group_dict
            if line.startswith('@SQ'):
                seq_dict={}
                for element in line.strip().split('\t'):
                    if element != '@SQ':
                        key,value = element.split(':')
                        if key=='LN':value=int(value)
                        seq_dict[key]=value
                if seq_dict.has_key('SN'):
                    name= seq_dict.get('SN')
                    if all_seq_dict.has_key(name) and not allow_collision:
                        raise StandardError("identical sequence dictionary name %s in two different entry and collision not allowed"%name)
                    all_seq_dict[name]=seq_dict

    return all_read_groups.values(), all_seq_dict.values()


def create_sequence_dict_from_contigs_file(contig_file):
    print "Read %s"%contig_file
    sequence_dictionary=[]
    all_names={}
    with open(contig_file) as open_file:
        for header,sequence in FastaReader(open_file):
            sequence_dictionary.append({"SN":header,"LN":len(sequence)})
            if all_names.has_key(header):
                raise StandardError("Duplicated reference name %s in %s"%(header,contig_file))
            all_names[header]=1
    return sequence_dictionary

def create_read_group_dict_from_read_group_text(readgroups):
    print "Read %s"%readgroups
    all_read_group_dicts=[]
    if isinstance(readgroups, list):
        open_rg = readgroups
    else:
        open_rg = open(readgroups)
    for read_group in open_rg:
        read_group=read_group.strip()
        read_group_dict={}
        for element in read_group.split('\t'):
            if element != '@RG':
                key,value = element.split(':')
                read_group_dict[key]=value
        if read_group_dict.has_key('ID'):
            all_read_group_dicts.append(read_group_dict)
    if isinstance(readgroups, str):
        open_rg.close()
    return all_read_group_dicts

def get_per_sample_read_group(read_group_dicts):
    per_sample_read_group=defaultdict(list)
    rgid2sample={}
    for read_group in read_group_dicts:
        if read_group.has_key('SM'):
            sample = read_group.get('SM')
            per_sample_read_group[sample].append(read_group)
            rgid2sample[read_group.get('ID')]=sample
        else:
            per_sample_read_group["no_name"].append(read_group)
            rgid2sample[read_group.get('ID')]="no_name"
    return per_sample_read_group, rgid2sample

def merge_bam_file_per_samples(all_bam_files, read_group_dicts,sequence_dictionary, output_bam_file):
    per_sample_open_bam_file={}
    read_group_dicts_per_sample,rgid2sample = get_per_sample_read_group(read_group_dicts)
    nb_record_per_sample=Counter()
    for sample in read_group_dicts_per_sample:
        output_bam_file_name,ext=os.path.splitext(output_bam_file)
        output_bam_file_for_sample=output_bam_file_name+'_%s.bam'%sample
        header={'HD': {'VN': '1.0','SO':'coordinate'},
        'SQ': sequence_dictionary,
        'RG': read_group_dicts_per_sample.get(sample)}
        print "Open %s"%output_bam_file_for_sample
        per_sample_open_bam_file[sample]=pysam.Samfile( output_bam_file_for_sample, "wb", header = header)
    nb_record=0
    for bam_file in all_bam_files:
        print "add record from %s"%bam_file

        with pysam.Samfile( bam_file, "rb" ) as open_bam:
            for record in open_bam.fetch( until_eof = True ):

                #Get the sample bam file associated with this read group
                sample=rgid2sample.get(record.opt('RG'))
                nb_record_per_sample[sample]+=1
                final_bam=per_sample_open_bam_file.get(sample)
                old_tid =  record.tid
                mate_tid = record.rnext

                rname= open_bam.getrname(record.tid)
                new_tid = final_bam.gettid(rname)
                record.tid=new_tid
                if record.is_paired and not record.mate_is_unmapped:
                    mate_rname = open_bam.getrname(record.rnext)
                    new_mate_tid= final_bam.gettid(mate_rname)
                    record.rnext=new_mate_tid

                final_bam.write(record)
                nb_record+=1
                if nb_record%1000000==0:
                    print "%s record added -- %s"%(nb_record, ', '.join(['%s:%s'%(s,n) for s ,n in nb_record_per_sample.iteritems()]))
        print "%s record added -- %s"%(nb_record, ', '.join(['%s:%s'%(s,n) for s ,n in nb_record_per_sample.iteritems()]))
    for sample in per_sample_open_bam_file:
        per_sample_open_bam_file.get(sample).close()

def merge_bam_files(all_bam_files, output_bam_file, contig_file=None, read_groups_file=None, per_sample=False):
    if contig_file and read_groups_file:
        read_group_dicts=create_read_group_dict_from_read_group_text(read_groups_file)
        sequence_dictionary = create_sequence_dict_from_contigs_file(contig_file)
    else:
        read_group_dicts,sequence_dictionary= get_readgroup_and_seq_dict_from_bam(all_bam_files)
    if per_sample:
        merge_bam_file_per_samples(all_bam_files, read_group_dicts,sequence_dictionary, output_bam_file)
    else:
        header={'HD': {'VN': '1.0','SO':'coordinate'},
                'SQ': sequence_dictionary,
                'RG': read_group_dicts}
        print "Created Header"
        with pysam.Samfile( output_bam_file, "wb", header = header) as final_bam:
            for bam_file in all_bam_files:
                print "add record from %s"%bam_file
                nb_record=0
                with pysam.Samfile( bam_file, "rb" ) as open_bam:
                    for record in open_bam.fetch( until_eof = True ):
                        nb_record+=1
                        old_tid =  record.tid
                        mate_tid = record.rnext

                        rname= open_bam.getrname(record.tid)
                        new_tid = final_bam.gettid(rname)
                        record.tid=new_tid
                        if record.is_paired and not record.mate_is_unmapped:
                            mate_rname = open_bam.getrname(record.rnext)
                            new_mate_tid= final_bam.gettid(mate_rname)
                            record.rnext=new_mate_tid

                        final_bam.write(record)
                print "%s record added"%nb_record
    return 0



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
    bam_files=[options.bam_files]
    if len(args)>0:
        bam_files.extend(args)
    code=merge_bam_files(bam_files, options.output_file, options.genome_file, options.read_group_file,options.per_sample)
    sys.exit(code)

def _prepare_optparser():
    """Prepare optparser object. New options will be added in this
    function first.
    """
    usage = """usage: %prog <-b bam_file> [ -o output_file]"""
    description = """This script will merge bam file that have completely independent reference into one. The bam file are assumed to be coordinate sorted"""
    
    optparser = OptionParser(version="None",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="show this help message and exit.")
    optparser.add_option("-b","--bam_files",dest="bam_files",type="string",
                         help="Path to all_bam files to merge. Default: %default")
    optparser.add_option("-o","--output_file",dest="output_file",type="string",
                         help="Path to the final merged file. Default: %default")
    optparser.add_option("-g","--genome_file",dest="genome_file",type="string",
                         help="Path to the file containing all contigs. Default: %default")
    optparser.add_option("-r","--read_group_file",dest="read_group_file",type="string",
                         help="Path to the file containing all read groups. Default: %default")
    optparser.add_option("--per_sample",dest="per_sample",action='store_true',default=False,
                         help="Set the script to write the bam file per sample. Default: %default")
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
    
    
