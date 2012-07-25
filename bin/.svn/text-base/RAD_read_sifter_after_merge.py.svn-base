#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 20 October 2011
@author: tcezard
'''
import sys,re
from utils import  utils_logging
import command_runner
from utils.FastaFormat import FastaReader
from IO_interface.samIterator import Sam_record
import utils
import logging
from optparse import OptionParser


def load_new_fasta(fasta_file):
    all_sequences={}
    open_file = open(fasta_file)
    fasta_reader = FastaReader(open_file)
    for fasta_record in fasta_reader:
        header,sequence = fasta_record
        name='_'.join(header.split('_')[:2])
        all_sequences[name]=(header,sequence)
    open_file.close()
    return all_sequences

def shift_reads(bam_file, fasta_file, output_sam_file):
    all_sequences=load_new_fasta(fasta_file)
    stream = utils.get_sam_stream(bam_file, options='-h')
    open_output= utils_logging.open_output_file(output_sam_file, pipe=True)
    open_output.write("@HD\tVN:1.0\tSO:unsorted\n")
    all_values = all_sequences.values()
    all_values.sort(key=lambda x: x[0])
    for header,sequence in all_values:
        open_output.write("@SQ\tSN:%s\tLN:%s\n"%(header,len(sequence)))
    #read the header to get the read groups
    for line in stream:
        if line.startswith('@'):
            if line.startswith('@RG'):
                open_output.write("%s\n"%(line.strip()))
        else:break
    
    sam_record=Sam_record(line)
    sam_record = process_one_record(sam_record, all_sequences)
    if sam_record:
        open_output.write(str(sam_record))
    for line in stream:
        sam_record=Sam_record(line)
        sam_record = process_one_record(sam_record, all_sequences)
        if sam_record:
            open_output.write(str(sam_record))
    open_output.close()
    
def process_one_record(sam_record, all_sequences):
    reference_name = sam_record.get_reference_name()
    name='_'.join(reference_name.split('_')[:2])
    fasta_record = all_sequences.get(name)
    discard_read=False
    if fasta_record:
        header, sequence = fasta_record
        to_add = get_addition_from_header(header)
        if sam_record.is_first_read():
            sam_record.set_reference_name(header)
            sam_record.set_mate_reference_name('=')
            if not sam_record.is_unmapped() and not sam_record.is_mate_unmapped():
                mate_position=sam_record.get_mate_pos()+to_add
                sam_record.set_mate_pos(mate_position)
                sam_record.set_insert_size(mate_position-sam_record.get_position())
        if sam_record.is_second_read():
            sam_record.set_reference_name(header)
            sam_record.set_mate_reference_name('=')
            if not sam_record.is_unmapped():
                position=sam_record.get_position()+to_add
                if position<1:
                    discard_read=True
                sam_record.set_position(position)
                sam_record.set_insert_size(sam_record.get_mate_pos()-position)
    else:
        discard_read=True
    if discard_read:
        return None
    else:
        return sam_record
            
            
def get_addition_from_header(header):
    sp_header = header.split('_')
    offset=int(sp_header[-2].strip('os'))
    cigar=sp_header[-1]
    all_cigar = re.findall('(\d+)([MXDI])', cigar)
    (count, type) =all_cigar[0]
    to_add=0
    if offset>0:
        to_add=-offset
    elif type=='D':
        to_add=int(count)
    return to_add

        
                        
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
    command_runner.set_command_to_run_localy()
    shift_reads(options.bam_file,options.fasta_file, options.output_file)

def _prepare_optparser():
    """Prepare optparser object. New options will be added in this
    function first.
    """
    usage = """usage: %prog <-f fast_file>"""
    description = """This script merge read1 and read2 contigs from RAD data."""
    
    optparser = OptionParser(version="None",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="show this help message and exit.")
    optparser.add_option("-f","--fasta_file",dest="fasta_file",type="string",
                         help="Path to the fasta file where the consensus are stored. Default: %default")
    optparser.add_option("-b","--bam_file",dest="bam_file",type="string",
                         help="Path to the bam file where the reads are stored. Default: %default")
    optparser.add_option("-o","--output_file",dest="output_file",type="string",
                         help="Path to the output file where the sam information. Default: %default")
    
    return optparser


def _verifyOption(options):
    """Check if the mandatory option are present in the options objects.
    @return False if any argument is wrong."""
    arg_pass=True
    
    if not options.fasta_file :
        logging.error("You must specify a fasta file.")
        arg_pass=False
    return arg_pass



if __name__=="__main__":
    main()
