#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 20 October 2011
@author: tcezard
'''
import sys, os, re
from utils import  utils_logging
import logging, threading
from optparse import OptionParser
import command_runner
from utils.FastaFormat import FastaReader
from utils.utils_commands import get_output_stream_from_command
from IO_interface.samIterator import Sam_record


 
def run_merger(consenus_r1,contig_r2,merged_fa,alignment_report):
    merger_bin="merger"
    command="%s -asequence %s -bsequence %s -outfile %s -outseq %s"%(merger_bin, consenus_r1,contig_r2, alignment_report, merged_fa)
    return_code = command_runner.run_command(command)
    return return_code


def parse_alignment_report(alignment_report):
    open_file = open(alignment_report)
    position=1
    read1_elements=[]
    align_elements=[]
    read2_elements=[]
    for line in open_file:
        #skip blank lines
        if len(line)<2:
            continue
        if line.startswith('#'):
            continue
        align_line=line[21:]
        
        if position==1:
            read1_elements.append(align_line.split()[0])
            position+=1
        elif position==2:
            align_elements.append(align_line[:-1])
            position+=1
        elif position==3:
            read2_elements.append(align_line.split()[0])
            position=1
    read1=''.join(read1_elements)
    align=''.join(align_elements)
    read2=''.join(read2_elements)
    cigar_elements=[]
    cigar_entry=[]
    curr_letter=""
    for i in range(len(read1)):
        if read1[i]=='-':
            letter='I'
        elif align[i]=='|':
            letter='M'
        elif align[i]=='.':
            letter='X'
        elif read2[i]=='-':
            letter='D'
        if curr_letter!=letter and len(cigar_entry)>0:
            cigar_elements.append(''.join(cigar_entry))
            cigar_entry=[]
        cigar_entry.append(letter)
        curr_letter=letter
    cigar_elements.append(''.join(cigar_entry))
    offset=0
    if cigar_elements[0].startswith('I'):
        offset = len(cigar_elements[0])
        cigar_elements=cigar_elements[1:]
    merged_sequence=[read1.replace('-','')]
    position=offset
    matches=0
    mismatches=0
    for element in cigar_elements:
        if element[0]=="M":
            position+=len(element)
            matches+=len(element)
        elif element[0]=="X":
            position+=len(element)
            mismatches+=len(element)
        elif element[0]=="D":
            position+=len(element)
        else:
            merged_sequence.append(read2[position:position+len(element)])
            position+=len(element)
    #print offset
    #print read1
    #print align
    #print read2
    #print ' '*offset+''.join(cigar_elements)
    #print ' '*offset+''.join(merged_sequence)
    return (matches, mismatches, offset, cigar_collapse(''.join(cigar_elements)), ''.join(merged_sequence) )

def cigar_collapse(cigar):
    curr_c=""
    count=0
    cigar_elements=[]
    for c in cigar:
        if c == curr_c:
            count+=1
        else:
            if count>0:
                cigar_elements.append('%d%s'%(count,curr_c))
            curr_c=c
            count=1
    cigar_elements.append('%d%s'%(count,curr_c))
    return  cigar_elements


def read_consensus(fasta_file, regex="consensus_(\d+)"):
    logging.info("load fasta file")
    all_fasta_entries={}
    open_fasta = open(fasta_file)
    fasta_reader = FastaReader(open_fasta)
    pattern = re.compile(regex)
    for header,sequence in fasta_reader:
        match = pattern.match(header)
        if not match:
            logging.warning("%s does'nt match the provided regex %s"%(header,regex))
        else:
            number = match.group(1)
            list_fasta = all_fasta_entries.get(number) 
            if not list_fasta:
                list_fasta=[]
                all_fasta_entries[number]=list_fasta
            list_fasta.append((header,sequence))

    logging.info("run emboss merger")
    for key in all_fasta_entries.keys():
        list_fasta=all_fasta_entries.get(key)
        if len(list_fasta)==2:
            header1,sequence1=list_fasta[0]
            header2,sequence2=list_fasta[1]
            if len(sequence1)<=len(sequence2):
                match = re.match(regex,header1)
                if match:
                    read1_name='consensus_%s_read1'%match.group(1)
                    read2_name='consensus_%s_read2'%match.group(1)
                read1_sequence=sequence1
                read2_sequence=sequence2
            else:
                match = re.match(regex,header1)
                if match:
                    read2_name='consensus_%s_read1'%match.group(1)
                    read1_name='consensus_%s_read2'%match.group(1)
                read1_sequence=sequence2
                read2_sequence=sequence1
            output_file = open(read1_name,'w')
            output_file.write(">"+read1_name+'\n')
            output_file.write(read1_sequence+'\n')
            output_file.close()
            
            output_file = open(read2_name,'w')
            output_file.write(">"+read2_name+'\n')
            output_file.write(read2_sequence+'\n')
            output_file.close()

            
            merge_2_contigs("consensus_%s"%key, read1_name , read2_name, output_dir='.')
            os.remove(read1_name)
            os.remove(read2_name)

def merge_2_contigs(consensus_name, contig_file1,contig_file2, output_dir):
    merged_fa = os.path.join(output_dir,"%s_putativemerged.fa"%(consensus_name))
    alignment_report= os.path.join(output_dir,"%s_alignment_report.txt"%(consensus_name))
    results_file=os.path.join(output_dir,"%s_merged.fa"%(consensus_name))
    
    return_code = run_merger(contig_file1,contig_file2, merged_fa, alignment_report)
    if return_code==0:
        (matches, mismatches, offset,
         cigar, merged_sequence)=parse_alignment_report(alignment_report)
        if matches+mismatches>10 and float(mismatches)/(matches+mismatches)<0.1:
            logging.info("Merging successful")
            output_file = open(results_file,'w')
            output_file.write(">%s_merged_os%s_%s\n"%(consensus_name, offset, ''.join(cigar)))
            output_file.write(merged_sequence.upper()+"\n")
            output_file.close()
        else:
            logging.info("Merging Failed")
            results_file=None
    else:
        logging.info("merger crashed with return code=%s"%(return_code))
        results_file=None
    if os.path.exists(merged_fa):
        os.remove(merged_fa)
    #os.remove(alignment_report)
    return results_file

def shift_read_2_reads(offset, cigar, bam_file):
    command= 'samtools view -h %s'%(bam_file)
    stream, process = get_output_stream_from_command(command)
    for line in stream:
        if line.startswith("@"):
            continue
        sam_record = Sam_record(line)
        if sam_record.is_first_read():
            pass
            #print sam_record
        else:
            print sam_record
        
                        
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
    read_consensus(options.fasta_file)

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

if __name__=="1__main__":
    alignment_report=sys.argv[1]
    parse_alignment_report(alignment_report)

if __name__=="1__main__":
    utils_logging.init_logging()
    command_runner.set_command_to_run_localy()
    fasta_file=sys.argv[1]
    read_consensus(fasta_file)
