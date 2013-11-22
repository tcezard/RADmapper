#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 20 October 2011
@author: tcezard
'''
import sys, os
from utils import  utils_logging
import logging
import command_runner
from utils.FastaFormat import FastaReader
from utils.utils_commands import get_output_stream_from_command
from IO_interface.samIterator import Sam_record


 
def run_merger(consenus_r1,contig_r2,merged_fa,alignment_report):
    merger_bin="merger"
    command="%s -asequence %s -bsequence %s -outfile %s -outseq %s"%(merger_bin, consenus_r1,contig_r2, alignment_report, merged_fa)
    return_code = command_runner.run_command(command)
    return return_code


def run_revseq(fasta_file,output_fastq):
    revseq_bin="revseq"
    command="%s -sequence %s -outseq %s"%(revseq_bin, fasta_file,output_fastq)
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



def merge_2_contigs(consensus_name, contig_file1,contig_file2, output_dir):
    contig_file2_revcomp = "%s_rev_com.fa"%(contig_file2)
    return_code = run_revseq(contig_file2,contig_file2_revcomp)

    merged_fa = os.path.join(output_dir,"%s_putativemerged.fa"%(consensus_name))
    alignment_report= os.path.join(output_dir,"%s_alignment_report.txt"%(consensus_name))
    results=None

    return_code = run_merger(contig_file1,contig_file2, merged_fa, alignment_report)

    matches1=mismatches1=offset1=matches2=mismatches2=offset2=0
    cigar1=""
    merged_sequence1=""
    cigar2=""
    merged_sequence2=""

    return_code = run_merger(contig_file1,contig_file2, merged_fa, alignment_report)
    if return_code==0:
        (matches1, mismatches1, offset1,
         cigar1, merged_sequence1)=parse_alignment_report(alignment_report)

    return_code = run_merger(contig_file1,contig_file2_revcomp, merged_fa, alignment_report)
    if return_code==0:
        (matches2, mismatches2, offset2,
         cigar2, merged_sequence2)=parse_alignment_report(alignment_report)

    if float(mismatches1)/(matches1+mismatches1) < float(mismatches2)/(matches2+mismatches2):
        matches=matches1
        mismatches=mismatches1
        offset=offset1
        cigar=cigar1
        merged_sequence=merged_sequence1
        orientation="forward"
    else:
        matches=matches2
        mismatches=mismatches2
        offset=offset2
        cigar=cigar2
        merged_sequence=merged_sequence2
        orientation="reverse"

    if matches+mismatches>15 and float(mismatches)/(matches+mismatches)<0.1:
        logging.info("Merging %s successful"%orientation)
        result_name="%s_merged_os%s_%s"%(consensus_name, offset, ''.join(cigar))
        result_seq = merged_sequence.upper()
        results=(result_name,result_seq)
    else:
        logging.info("Merging Failed")
        results=None


    if os.path.exists(merged_fa):
        os.remove(merged_fa)
        os.remove(alignment_report)
    return results

        
                        

if __name__=="1__main__":
    alignment_report=sys.argv[1]
    parse_alignment_report(alignment_report)

if __name__=="1__main__":
    utils_logging.init_logging()
    command_runner.set_command_to_run_localy()
    fasta_file=sys.argv[1]
    read_consensus(fasta_file)
