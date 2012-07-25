#!/usr/bin/env python
'''
Created on 4 Feb 2010

@author: tcezard
'''
import os, sys, logging
from optparse import OptionParser
import utils
from utils import utils_param, utils_logging, utils_commands
from utils.parameters import Config_file_error
from IO_interface.samIterator import Sam_record
import itertools
import command_runner
import math

"""Method that is used in the computation of estimated library size.
    extracted from picard markduplicate.
    """
def f(x,  c,  n):
    return c/x - 1 + math.exp(-n/x);

"""Estimates the size of a library based on the number of paired end molecules observed
   and the number of unique pairs ovserved.
     
   Based on the Lander-Waterman equation that states:
       C/X = 1 - exp( -N/X )
   where
       X = number of distinct molecules in library
       N = number of read pairs
       C = number of distinct fragments observed in read pairs

extracted from picard markduplicate.
     """
def estimate_library_size(nb_read_pairs, nb_unique_read_pairs):
    read_pair_duplicates = nb_read_pairs - nb_unique_read_pairs

    if nb_read_pairs > 0 and read_pair_duplicates > 0:
        n = nb_read_pairs
        c = nb_unique_read_pairs

        m = 1.0
        M = 100.0

        if c >= n or f(m*c, c, n) < 0:
            raise ValueError("Invalid values for pairs and unique pairs: "+ n + ", " + c)

        while f(M*c, c, n) >= 0:
            M *= 10.0

        for i in range(40):
            r = (m+M)/2.0
            u = f( r * c, c, n )
            if u == 0 :
                break
            elif u > 0: 
                m = r
            elif u < 0:
                M = r

        return c * (m+M)/2.0
    else:
        return None;

def mark_duplicate(bam_file,  distance_threshold=5):
    try:
        pipeline_param=utils_param.get_pipeline_parameters()
        samtools_dir=pipeline_param.get_samtools_dir()
        picard_dir=pipeline_param.get_picard_dir()
    except Config_file_error, e:
        logging.warning("You'll need to have samtools in your path")
        samtools_dir=''
        picard_dir=None

    total_nb_dups=0
    total_nb_uniqs=0
    nb_fragment=0    
    samtools_bin=os.path.join(samtools_dir,'samtools')
    tmp,ext = os.path.splitext(bam_file)
    
    command = "%s view -h %s"%(samtools_bin, bam_file)
    input_stream, process = utils_commands.get_output_stream_from_command(command)
    tmp_bam_file = output_bam_file = tmp + '_mrk_dup.bam.tmp'
    command = "%s view -bS - >%s"%(samtools_bin, tmp_bam_file)
    output_stream, process = utils_commands.get_input_stream_from_command(command)
    nb_reference=0
    current_reference=None
    first_reads={}
    second_reads={}
    for line in input_stream:
        if line.startswith("@"):
            output_stream.write(line)
            continue
        sam_record = Sam_record(line)
        if sam_record.get_reference_name()!=current_reference and not current_reference is None:
            #process this consensus
            if current_reference!='*':
                
                nb_dups, nb_uniq = find_duplicates(first_reads,second_reads, distance_threshold)
                total_nb_uniqs+=nb_uniq
                total_nb_dups+=nb_dups
                nb_fragment+=len(second_reads)
            output_reads(output_stream, first_reads, second_reads)
            first_reads={}
            second_reads={}
        if sam_record.is_second_read():
            second_reads[sam_record.get_query_name()]=sam_record
        else:
            first_reads[sam_record.get_query_name()]=sam_record
        nb_reference+=1
        if nb_reference%1000==0:
            print "process %s consensus"%nb_reference
        current_reference = sam_record.get_reference_name()
    if sam_record.get_reference_name()!=current_reference and not current_reference is None:
        #process this consensus
        if current_reference!='*':
            nb_dups = find_duplicates(first_reads,second_reads, distance_threshold)
            total_nb_dups+=nb_dups
            nb_fragment+=len(second_reads)
        output_reads(output_stream, first_reads, second_reads)
    library_size = estimate_library_size(nb_fragment, total_nb_uniqs)
    print "%s fragments"%(nb_fragment)
    print "%s (%.2f%%) duplicates"%(total_nb_dups,float(total_nb_dups)/nb_fragment*100)
    print "nb unique=%d"%(total_nb_uniqs)
    print "library size=%d"%round(library_size,0)
    print "Sort the new bam file"
    output_stream.flush()
    output_stream.close()
    if picard_dir:
        output_bam_file = tmp + '_mrk_dup.bam'
        return_code = utils.sort_bam_file_per_coordinate(picard_dir, tmp_bam_file, output_bam_file, overwrite=True, validation_stringency="SILENT")
    else:
        output_bam_file = tmp + '_mrk_dup'
        command='%s sort %s %s'%(samtools_bin, tmp_bam_file, output_bam_file)
        return_code = command_runner.run_command( command)
    if return_code==0:
        command_runner.run_command( 'rm -f %s'%(tmp_bam_file))

def find_duplicates(first_reads,second_reads, distance_threshold):
    uniq_second_sequences={}
    all_second_reads=second_reads.values()
    nb_duplicate=0
    if len(all_second_reads)>0:
        uniq_second_sequences[all_second_reads[0].get_query_sequence()]=[all_second_reads[0]]
        
        for sam_record in all_second_reads[1:]:
            uniq=True
            for unique in uniq_second_sequences.keys():
                if hamming1(unique,sam_record.get_query_sequence()) <= distance_threshold:
                    uniq=False
                    break
            if uniq:
                uniq_second_sequences[sam_record.get_query_sequence()]=[sam_record]
            else:
                uniq_second_sequences[unique].append(sam_record)
        
        for duplicated_second_sequences in uniq_second_sequences.values():
            if len(duplicated_second_sequences) > 1:
                duplicated_second_sequences.sort(cmp=compare_quality, reverse=True)
                for sequence2 in duplicated_second_sequences[1:]:
                    sequence2.set_duplicate_flag()
                    read1 = first_reads.get(sequence2.get_query_name())
                    if read1:
                        read1.set_duplicate_flag()
                        nb_duplicate+=1
                    else:
                        logging.warning("%s has a read 2 but no read1."%sequence2.get_query_name())
    return nb_duplicate, len(uniq_second_sequences)

def output_reads(output_stream, first_reads, second_reads):
    for sam_record in second_reads.values():
        output_stream.write(str(sam_record))
    for sam_record in first_reads.values():
        output_stream.write(str(sam_record))
    

def compare_quality(sam_rec1,sam_rec2):
    return sum(itertools.imap(ord,sam_rec1.get_query_quality())) - sum(itertools.imap(ord,sam_rec2.get_query_quality()))
    
def hamming1(str1, str2):
    return sum(itertools.imap(str.__ne__, str1, str2))



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
    mark_duplicate(options.input_bam)
        
def _prepare_optparser():
    """Prepare optparser object. New options will be added in this
    function first.
    """
    usage = """usage: %prog -b <bam_file>"""
    description = """This script will remove duplicate from a bam file based on the first read reference and the second read sequence."""
    
    prog_version=utils.getWtss_version()
    optparser = OptionParser(version="%prog of wtss_pipeline v"+prog_version,description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="show this help message and exit.")
    optparser.add_option("-b","--input_bam",dest="input_bam",type="string",
                         help="Path to input bam file. Default: %default")
    return optparser


def _verifyOption(options):
    """Check if the mandatory option are present in the options objects.
    @return False if any argument is wrong."""
    arg_pass=True
    
    if not options.input_bam or not os.path.exists(options.input_bam):
        logging.error("You must specify a valid input file.")
        arg_pass=False
    return arg_pass



if __name__=="__main__":
    main()
