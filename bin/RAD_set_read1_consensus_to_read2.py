#!/usr/bin/env python
'''
Created on 4 Feb 2010

@author: tcezard
'''
import os, sys, logging
from optparse import OptionParser
import utils
from utils import utils_param, utils_logging, compare_version_number,\
    get_bwa_version, longest_common_substr_from_start, utils_commands
import command_runner
from utils.GenomeLoader import GenomeLoader
from utils.parameters import Config_file_error
from IO_interface.samIterator import Sam_record



def set_read1_consensus_to_read2_old(input_bam_file, output_bam_file):
    try:
        pipeline_param=utils_param.get_pipeline_parameters()
        samtools_dir=pipeline_param.get_samtools_dir()
    except Config_file_error, e:
        #logging.exception('Config_file_error:')
        logging.critical("You need to have the environment variable properly set to use that script")
        return False
        
    
    samtools_bin=os.path.join(samtools_dir,'samtools')
    name, ext = os.path.splitext(output_bam_file)
    if ext=='.bam':
        output_bam_file=name
    #change_consensus_on_read2
    command ="%s view -h %s "%(samtools_bin,input_bam_file)
    logging.info(command)
    input_stream,process_input = utils_commands.get_output_stream_from_command(command)
    command ="%s view -bS - | %s sort - %s"%(samtools_bin,  samtools_bin, output_bam_file)
    logging.info(command)
    output_stream,process_output= utils_commands.get_input_stream_from_command(command)
    
    #get the header
    line = input_stream.readline()
    while line.startswith("@"):
        output_stream.write(line)
        line = input_stream.readline()
    
    while line:
        read1=Sam_record(line)
        line = input_stream.readline()
        read2=Sam_record(line)
        if read1.get_query_name() == read2.get_query_name():
            if read1.is_second_read() and read2.is_first_read():
                tmp = read1
                read1=read2
                read2=tmp
            read2.set_reference_name(read1.get_reference_name())
            output_stream.write(str(read1))
            output_stream.write(str(read2))
        else:
            logging.critical("bam file is not sorted by read name")
            input_stream.close()
            output_stream.close()
            #os.remove(output_bam_file+'.bam')
            return
        line = input_stream.readline()
        
    return_code=process_input.wait()
    print return_code
    if return_code!=0:
        sys.exit(return_code)
    return_code=process_output.wait()
    output_stream.close()
    print return_code
    if return_code!=0:
        sys.exit(return_code)
    
    
def set_read1_consensus_to_read2(input_stream, output_stream):
    
    #get the header
    line = input_stream.readline()
    while line.startswith("@"):
        output_stream.write(line)
        line = input_stream.readline()
    prev_read=Sam_record(line)
    for line in input_stream:
        read=Sam_record(line)
        if prev_read and read.get_query_name() == prev_read.get_query_name():
            if read.is_second_read() and prev_read.is_first_read():
                read1=prev_read
                read2=read
            else:
                read2=prev_read
                read1=read
            read2.set_reference_name(read1.get_reference_name())
            read2.set_unmapped_flag(False)
            read2.set_position(1)
            read2.set_cigar_string("%sM"%len(read2.get_query_sequence()))
            output_stream.write(str(read1))
            output_stream.write(str(read2))
            prev_read=None
        elif prev_read:
            output_stream.write(str(prev_read))
            prev_read=read
        else:
            prev_read=read
        
    #input_stream.close()
    #output_stream.close()


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
    set_read1_consensus_to_read2(sys.stdin, sys.stdout)
   
def _prepare_optparser():
    """Prepare optparser object. New options will be added in this
    function first.
    """
    usage = """usage: %prog """
    description = """"""
    
    optparser = OptionParser(description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="show this help message and exit.")
    
    return optparser


def _verifyOption(options):
    """Check if the mandatory option are present in the options objects.
    @return False if any argument is wrong."""
    arg_pass=True
    
    return arg_pass



if __name__=="__main__":
    main()
