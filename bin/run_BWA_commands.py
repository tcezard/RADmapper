#!/usr/bin/env python
'''
Created on 4 Feb 2010

@author: tcezard
'''
import os, sys, logging
from optparse import OptionParser
import utils
from utils import utils_param, utils_logging, compare_version_number,\
    get_bwa_version, longest_common_substr_from_start
import command_runner
from utils.GenomeLoader import GenomeLoader
from utils.parameters import Config_file_error

PREPARE_GENOME='prepare_genome'
ALIGN_READS='align_reads'
RUN_TYPE=[PREPARE_GENOME,ALIGN_READS]
ANALYSIS_DIGITAL_TRANSC='digital_transc'
ANALYSIS_RNA_SEQ='rna_seq'
ANALYSIS_RAD_SEQ='rad'
ANALYSIS_TYPE=[ANALYSIS_DIGITAL_TRANSC,ANALYSIS_RNA_SEQ,ANALYSIS_RAD_SEQ]



def prepare_genome(genome_file,color_space=False):
    run_fine=True
    pipeline_param=utils_param.get_pipeline_parameters()
    BWA_dir=pipeline_param.get_bwa_dir()
    BWA_bin=os.path.join(BWA_dir,'bwa')
    genome_loader = GenomeLoader(genome_file=genome_file)
    length=0
    for fasta_rec in genome_loader:
        header, sequence = fasta_rec
        length+=len(sequence)
        if length>1000000000:
            break
    genome_loader.close()
    #Following recommendation set the indexing algorithm to is if genome is <10M
    if length>1000000000:
        a_option='bwtsw'
    else:
        a_option='is'
    
    #Create the indexes
    if color_space:
        command='%s index -c -a %s %s'%(BWA_bin, a_option, genome_file)
    else: 
        command='%s index -a %s %s'%(BWA_bin, a_option, genome_file)
    command_runner.run_command(command)
    return run_fine
    
    
def run_BWA_Command(genome_file, fastq_file1, fastq_file2=None, output_dir=None, sample_name=None,
                    clean_up=True, sort=False, thread=1, analysis_type=None, read_group=None, illumina=False):
    run_fine=True
    try:
        pipeline_param=utils_param.get_pipeline_parameters()
        BWA_dir=pipeline_param.get_bwa_dir()
        samtools_dir=pipeline_param.get_samtools_dir()
        picard_dir=pipeline_param.get_picard_dir()
    except Config_file_error, e:
        logging.exception('Config_file_error:')
        logging.warning("You'll need to have bwa and samtools in your path")
        BWA_dir=''
        samtools_dir=''
        picard_dir=None
        
    files_and_dir=[]
    if fastq_file1.endswith('.gz'):
        fastq_file1_unzip,ext = os.path.splitext(fastq_file1)
        command = 'gunzip -c %s > %s'%(fastq_file1,fastq_file1_unzip)
        return_code = command_runner.run_command(command)
        if return_code is not 0:
            run_fine = False
        fastq_file1=fastq_file1_unzip
        files_and_dir.append(fastq_file1)
    if fastq_file2 and fastq_file2.endswith('.gz'):
        fastq_file2_unzip,ext = os.path.splitext(fastq_file2)
        command = 'gunzip -c %s > %s'%(fastq_file2,fastq_file2_unzip)
        return_code = command_runner.run_command(command)
        if return_code is not 0:
            run_fine = False 
        fastq_file2=fastq_file2_unzip
        files_and_dir.append(fastq_file2)
    
    #Get the sample name
    if not sample_name:
        if fastq_file2:
            fastq_common = longest_common_substr_from_start(fastq_file1,fastq_file2)
            #remove trailing underscore _ and get the base name
            sample_name = os.path.basename(fastq_common.rstrip('_'))
        else:
            tmp,ext=os.path.splitext(os.path.basename(fastq_file1))
            tmp=tmp.rstrip('1')
            sample_name = os.path.basename(tmp.rstrip('_'))
    
    BWA_bin=os.path.join(BWA_dir,'bwa')
    # Check bwa version do not allow read group before version 0.9
    if compare_version_number(get_bwa_version(BWA_bin), "0.5.9") <0:
        logging.warning("Your version of bwa does not support the read group. Get version 0.5.9 or later to use this function.")
        read_group_command=''
    else:
        if read_group:
            read_group_command='-r "%s"'%(read_group)
        elif read_group is None:
            read_group_element=[]
            read_group_element.append("@RG")
            read_group_element.append("ID:%s"%sample_name)
            read_group_element.append("LB:%s"%sample_name)
            read_group_element.append("CN:The Genepool")
            read_group_element.append("PL:ILLUMINA")
            read_group_element.append("SM:%s"%sample_name)
            read_group_command= '-r "%s"'%('\\t'.join(read_group_element))
        else:
            read_group_command=''
    
    samtools_bin=os.path.join(samtools_dir,'samtools')
    
    if output_dir is None: 
        output_dir=os.path.dirname(fastq_file1)
    fastq_name, ext=os.path.splitext(os.path.basename(fastq_file1))
    sai_file1='%s.sai'%os.path.join(output_dir,fastq_name)
    illumina_str=""
    if illumina:
        illumina_str=" -I "
    command='%s aln %s -t %s %s %s > %s'%(BWA_bin, illumina_str,thread, genome_file, fastq_file1, sai_file1)
    if analysis_type == ANALYSIS_DIGITAL_TRANSC:
        #disable indels
        command='%s aln %s -o 0 -t %s %s %s > %s'%(BWA_bin, illumina_str, thread, genome_file, fastq_file1, sai_file1)
        
    return_code = command_runner.run_command(command)
    if return_code is not 0:
        run_fine = False 
    files_and_dir.append(sai_file1)
    
    if fastq_file2:
        fastq_name, ext=os.path.splitext(os.path.basename(fastq_file2))
        sai_file2='%s.sai'%os.path.join(output_dir,fastq_name)
        command='%s aln %s -t %s %s %s > %s'%(BWA_bin, illumina_str, thread, genome_file, fastq_file2, sai_file2)
        return_code = command_runner.run_command(command)
        files_and_dir.append(sai_file2)
        if return_code is not 0:
            run_fine = False
    
    bam_file=os.path.join(output_dir, sample_name+'.bam')
    
    if fastq_file2:
        if analysis_type == ANALYSIS_RNA_SEQ:
            command='%s sampe -a 20000 %s %s %s %s %s %s | %s view -bS - > %s'%(BWA_bin, read_group_command, genome_file, sai_file1, sai_file2,
                                                                                fastq_file1, fastq_file2, samtools_bin, bam_file)
        elif analysis_type == ANALYSIS_RAD_SEQ:
            command='%s sampe -A %s %s %s %s %s %s | %s view -bS - > %s'%(BWA_bin, read_group_command, genome_file, sai_file1, sai_file2, fastq_file1, 
                                                                          fastq_file2, samtools_bin,  bam_file)
        else:
            command='%s sampe %s %s %s %s %s %s | %s view -bS - > %s'%(BWA_bin, read_group_command, genome_file, sai_file1, sai_file2, fastq_file1, 
                                                                       fastq_file2, samtools_bin, bam_file)
    else:
        command='%s samse %s %s %s %s | %s view -bS - > %s'%(BWA_bin, read_group_command, genome_file, sai_file1, fastq_file1, samtools_bin, 
                                                             bam_file)
    return_code = command_runner.run_command( command)
    
    if return_code is not 0:
        run_fine = False
    
    
    if sort:
        files_and_dir.append(bam_file)
        if picard_dir:
            sorted_bam_file=os.path.join(output_dir,sample_name+'_sorted.bam')
            return_code = utils.sort_bam_file_per_coordinate(picard_dir, bam_file, sorted_bam_file, overwrite=True)
        else:
            sorted_bam_file=os.path.join(output_dir,sample_name+'_sorted')
            command='%s sort %s %s'%(samtools_bin, bam_file, sorted_bam_file)
            return_code = command_runner.run_command( command)
        if return_code is not 0:
            run_fine = False
    
            
    if run_fine and clean_up:
        return_code = remove_file(files_and_dir)
        if return_code is not 0:
            run_fine = False
    
    return run_fine
        
def remove_file(files_and_dir):
    return command_runner.run_command( 'rm -fr %s'%(' '.join(files_and_dir)))


def check_genome_index(genome_file):
    genome_index_valid=True
    extensions=['', '.amb', '.ann','.bwt','.sa', '.pac']
    for ext in extensions:
        if not os.path.exists(genome_file+ext):
            logging.error("%s doesn't exist"%(genome_file+ext))
            genome_index_valid=False
    return genome_index_valid


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
    if options.print_commands:
        utils_logging.change_log_stdout_to_log_stderr()
    else:
        command_runner.set_command_to_run_localy()
    run_fine=run_BWA_Command(options.genome_file, options.fastq_file1, options.fastq_file2,
                             options.output_dir, options.name,analysis_type=options.analysis,
                             sort=options.sort, thread=options.thread, read_group=options.read_group,
                             illumina=options.illumina)
    if run_fine:
        logging.info('Run completed')
    else:
        logging.error('Run Failed')
        sys.exit(1)
        
def _prepare_optparser():
    """Prepare optparser object. New options will be added in this
    function first.
    """
    usage = """usage: %prog <-g genome_fasta> <-1 first fastq file> [ -2 second fastq file -n sample_name]"""
    description = """This script will align read in Sanger fastq format to a reference genome and create a bam file. using bwa and samtools"""
    
    prog_version=utils.getWtss_version()
    optparser = OptionParser(version="%prog of wtss_pipeline v"+prog_version,description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="show this help message and exit.")
    optparser.add_option("-g","--genome_file",dest="genome_file",type="string",
                         help="Path to a fasta file where the genome is located. Default: %default")
    optparser.add_option("-1","--fastq1",dest="fastq_file1",type="string",
                         help="Path to the first fastq file where the first reads are. This file is mandatory. Default: %default")
    optparser.add_option("-2","--fastq2",dest="fastq_file2",type="string",
                         help="Path to the second fastq file where the second reads are. This file is optional. Default: %default")
    optparser.add_option("-o","--output_dir",dest="output_dir",type="string",
                         help="The path to the directory where the results will be output. If not set, the results will be put in the same folder as fastq1. Default: %default")
    optparser.add_option("-n","--name",dest="name",type="string",
                         help="The name of the sample currently being aligned. Default: %default")
    optparser.add_option("-s", "--sort",dest="sort",action='store_true',default=False,
                         help="Sort the bam file by coordinates at the end of the alignment. Default: %default")
    optparser.add_option("-t", "--thread",dest="thread",type='int',default=1,
                         help="Number of thread used by the alignment algorithm. Default: %default")
    optparser.add_option("--illumina",dest="illumina",action='store_true',default=False,
                         help="the fastq file are in illumina 1.3-1.6 fastq format. Default: %default")
    optparser.add_option("--print",dest="print_commands",action='store_true',default=False,
                         help="Print the command instead of running them. Default: %default")
    help_rna_seq="%s --> increase the maximum insert size to 20kb.\n"%(ANALYSIS_RNA_SEQ)
    help_digit_transc="%s --> prevent any gap in the tag.\n"%(ANALYSIS_DIGITAL_TRANSC)
    
    optparser.add_option("--analysis",dest="analysis",type="string",default=None,
                         help="Set analysis specific parameters:\n"+help_rna_seq+help_digit_transc+"Default: %default")
    optparser.add_option("-r", "--readgroup",dest="read_group",type="string",help="Set read group for SAM file:\n"+"Example: RG\tID:uid\tSM:sample\tPL:Illumina\n")
    
    return optparser


def _verifyOption(options):
    """Check if the mandatory option are present in the options objects.
    @return False if any argument is wrong."""
    arg_pass=True
    
    if not options.genome_file :
        logging.error("You must specify a genome fasta file.")
        arg_pass=False
    elif not os.path.exists(options.genome_file):
        logging.error("Genome fasta file not found. You must specify an existing genome fasta file.")
        arg_pass=False
    elif not check_genome_index(options.genome_file):
        logging.error("Genome fasta file is not indexed properly. You must index the genome first.")
        arg_pass=False
    #if not options.read_group:
    #    logging.error("You must specify a read group with ID, SM and PL attributes.")
    #    arg_pass=False
    if not options.fastq_file1:
        logging.error("You must specify at least one fastq file.")
        arg_pass=False
    elif not os.path.exists(options.fastq_file1):
        logging.error("fastq1 file %s not found. You must specify an existing file."%options.fastq_file1)
        arg_pass=False
    if options.fastq_file2:
        if not os.path.exists(options.fastq_file2):
            logging.error("fastq2 file %s not found. You must specify an existing file."%options.fastq_file2)
            arg_pass=False
    if options.output_dir:
        if not os.path.exists(options.output_dir):
            logging.error("output directory %s not found. You must specify an existing directory."%options.output_dir)
            arg_pass=False
    return arg_pass



if __name__=="__main__":
    main()
