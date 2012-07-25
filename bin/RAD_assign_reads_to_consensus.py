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



def run_BWA_Command_for_RAD(genome_file, fastq_file1, fastq_file2=None, output_dir=None, sample_name=None,
                    clean_up=True, thread=1, read_group=None, illumina=False, fifo=False):
    run_fine=True
    try:
        pipeline_param=utils_param.get_pipeline_parameters()
        BWA_dir=pipeline_param.get_bwa_dir()
        samtools_dir=pipeline_param.get_samtools_dir()
        picard_dir=pipeline_param.get_picard_dir()
        #get the path to the current script to infer the path to RAD_set_read1_consensus_to_read2.py
        script_dir = os.path.dirname(os.path.realpath(__file__))
        path_to_script = os.path.join(script_dir,"RAD_set_read1_consensus_to_read2.py")
        if not os.path.exists(path_to_script):
            logging.warning("Can't find RAD_set_read1_consensus_to_read2.py which is required")
            path_to_script = "RAD_set_read1_consensus_to_read2.py"
            
    except Config_file_error, e:
        #logging.exception('Config_file_error:')
        logging.critical("You need to have the environment variable properly set to use that script")
        return False
        
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
            tmp, ext=os.path.splitext(os.path.basename(fastq_file1))
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
            for element in read_group.split('\\t'):
                if element.startswith("ID"):
                    rgid=element.split(":")[1]
                elif element.startswith("LB"):
                    libid=element.split(":")[1]
                elif element.startswith("SM"):
                    smid=element.split(":")[1]
                
        elif read_group is None:
            rgid=sample_name
            libid=sample_name
            smid=sample_name
            read_group_element=[]
            read_group_element.append("@RG")
            read_group_element.append("ID:%s"%sample_name)
            read_group_element.append("LB:%s"%sample_name)
            read_group_element.append("CN:The Genepool")
            read_group_element.append("PL:ILLUMINA")
            read_group_element.append("SM:%s"%sample_name)
            read_group_command= '-r "%s"'%('\\t'.join(read_group_element))
            
    
    samtools_bin=os.path.join(samtools_dir,'samtools')
    
    if output_dir is None: 
        output_dir=os.path.dirname(fastq_file1)
    fastq_name, ext=os.path.splitext(os.path.basename(fastq_file1))
    sai_file1='%s.sai'%os.path.join(output_dir,fastq_name)
    illumina_str=""
    if illumina:
        illumina_str=" -I "
    command='%s aln %s -t %s %s %s > %s'%(BWA_bin, illumina_str,thread, genome_file, fastq_file1, sai_file1)
        
    return_code = command_runner.run_command(command)
    if return_code is not 0:
        run_fine = False 
    files_and_dir.append(sai_file1)
    if not fastq_file2:
        #only one end so just run get the sorted bam file
        bam_file=os.path.join(output_dir, sample_name)
        command="""%s samse %s %s %s %s | %s view -bS - > %s"""%(BWA_bin, read_group_command, genome_file, sai_file1, fastq_file1, samtools_bin, bam_file )
        return_code = command_runner.run_command( command)
        if return_code is not 0:
            run_fine = False 
        files_and_dir.append(bam_file)
        sorted_bam_file=os.path.join(output_dir, sample_name+'_sorted')
        return_code = utils.sort_bam_file_per_coordinate(picard_dir, bam_file, sorted_bam_file, overwrite=True)
        if return_code is not 0:
            run_fine = False 
    else:
        #Need to add to the first read the flag info specific to the paired end
        sam_file=os.path.join(output_dir, sample_name+'.sam')
        if fifo:
            command="mkfifo %s"%sam_file
            if os.path.exists(sam_file):
                os.remove(sam_file)
            return_code = command_runner.run_command( command)
        command="""%s samse %s %s %s %s """%(BWA_bin, read_group_command,genome_file,sai_file1,fastq_file1)
        command+="""| awk '{if (/^@/){ print }else{ $2+=73; printf "%s",$1; for (i=2;i<=NF;i++){printf "\\t%s",$i}; printf "\\n" }}'"""
        command+=""" > %s"""%( sam_file )
        if fifo:
            command+=" &"
        return_code = command_runner.run_command( command)
        if return_code is not 0:
            run_fine = False 
        files_and_dir.append(sam_file)
    
        #Now convert the second read to sam
        fastq_name, ext=os.path.splitext(os.path.basename(fastq_file2))
        sam_file2='%s.sam'%os.path.join(output_dir,fastq_name)
        if fifo:
            command="mkfifo %s"%sam_file2
            if os.path.exists(sam_file2):
                os.remove(sam_file2)
            return_code = command_runner.run_command( command)
        
        fastqToSam_jar=os.path.join(picard_dir,"FastqToSam.jar")
        if illumina:
            qual="Illumina"
        else:
            qual="Standard"
        command="java -Xmx2G -jar %s F1=%s O=%s RG=%s LB=%s SM=%s, QUALITY_FORMAT=%s"%(fastqToSam_jar,fastq_file2,sam_file2,rgid,libid,smid,qual)
        if fifo:
            command+=" &"
        return_code = command_runner.run_command(command)
        if return_code is not 0:
            run_fine = False 
        files_and_dir.append(sam_file2)
        sam_file2_tmp='%s.sam.tmp'%os.path.join(output_dir,fastq_name)
        if fifo:
            if os.path.exists(sam_file2_tmp):
                os.remove(sam_file2_tmp)
            command="mkfifo %s"%sam_file2_tmp
            return_code = command_runner.run_command( command)
        command = """awk '{if (/^@/){}else{if (substr($1,length($1)-1,2)=="/2"){$1=substr($1,1,length($1)-2)}; $2+=129; printf "%s",$1; for (i=2;i<=NF;i++){printf "\\t%s",$i}; printf "\\n"}}' """
        command+=""" %s > %s"""%(sam_file2,sam_file2_tmp)
        if fifo:
            command+=" &"
        return_code = command_runner.run_command(command)
        if return_code is not 0:
            run_fine = False 
        files_and_dir.append(sam_file2_tmp)
        
        name_sorted_bam_file=os.path.join(output_dir, sample_name+'_sorted_name')
        command = "cat %s %s | %s view -bS - | %s sort -n - %s"%(sam_file, sam_file2_tmp, samtools_bin, samtools_bin, name_sorted_bam_file)
        return_code = command_runner.run_command(command)
        if return_code is not 0:
            run_fine = False
        files_and_dir.append(name_sorted_bam_file+'.bam')
        
        sorted_bam_file=os.path.join(output_dir, sample_name+'_sorted')
        #change_consensus_on_read2
        command ="%s view -h %s | "%(samtools_bin,name_sorted_bam_file+'.bam')
        command +="python %s | "%path_to_script
        command +="%s view -bS - | %s sort - %s"%(samtools_bin,  samtools_bin, sorted_bam_file)
        return_code = command_runner.run_command(command)
        if return_code is not 0:
            run_fine = False
    
        sorted_bam_file=sorted_bam_file+".bam"
        command="%s index %s"%(samtools_bin, sorted_bam_file)
        return_code = command_runner.run_command(command)
        
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
    #extensions=['','.rbwt', '.amb', '.rpac', '.ann', '.rsa','.bwt','.sa', '.pac']
    extensions=['','.amb', '.ann', '.bwt','.sa', '.pac']
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
    run_fine=run_BWA_Command_for_RAD(options.genome_file, options.fastq_file1, options.fastq_file2,
                             options.output_dir, options.name, thread=options.thread, read_group=options.read_group,
                             illumina=options.illumina, fifo=options.fifo)
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
    optparser.add_option("-t", "--thread",dest="thread",type='int',default=1,
                         help="Number of thread used by the alignment algorithm. Default: %default")
    optparser.add_option("--illumina",dest="illumina",action='store_true',default=False,
                         help="the fastq file are in illumina 1.3-1.6 fastq format. Default: %default")
    optparser.add_option("--print",dest="print_commands",action='store_true',default=False,
                         help="Print the command instead of running them. Default: %default")
    optparser.add_option("--fifo",dest="fifo",action='store_true',default=False,
                         help="Use fifo to avoid writing big temp file to the disk. Default: %default")
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
