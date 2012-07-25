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
from RAD_assemble_read2 import run_all_fastq_files


def run_smalt(consensus_file, read1_fastq, read2_fastq, **kwarg):    
    index1='%s.sma'%consensus_file
    index2='%s.smi'%consensus_file
    command='rm -rf %s'%index1
    command='rm -rf %s'%index2
    if os.path.exists(index1):
        return_code = command_runner.run_command(command)
    if os.path.exists(index2):
        return_code = command_runner.run_command(command)
    
    command = "smalt index %s %s"%(consensus_file, consensus_file)
    return_code = command_runner.run_command(command)
    name=longest_common_substr_from_start(read1_fastq, read2_fastq).rstrip('_')
    sam_file=name+'.sam'
    command = "smalt map -f samsoft -o %s %s %s %s"%(sam_file, consensus_file, read1_fastq, read2_fastq)
    return_code = command_runner.run_command(command)
    
    return sam_file;

def correct_smalt_sam_file(sam_file, readgroup_file):
    open_file = open(readgroup_file)
    all_read_groups=[]
    for line in open_file:
        all_read_groups.append(line.strip())
    open_file.close()
    
    first_read=True
    open_sam = open(sam_file)
    tmp,ext = os.path.splitext(sam_file)
    corrected_sam_file = tmp+"_corrected.sam"
    open_corrected_sam = open(corrected_sam_file, 'w')
    
    for line in open_sam:
        if line.startswith('@'):
            if line.startswith('@HD'):
                open_corrected_sam.write("@HD\tVN:1.3\tSO:unsorted\n")
            else:
                open_corrected_sam.write(line.strip()+'\n')
        else:
            if first_read:
                open_corrected_sam.write('\n'.join(all_read_groups)+'\n')
                first_read=False
            sp_line=line.strip().split()
            if sp_line[2] != "*":
                if int(sp_line[1]) & 8 == 8  and int(sp_line[1]) & 1024 == 1024:
                    #if my mate is unmapped and I'm a duplicate
                    #remove duplicate flag
                    sp_line[1]=str(int(sp_line[1])-1024)
                match=re.match('(.+)RGID:(.+)',sp_line[0])
                if match:
                    sp_line[0]=match.group(1)
                    sp_line.append("RG:Z:%s"%match.group(2))
                    open_corrected_sam.write('\t'.join(sp_line)+'\n')
    open_sam.close()
    open_corrected_sam.close()
    return corrected_sam_file
    



def run_alignment(consensus_file, read1_fastq, read2_fastq, readgroup_file):
    try:
        pipeline_param=utils_param.get_pipeline_parameters()
        picard_dir=pipeline_param.get_picard_dir()
        GATK_dir=pipeline_param.get_gatk_dir()
    except Config_file_error, e:
        logging.exception('Config_file_error:')
        sys.exit(1)

    #name=longest_common_substr_from_start(read1_fastq, read2_fastq).rstrip('_')
   
    #os.path.join(name+'_sorted_mrk_dup.bam')
 
    sam_file = run_smalt(consensus_file, read1_fastq, read2_fastq)
    corrected_sam_file = correct_smalt_sam_file(sam_file, readgroup_file)
    name, ext = os.path.splitext(corrected_sam_file)


    output_bam=os.path.join(name+"_sorted.bam")
    sort_bam_file_per_coordinate(picard_dir, input_bam=corrected_sam_file, output_bam=output_bam, overwrite=True,
                                 CREATE_INDEX="true")
    
    mark_dups_jar = os.path.join(picard_dir, 'MarkDuplicates.jar')
    mark_dups_bam = os.path.join(name+'_sorted_mrk_dup.bam')
    mark_dups_metric = os.path.join(name+'_sorted_mrk_dup.metric')
    command = 'java -Xmx5G -jar %s I=%s O=%s METRICS_FILE=%s  VALIDATION_STRINGENCY=LENIENT MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=100 CREATE_INDEX=true'%(mark_dups_jar, output_bam,
                                                                                                                                     mark_dups_bam, mark_dups_metric)
    command_runner.run_command(command)
    fixed_bam=os.path.join(name+'_sorted_mrk_dup_fixed.bam')
    command="""samtools view -h %s | 
awk '{if (and($2,8)==8 && and($2,1024)==1024){$2=$2-1024}; if (!/^@/){$2+=2} print $0; }' |
sed 's/ /\\t/g' | sed 's/The\\tGenepool/The Genepool/' | samtools view -Sb - > %s""".replace("\n"," ")
    command=command%(mark_dups_bam,fixed_bam)
    command_runner.run_command(command)
    
    command="samtools index %s"%fixed_bam
    command_runner.run_command(command)
    
    GATK_jar=os.path.join(GATK_dir,"GenomeAnalysisTK.jar")
    gatk_all_sites_vcf=os.path.join(name+'_sorted_mrk_dup_fixed_GATK_all_sites.vcf')
    command="java -Xmx1G -jar %s -T UnifiedGenotyper -nt 1 -I %s -R %s -o %s -out_mode EMIT_ALL_CONFIDENT_SITES --downsampling_type NONE -glm BOTH"
    command=command%(GATK_jar, fixed_bam, consensus_file, gatk_all_sites_vcf)    
    command_runner.run_command(command)
    
    gatk_var_sites_vcf=os.path.join(name+'_sorted_mrk_dup_fixed_GATK_var_sites.vcf')
    command="""awk '{if (/^#/ || $5!="."){print}}' %s > %s"""%(gatk_all_sites_vcf, gatk_var_sites_vcf)
    command_runner.run_command(command)
    
    gatk_var_sites_filter_vcf=os.path.join(name+'_sorted_mrk_dup_fixed_GATK_var_sites_filterd20q60.vcf')
    command = "vcfutils.pl varFilter -d 20 %s | awk '{if (/^#/ || $6>60){print}}' > %s"%(gatk_var_sites_vcf, gatk_var_sites_filter_vcf)
    command_runner.run_command(command)

    gatk_var_sites_filter_phased_vcf=os.path.join(name+'_sorted_mrk_dup_fixed_GATK_var_sites_filterd20q60_phased.vcf')
    command="java -Xmx1G -jar %s -T ReadBackedPhasing -R %s -I %s --variant %s -L %s -o %s --phaseQualityThresh 20.0"
    command=command%(GATK_jar, consensus_file, fixed_bam, gatk_var_sites_filter_vcf, gatk_var_sites_filter_vcf, gatk_var_sites_filter_phased_vcf)
    command_runner.run_command(command)

    logging.info("\n")


def run_all_fastq_files(directory,readgroup_file):
    directory=os.path.abspath(directory)
    all_dirs = glob(os.path.join(directory,'consensus*_dir'))

    for dir in all_dirs:
        name=os.path.basename(dir)[:-len("_dir")]
        
        consensus_file=os.path.join(dir,'best_assembly.fa')
        read1_fastq=os.path.join(dir,name+"_1.fastq")
        read2_fastq=os.path.join(dir,name+"_2.fastq")
        run_alignment(consensus_file, read1_fastq, read2_fastq, readgroup_file)
        

    
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
    start_time = time.time()
    run_all_fastq_files(options.consensus_dir, options.readgroup_file)
    logging.info("Elapsed time:%.1f seconds"%(time.time()-start_time))

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
    optparser.add_option("-r","--readgroup",dest="readgroup_file",type="string",default=None,
                         help="The name of the assembler that will be used on the fastq files. Default: %default")
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
    
