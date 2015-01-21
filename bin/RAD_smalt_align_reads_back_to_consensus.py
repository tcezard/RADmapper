#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 19 March  2012
@author: tcezard
'''
import sys
import os
import logging
import re
from optparse import OptionParser
from glob import glob
import time

from utils import utils_logging, longest_common_substr_from_start, utils_param, sort_bam_file_per_coordinate
import command_runner
from utils.parameters import Config_file_error
from utils import DNA_tools


def clean_up_prev_run(name):
    to_rm = ['.sam',
             '_corrected.sam',
             '_corrected_sorted.bai',
             '_corrected_sorted.bai',
             '_corrected_sorted.bam',
             '_corrected_sorted_mrk_dup.bai',
             '_corrected_sorted_mrk_dup.bam',
             '_corrected_sorted_mrk_dup_fixed.bam',
             '_corrected_sorted_mrk_dup_fixed.bam.bai',
             '_corrected_sorted_mrk_dup_fixed.bam.stat',
             '_corrected_sorted_mrk_dup_fixed_GATK_all_sites.vcf',
             '_corrected_sorted_mrk_dup_fixed_GATK_all_sites.vcf.idx',
             '_corrected_sorted_mrk_dup_fixed_GATK_var_sites_filterd20q60.vcf',
             '_corrected_sorted_mrk_dup_fixed_GATK_var_sites_filterd20q60.vcf.idx',
             '_corrected_sorted_mrk_dup_fixed_GATK_var_sites.vcf',
             '_corrected_sorted_mrk_dup_fixed_GATK_var_sites.vcf.idx',
             '_corrected_sorted_mrk_dup_fixed_GATK_var_sites_filterd20q60_phased.vcf',
             '_corrected_sorted_mrk_dup_fixed_GATK_var_sites_filterd20q60_phased.vcf.idx',
             '_corrected_sorted_mrk_dup_fixed_samtools.vcf',
             '_corrected_sorted_mrk_dup_fixed_samtools_filterd20q60.vcf',
             '_corrected_sorted_mrk_dup.metric',
             '.sam',
             '_summary_stat.txt']
    logging.info('Clean up previous Run')
    for val in to_rm:
        file_name = name + val
        if os.path.exists(file_name):
            logging.debug('remove %s' % file_name)
            os.remove(file_name)
        else:
            logging.debug('%s does not exist' % file_name)


def reverse_complement(fastq_file):
    rev_comp_file_name = fastq_file + '.rev_comp'
    with open(rev_comp_file_name, 'w') as open_rev_com:
        line_number = 0
        with open(fastq_file) as open_fastq:
            for line in open_fastq:
                line_number += 1
                if line_number % 4 == 2:
                    line = DNA_tools.rev_complements(line.strip()) + '\n'
                elif line_number % 4 == 0:
                    line = line.strip()[::-1] + '\n'
                open_rev_com.write(line)
    return rev_comp_file_name


def run_smalt(consensus_file, read1_fastq, read2_fastq, **kwarg):
    index1 = '%s.sma' % consensus_file
    command = 'rm -rf %s' % index1
    if os.path.exists(index1):
        return_code = command_runner.run_command(command)
    index2 = '%s.smi' % consensus_file
    command = 'rm -rf %s' % index2
    if os.path.exists(index2):
        return_code = command_runner.run_command(command)
    index3 = '%s.fai' % consensus_file
    command = 'rm -rf %s' % index3
    if os.path.exists(index3):
        return_code = command_runner.run_command(command)

    command = "smalt index %s %s" % (consensus_file, consensus_file)
    return_code = command_runner.run_command(command)
    name = longest_common_substr_from_start(read1_fastq, read2_fastq).rstrip('_')
    read2_fastq_rev_comp = reverse_complement(read2_fastq)

    sam_file = name + '.sam'
    command = "smalt map -f samsoft -o %s %s %s %s" % (sam_file, consensus_file, read1_fastq, read2_fastq_rev_comp)
    return_code = command_runner.run_command(command)

    return sam_file;


def read_readgroup_file(readgroup_file):
    open_file = open(readgroup_file)
    all_read_groups = []
    for line in open_file:
        all_read_groups.append(line.strip())
    open_file.close()
    return all_read_groups


def correct_smalt_sam_file(sam_file, all_read_groups):
    first_read = True
    open_sam = open(sam_file)
    tmp, ext = os.path.splitext(sam_file)
    corrected_sam_file = tmp + "_corrected.sam"
    open_corrected_sam = open(corrected_sam_file, 'w')

    for line in open_sam:
        if line.startswith('@'):
            if line.startswith('@HD'):
                open_corrected_sam.write("@HD\tVN:1.3\tSO:unsorted\n")
            else:
                open_corrected_sam.write(line.strip() + '\n')
        else:
            if first_read:
                open_corrected_sam.write('\n'.join(all_read_groups) + '\n')
                first_read = False
            sp_line = line.strip().split()
            match = re.match('(.+)RGID:(.+)', sp_line[0])
            if match:
                sp_line[0] = match.group(1)
                sp_line.append("RG:Z:%s" % match.group(2))
                open_corrected_sam.write('\t'.join(sp_line) + '\n')
    open_sam.close()
    open_corrected_sam.close()
    return corrected_sam_file


def run_alignment(consensus_file, read1_fastq, read2_fastq, all_read_groups, snp_call=False):
    try:
        pipeline_param = utils_param.get_pipeline_parameters()
        picard_dir = pipeline_param.get_picard_dir()
        # GATK_dir=pipeline_param.get_gatk_dir()
        samtools_dir = pipeline_param.get_samtools_dir()
    except Config_file_error, e:
        logging.exception('Config_file_error:')
        sys.exit(1)

    # Cleanup previous run if it exist
    name = longest_common_substr_from_start(read1_fastq, read2_fastq).rstrip('_')
    clean_up_prev_run(name)

    file_to_remove = []
    sam_file = run_smalt(consensus_file, read1_fastq, read2_fastq)
    file_to_remove.append(sam_file)
    corrected_sam_file = correct_smalt_sam_file(sam_file, all_read_groups)
    file_to_remove.append(corrected_sam_file)
    name, ext = os.path.splitext(corrected_sam_file)

    output_bam = os.path.join(name + "_sorted.bam")
    sort_bam_file_per_coordinate(picard_dir, input_bam=corrected_sam_file, output_bam=output_bam, overwrite=True,
                                 CREATE_INDEX="true")

    file_to_remove.append(output_bam)

    mark_dups_jar = os.path.join(picard_dir, 'MarkDuplicates.jar')
    mark_dups_bam = os.path.join(name + '_sorted_mrk_dup.bam')
    mark_dups_metric = os.path.join(name + '_sorted_mrk_dup.metric')
    command = 'java -Xmx5G -jar %s I=%s O=%s METRICS_FILE=%s  VALIDATION_STRINGENCY=LENIENT MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=100 CREATE_INDEX=true' % (
    mark_dups_jar, output_bam,
    mark_dups_bam, mark_dups_metric)
    command_runner.run_command(command)
    file_to_remove.append(mark_dups_bam)
    fixed_bam = os.path.join(name + '_sorted_mrk_dup_fixed.bam')
    #This command remove the duplicate flag when a read is mapped and its mate isn't
    #It also remove the unmapped read from the bam file as this prevent the merging for some reason !!
    command = """samtools view -h %s |
awk 'BEGIN{OFS="\t"}{if (and($2,8)==8 && and($2,1024)==1024){$2=$2-1024}; if(and($2,4)!=4){print $0} }' |
samtools view -Sb - > %s""".replace("\n", " ")
    command = command % (mark_dups_bam, fixed_bam)
    command_runner.run_command(command)

    command = "samtools index %s" % fixed_bam
    command_runner.run_command(command)

    if snp_call:
        SNP_call_with_samtools(samtools_dir, name, bam_file=fixed_bam, ref_file=consensus_file)
        logging.info("\n")


def SNP_call_with_samtools(samtools_dir, name, bam_file, ref_file):
    samtools_bin = os.path.join(samtools_dir, "samtools")
    bcftools_bin = os.path.join(samtools_dir, "bcftools/bcftools")
    if not os.path.exists(bcftools_bin):
        bcftools_bin = os.path.join(samtools_dir, "bcftools")
    samtools_raw_vcf = os.path.join(name + '_sorted_mrk_dup_fixed_samtools.vcf')
    command = "%s mpileup -d 2000 -ADESuf %s %s | %s view -gv - > %s"
    command = command % (samtools_bin, ref_file, bam_file, bcftools_bin, samtools_raw_vcf)
    command_runner.run_command(command)

    samtools_raw_filtered = os.path.join(name + '_sorted_mrk_dup_fixed_samtools_filterd20q60.vcf')
    command = "vcfutils.pl varFilter -d 20 %s | awk '{if (/^#/ || $6>60){print}}' > %s" % (
    samtools_raw_vcf, samtools_raw_filtered)
    command_runner.run_command(command)


def run_all_fastq_files(directory, readgroup_file=None, all_read_groups=None, snp_call=False):
    """This function will align all the reads back to the best assembly for each directory"""
    if readgroup_file:
        all_read_groups = read_readgroup_file(readgroup_file)
    directory = os.path.abspath(directory)
    all_dirs = glob(os.path.join(directory, '*_dir'))

    for sub_dir in all_dirs:
        name = os.path.basename(sub_dir)[:-len("_dir")]
        consensus_file = os.path.join(sub_dir, 'best_assembly.fa')
        read1_fastq = os.path.join(sub_dir, name + "_1.fastq")
        read2_fastq = os.path.join(sub_dir, name + "_2.fastq")
        run_alignment(consensus_file, read1_fastq, read2_fastq, all_read_groups, snp_call)


def main():
    # initialize the logging
    utils_logging.init_logging(logging.INFO)
    #utils_logging.init_logging(logging.CRITICAL)
    #Setup options
    optparser = _prepare_optparser()
    (options, args) = optparser.parse_args()
    #verify options
    arg_pass = _verifyOption(options)
    if not arg_pass:
        logging.warning(optparser.get_usage())
        logging.critical("Non valid arguments: exit")
        sys.exit(1)
    if options.debug:
        utils_logging.init_logging(logging.DEBUG)
    if not options.print_command:
        command_runner.set_command_to_run_localy()
    start_time = time.time()
    run_all_fastq_files(options.consensus_dir, readgroup_file=options.readgroup_file, snp_call=options.snp_call)
    logging.info("Elapsed time:%.1f seconds" % (time.time() - start_time))


def _prepare_optparser():
    """Prepare optparser object. New options will be added in this
    function first.
    """
    usage = """usage: %prog <-b bam_file> [ -o output_file]"""
    description = """This script will take aligned RAD read to the consensuses and calculate per consensus coverage."""

    optparser = OptionParser(version="None", description=description, usage=usage, add_help_option=False)
    optparser.add_option("-h", "--help", action="help", help="show this help message and exit.")
    optparser.add_option("-d", "--consensus_dir", dest="consensus_dir", type="string",
                         help="Path to a directory containing fastq file (only extension .fastq will be processed). Default: %default")
    optparser.add_option("-r", "--readgroup", dest="readgroup_file", type="string", default=None,
                         help="The name of the assembler that will be used on the fastq files. Default: %default")
    optparser.add_option("-s", "--snp_call", dest="snp_call", action='store_true', default=False,
                         help="Run the SNPs call with GATK at the end of the alignment. Default: %default")
    optparser.add_option("--print", dest="print_command", action='store_true', default=False,
                         help="print the commands instead of running them. Default: %default")
    optparser.add_option("--debug", dest="debug", action='store_true', default=False,
                         help="Output debug statment. Default: %default")
    return optparser


def _verifyOption(options):
    """Check if the mandatory option are present in the options objects.
    @return False if any argument is wrong."""
    arg_pass = True

    if not options.consensus_dir or not os.path.isdir(options.consensus_dir):
        logging.error("You need to specify a valid directory with -d")
        arg_pass = False
    if not options.readgroup_file or not os.path.isfile(options.readgroup_file):
        logging.error("You need to specify a valid readgroup file with -r")
        arg_pass = False
    return arg_pass


if __name__ == "__main__":
    main()
