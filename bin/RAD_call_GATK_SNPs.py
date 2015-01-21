#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 19 March  2012
@author: tcezard
'''
import sys
import os
import logging
from optparse import OptionParser
from glob import glob

from RAD_merge_results import concatenate_file
from utils import utils_logging, utils_param
import command_runner
from utils.parameters import Config_file_error


def SNP_call_with_GATK(GATK_dir, Picard_dir, samtools_bin, name, bam_file, ref_file):
    CreateSeqDict_jar = os.path.join(Picard_dir, "CreateSequenceDictionary.jar")
    base, ext = os.path.splitext(ref_file)
    dict_file = base + '.dict'

    command = "java -jar %s R=%s O=%s" % (CreateSeqDict_jar, ref_file, dict_file)
    command_runner.run_command(command)

    command = "%s faidx %s" % (samtools_bin, ref_file)
    command_runner.run_command(command)
    base_dir = os.path.dirname(bam_file)
    GATK_jar = os.path.join(GATK_dir, "GenomeAnalysisTK.jar")
    gatk_all_sites_vcf = os.path.join(base_dir, name + '_sorted_mrk_dup_fixed_GATK_all_sites.vcf')
    command = "java -Xmx1G -jar %s -T UnifiedGenotyper -nt 1 -I %s -R %s -o %s -out_mode EMIT_ALL_CONFIDENT_SITES --downsampling_type NONE -glm BOTH"
    command = command % (GATK_jar, bam_file, ref_file, gatk_all_sites_vcf)
    command_runner.run_command(command)

    gatk_var_sites_vcf = os.path.join(base_dir, name + '_sorted_mrk_dup_fixed_GATK_var_sites.vcf')
    command = """awk '{if (/^#/ || $5!="."){print}}' %s > %s""" % (gatk_all_sites_vcf, gatk_var_sites_vcf)
    command_runner.run_command(command)

    gatk_var_sites_filter_vcf = os.path.join(base_dir, name + '_sorted_mrk_dup_fixed_GATK_var_sites_filterd20q60.vcf')
    command = "vcfutils.pl varFilter -d 20 %s | awk '{if (/^#/ || $6>60){print}}' > %s" % (
    gatk_var_sites_vcf, gatk_var_sites_filter_vcf)
    command_runner.run_command(command)

    gatk_var_sites_filter_phased_vcf = os.path.join(base_dir,
                                                    name + '_sorted_mrk_dup_fixed_GATK_var_sites_filterd20q60_phased.vcf')
    command = "java -Xmx1G -jar %s -T ReadBackedPhasing -R %s -I %s --variant %s -L %s -o %s --phaseQualityThresh 20.0"
    command = command % (GATK_jar, ref_file, bam_file, gatk_var_sites_filter_vcf, gatk_var_sites_filter_vcf,
                         gatk_var_sites_filter_phased_vcf)
    command_runner.run_command(command)


def merge_snps_files(directory):
    """This function merge snps files from a single directory"""
    return_code = 0
    all_vcf_files = glob(os.path.join(directory, '*_dir', '*_phased.vcf'))
    output_file_body = os.path.join(directory, '%s_snps_files.vcf.body' % len(all_vcf_files))
    output_file_body = concatenate_file(all_vcf_files, output_file_body, filter="^#")
    if output_file_body:
        return_code = 0
    output_file_header = os.path.join(directory, '%s_snps_files.vcf.header' % len(all_vcf_files))
    command = 'grep  "^#" %s > %s ' % (all_vcf_files[0], output_file_header)
    if return_code == 0:
        return_code = command_runner.run_command(command)
    output_file = os.path.join(directory, '%s_phased_snps_files.vcf' % len(all_vcf_files))
    command = 'cat %s %s > %s ' % (output_file_header, output_file_body, output_file)
    if return_code == 0:
        return_code = command_runner.run_command(command)
    command = 'rm %s %s' % (output_file_header, output_file_body)
    if return_code == 0:
        return_code = command_runner.run_command(command)
    return return_code


def run_all_snps_call(directory):
    directory = os.path.abspath(directory)
    all_dirs = glob(os.path.join(directory, '*_dir'))
    GATK_dir = None
    try:
        pipeline_param = utils_param.get_pipeline_parameters()
        picard_dir = pipeline_param.get_picard_dir()
        GATK_dir = pipeline_param.get_gatk_dir()
        samtools_dir = pipeline_param.get_samtools_dir()
        samtools_bin = os.path.join(samtools_dir, 'samtools')
    except Config_file_error, e:
        logging.exception('Config_file_error:')
        sys.exit(1)

    for sub_dir in all_dirs:
        name = os.path.basename(sub_dir)[:-len("_dir")]
        consensus_file = os.path.join(sub_dir, 'best_assembly.fa')
        read1_fastq = os.path.join(sub_dir, name + "_1.fastq")
        read2_fastq = os.path.join(sub_dir, name + "_2.fastq")
        bam_file = os.path.join(sub_dir, name + "_corrected_sorted_mrk_dup_fixed.bam")
        SNP_call_with_GATK(GATK_dir, picard_dir, samtools_bin, name, bam_file, ref_file=consensus_file)
    merge_snps_files(directory)


def main():
    # initialize the logging
    utils_logging.init_logging(logging.INFO)
    #Setup options
    optparser = _prepare_optparser()
    (options, args) = optparser.parse_args()
    #verify options
    arg_pass = _verifyOption(options)
    if not arg_pass:
        logging.warning(optparser.get_usage())
        logging.critical("Non valid arguments: exit")
        sys.exit(1)
    command_runner.set_command_to_run_localy()
    if options.debug:
        utils_logging.init_logging(logging.DEBUG)
    run_all_snps_call(options.consensus_dir)


def _prepare_optparser():
    """Prepare optparser object. New options will be added in this
    function first.
    """
    usage = """usage: %prog <-b bam_file> [ -o output_file]"""
    description = """This script will take RAD reads aligned to there read1+read2 consensus and call SNP with GATK."""

    optparser = OptionParser(version="None", description=description, usage=usage, add_help_option=False)
    optparser.add_option("-h", "--help", action="help", help="show this help message and exit.")
    optparser.add_option("-d", "--consensus_dir", dest="consensus_dir", type="string",
                         help="Path to a directory containing fastq file (only extension .fastq will be processed). Default: %default")
    optparser.add_option("--debug", dest="debug", action='store_true', default=False,
                         help="Output debug statment. Default: %default")
    return optparser


def _verifyOption(options):
    """Check if the mandatory option are present in the options objects.
    @return False if any argument is wrong."""
    arg_pass = True

    return arg_pass


if __name__ == "__main__":
    main()
    
    