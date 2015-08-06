#!/usr/bin/env python
'''
Created on Mar 9, 2011

@author: tcezard
'''
import sys
import logging
from optparse import OptionParser
import os
from collections import Counter, defaultdict
import getpass
import time
import re

from utils import utils_logging
from IO_interface import vcfIO
from utils import atgc2iupac


def read_pop_file(pop_file):
    sample2pop = {}
    pop2sample = defaultdict(list)
    with open(pop_file) as open_file:
        for line in open_file:
            sp_line = line.strip().split()
            if len(sp_line) > 1:
                pop = sp_line[1]
                sample = sp_line[0]
                sample2pop[sample] = pop
                pop2sample[pop].append(sample)
    return sample2pop, pop2sample



"""1st column Id for each family
2nd column individual Id
3rd Sire id
4th Dam id
5th sex (doesn't play a role just ensure your parents are coded as male
and female) 1 for male 2 for female
6th just put 0
7th onwards markers (both alleles are coded separated by space)"""

def vcf_to_lepmap(vcf_file, output_file, sex_file, mother_name, father_name, family_name,
                  genotype_quality_threshold=20, max_prop_missing=.5):
    file_handle = utils_logging.open_input_file(vcf_file, pipe=False)
    reader = vcfIO.VcfReader(file_handle)
    all_samples_in_file = reader.get_sample_names()
    if sex_file:
        sample2sex, sex2sample = read_pop_file(sex_file)
        all_samples_in_file = sample2sex.keys()
    if mother_name in all_samples_in_file: all_samples_in_file.pop(mother_name)
    if father_name in all_samples_in_file: all_samples_in_file.pop(father_name)

    max_missing = int(len(all_samples_in_file) * max_prop_missing)
    valid_genotypes = ['0/1 0/0', 
                      '0/1 0/1',
                      '0/1 1/1',
                      '0/0 0/1',
                      '0/0 1/1',
                      '1/1 0/0',
                      '1/1 0/1']
    all_lines = {}
    for sample in all_samples_in_file:
        if sample2sex.get(sample) == 'M': sex="1"
        else: sex="2"
        all_lines[sample] = [family_name, sample, father_name, mother_name, sex, "0"]
    all_lines[father_name] = [family_name, father_name, "0", "0", "1", "0"]
    all_lines[mother_name] = [family_name, mother_name, "0", "0", "2", "0"]
    all_samples_in_file.append(father_name)
    all_samples_in_file.append(mother_name)
    nb_lines = 0
    nb_sequence = 0
    count_more_than_one_allele = 0
    count_indel = 0
    count_too_many_missing = 0
    count_non_polymorphic = 0
    nb_missing_per_sample = Counter()
    count_invalid_parent_geno = 0
    for vcf_records in reader:
        nb_lines += 1
        if nb_lines % 10000 == 0:
            sys.stdout.write('.')
        ref_base = vcf_records.get_reference_base()
        alt_bases = vcf_records.get_alt_bases()
        if len(alt_bases) > 1:
            count_more_than_one_allele += 1
            continue
        if vcf_records.is_indel():
            count_indel += 1
            continue

        nb_missing = 0
        all_chars = []
        all_codes = set()
        parent_geno = []
        for sample in all_samples_in_file:
            gt = vcf_records.get_genotype(sample)
            gq = vcf_records.get_genotype_quality(sample)
            if gt and gq > genotype_quality_threshold:
                value1, value2 = re.split('[/|]', gt)
                code = '%s %s'%(int(value1)+1, int(value2)+1)
            else:
                nb_missing += 1
                nb_missing_per_sample[sample] += 1
                code = '0 0'
            all_chars.append(code)
            all_codes.add(code)
            if sample in [father_name, mother_name] and gq>genotype_quality_threshold:
                parent_geno.append(gt)
        #if len(parent_geno)!=2 or not ' '.join(parent_geno) in valid_genotypes:
        #    count_invalid_parent_geno += 1
        #    continue
        if len(all_codes) == 1:
            count_non_polymorphic += 1
            continue
        if nb_missing <= max_missing:
            nb_sequence += 1
            for i, sample in enumerate(all_samples_in_file):
                all_lines[sample].append(all_chars[i])
        else:
            count_too_many_missing += 1
    if count_invalid_parent_geno:
        logging.warning("%s snps remove because they missing or non informative parental genotypes" % (count_invalid_parent_geno))
    if count_more_than_one_allele:
        logging.warning("%s snps remove because they had more than 2 alleles" % (count_more_than_one_allele))
    if count_indel:
        logging.warning("%s indels removed" % (count_indel))
    if count_non_polymorphic:
        logging.warning(
            "%s snps removed because no polymorphism was found between populations" % (count_non_polymorphic))
    if count_too_many_missing:
        logging.warning("%s snps removed because >%s missing samples" % (count_too_many_missing, max_missing))
    logging.info("%s snps output in Lepmap format" % (nb_sequence))
    for sample in nb_missing_per_sample:
        logging.info("%s markers missing in %s" % (nb_missing_per_sample.get(sample), sample))
    with open(output_file, 'w') as open_output:
        for sample in all_lines:
            open_output.write('%s\n'%('\t'.join(all_lines.get(sample))))


def main():
    # initialize the logging
    utils_logging.init_logging()
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

    vcf_to_lepmap(options.input_vcf_file, options.output_file, options.sex_file, options.mother_name,
                  options.father_name, options.family_name, options.geno_qual_threshold, options.max_prop_missing)


def _prepare_optparser():
    """Prepare optparser object. New options will be added in this
    function first.
    """
    usage = """usage: %prog"""
    description = """"""

    optparser = OptionParser(description=description, usage=usage, add_help_option=False)
    optparser.add_option("-h", "--help", action="help", help="show this help message and exit.")
    optparser.add_option("-i", "--input_vcf_file", dest="input_vcf_file", type="string",
                         help="Path to input vcf file where the SNPs are located. Default: %default")
    optparser.add_option("-o", "--output_file", dest="output_file", type="string",
                         help="Path to the file that reformated output. Default: %default")
    optparser.add_option("-s", "--sex_file", dest="sex_file", type="string", default=None,
                         help="Path to the file that contains sex information. Default: %default")
    optparser.add_option("-m", "--mother_name", dest="mother_name", type="string", default=None,
                         help="name of the mother for that family. Default: %default")
    optparser.add_option("-f", "--father_name", dest="father_name", type="string", default=None,
                         help="name of the father for that family. Default: %default")
    optparser.add_option("-n", "--family_name", dest="family_name", type="string", default=None,
                         help="Name of the family. Default: %default")
    optparser.add_option("-g", "--geno_qual_threshold", dest="geno_qual_threshold", type="int", default=20,
                         help="The genotype quality threshold above which genotypes will be used. Default: %default")
    optparser.add_option("-x", "--max_prop_missing", dest="max_prop_missing", type="float", default=.5,
                         help="The maximum of missing samples across all populations. Default: %default")
    optparser.add_option("--phased", dest="phased", action="store_true", default=False,
                         help="Use Phasing information (if available) to create larger markers. Default: %default")
    optparser.add_option("--debug", dest="debug", action="store_true", default=False,
                         help="Set the verbosity to debug mode. Default: %default")
    return optparser


def _verifyOption(options):
    """Check if the mandatory option are present in the options objects.
    @return False if any argument is wrong."""
    arg_pass = True

    if not options.input_vcf_file or not os.path.exists(options.input_vcf_file):
        logging.error("You must specify a valid input file.")
        arg_pass = False
    if not options.output_file or not os.path.exists(os.path.dirname(os.path.abspath(options.output_file))):
        logging.error("You must specify a valid output file.")
        arg_pass = False
    return arg_pass


if __name__ == "__main__":
    main()   

