#!/usr/bin/env python
'''
Created on Mar 9, 2011

@author: tcezard
'''
from utils import utils_logging
from IO_interface import vcfIO
import sys
import logging
from utils import atgc2iupac
from optparse import OptionParser
import os
from collections import Counter,defaultdict
import getpass
import time
import re

OUTPUT_TYPE_PHYLIP='phylip'
OUTPUT_TYPE_STRUCTURE='structure'
OUTPUT_TYPE_GENEPOP='genepop'
OUTPUT_TYPE=[OUTPUT_TYPE_PHYLIP,OUTPUT_TYPE_STRUCTURE,OUTPUT_TYPE_GENEPOP]

def read_pop_file(pop_file):
    sample2pop={}
    pop2sample=defaultdict(list)
    with open(pop_file) as open_file:
        for line in open_file:
            sp_line = line.strip().split()
            if len(sp_line)>1:
                pop = sp_line[1]
                sample = sp_line[0]
                sample2pop[sample]=pop
                pop2sample[pop].append(sample)
    return sample2pop, pop2sample


def vcf_to_phylip(vcf_file, output_file, pop_file=None, genotype_quality_threshold=20 , max_prop_missing = .5, phased=False):
    
    file_handle = utils_logging.open_input_file(vcf_file, pipe=False)
    if phased:
        #Not useful for phylip format
        reader  = vcfIO.VcfReader(file_handle)
    else:
        reader  = vcfIO.VcfReader(file_handle)
    all_samples_in_file=reader.get_sample_names()
    if pop_file:
        sample2pop, pop2sample = read_pop_file(pop_file)
        all_samples_in_file=sample2pop.keys()
    
    max_missing =int( len(all_samples_in_file) * max_prop_missing)

    all_lines={}
    for sample in all_samples_in_file:
        all_lines[sample]=[]
    nb_sample = len(all_samples_in_file)
    nb_sequence = 0
    count_more_than_one_allele=0
    count_indel=0
    count_too_many_missing=0
    count_non_polymorphic=0
    nb_missing_per_sample=Counter()
    for vcf_records in reader:
        ref_base=vcf_records.get_reference_base()
        alt_bases=vcf_records.get_alt_bases()
        array = [ref_base]
        array.extend(alt_bases)
        if len(alt_bases)>1:
            count_more_than_one_allele+=1
            continue
        if vcf_records.is_indel():
            count_indel+=1
            continue
        
        nb_missing = 0
        set_missing=set() 
        all_chars = []
        all_codes=set()
        for sample in all_samples_in_file:
            gt = vcf_records.get_genotype(sample)
            gq = vcf_records.get_genotype_quality(sample)
            if gt and gq>genotype_quality_threshold:
                value1,value2 = re.split('[/|]',gt)
                code = atgc2iupac(list(set([array[int(value1)], array[int(value2)]])))
                all_codes.add(code)
            else:
                nb_missing+=1
                set_missing.add(sample)
                code = '?'
            all_chars.append(code)
        if len(all_codes)==1:
            count_non_polymorphic+=1
            continue
        if nb_missing <= max_missing:
            nb_sequence+=1
            
            for s in set_missing: nb_missing_per_sample[s]+=1
            for i,sample in enumerate(all_samples_in_file):
                all_lines[sample].append(all_chars[i])
        else:
            count_too_many_missing+=1
    if count_more_than_one_allele:
        logging.warning("%s snps remove because they had more than 2 alleles"%(count_more_than_one_allele))
    if count_indel:
        logging.warning("%s indels removed"%(count_indel))
    if count_non_polymorphic:
        logging.warning("%s snps removed because no polymorphism was found between populations"%(count_non_polymorphic))
    if count_too_many_missing:
        logging.warning("%s snps removed because >%s missing samples"%(count_too_many_missing, max_missing))
    logging.info("%s snps output in %s format"%(nb_sequence, OUTPUT_TYPE_PHYLIP))
    for sample in nb_missing_per_sample:
        logging.info("%s markers missing in %s"%(nb_missing_per_sample.get(sample), sample))
    #Output the phylip formated file
    with open(output_file,'w') as open_output:
        open_output.write('%s %s\n'%(nb_sample, nb_sequence))
        for sample in all_samples_in_file:
            open_output.write("%s %s\n"%(sample,''.join(all_lines.get(sample))) )
        


def vcf_to_structure(vcf_file, output_file, pop_file, genotype_quality_threshold=20 , max_prop_missing = .5, phased=False):
    
    file_handle = utils_logging.open_input_file(vcf_file, pipe=False)
    if phased:
        reader  = vcfIO.PhasedVcfReader(file_handle)
    else:
        reader  = vcfIO.VcfReader(file_handle)
    all_samples_in_file=reader.get_sample_names()
    
    
    if pop_file:
        sample2pop, pop2sample = read_pop_file(pop_file)
        pops = pop2sample.keys()
        all_samples = sample2pop.keys()
    else:
        pop = 'dummy_pop'
        pops = [pop]
        all_samples=all_samples_in_file
        pop2sample={pop:[all_samples_in_file]}
        sample2pop={}
        for sample in all_samples:sample2pop[sample]=pop
    
    sample_errors =  set(all_samples).difference(set(all_samples_in_file))
    if len(sample_errors)>0:
        logging.critical('%s samples (%s) from the population file not found in the vcf file'%(len(sample_errors), ', '.join(sample_errors)))
        return -2
    
    all_lines={}
    headers = []
    for sample in all_samples:
        all_lines[sample+'1']=[]
        all_lines[sample+'2']=[]
    nb_sample = len(all_samples)
    max_missing =int( len(all_samples) * max_prop_missing)
    nb_sequence = 0
    count_more_than_one_allele=0
    count_indel=0
    count_too_many_missing=0
    count_non_polymorphic=0
    for vcf_records in reader:
        ref_base=vcf_records.get_reference_base()
        alt_bases=vcf_records.get_alt_bases()
        
        if len(alt_bases)>1:
            count_more_than_one_allele+=1
            continue
        if vcf_records.is_indel():
            count_indel+=1
            continue
        
        nb_missing = 0 
        all_alleles1 = [] 
        all_alleles2 = []
        all_genotypes=set()
        for sample in all_samples:
            gt = vcf_records.get_genotype(sample)
            gq = vcf_records.get_genotype_quality(sample)
            if gt and gq>genotype_quality_threshold:
                allele1, allele2 = gt.split('/')
                all_genotypes.add(gt)
            else:
                nb_missing+=1
                allele1, allele2 =('-9','-9')
            all_alleles1.append(allele1)
            all_alleles2.append(allele2)
        
        if len(all_genotypes)==1:
            count_non_polymorphic+=1
            continue
        
        if nb_missing <= max_missing:
            nb_sequence+=1
            variant_name = '%s:%s'%(vcf_records.get_reference(),vcf_records.get_position())
            headers.append(variant_name)
            for i,sample in enumerate(all_samples):
                all_lines[sample+'1'].append(all_alleles1[i])
                all_lines[sample+'2'].append(all_alleles2[i])
        else:
            count_too_many_missing+=1
    if count_more_than_one_allele:
        logging.warning("%s snps remove because they had more than 2 alleles"%(count_more_than_one_allele))
    if count_indel:
        logging.warning("%s indels removed"%(count_indel))
    if count_non_polymorphic:
        logging.warning("%s snps removed because no polymorphism was found between populations"%(count_non_polymorphic))
    if count_too_many_missing:
        logging.warning("%s snps removed because >%s missing samples"%(count_too_many_missing, max_missing))
    logging.info("%s samples and %s SNPs output"%(nb_sample, nb_sequence))
    with open(output_file,'w') as open_file:
        open_file.write('\t%s\n'%('\t'.join(headers)))
        for sample in all_samples:
            open_file.write("%s\t%s\t%s\n"%(sample, pops.index(sample2pop.get(sample))+1, '\t'.join(all_lines.get(sample+'1')) ) )
            open_file.write("%s\t%s\t%s\n"%(sample, pops.index(sample2pop.get(sample))+1, '\t'.join(all_lines.get(sample+'2')) ) )
        
        
def vcf_2_genepop(vcf_file, output_file, pop_file, genotype_quality_threshold=20 , max_prop_missing = .5, phased=False):
        
    file_handle = utils_logging.open_input_file(vcf_file, pipe=False)
    if phased:
        reader  = vcfIO.PhasedVcfReader(file_handle)
    else:
        reader  = vcfIO.VcfReader(file_handle)
    all_samples_in_file=reader.get_sample_names()
    
    if pop_file:
        sample2pop, pop2samples = read_pop_file(pop_file)
        all_samples = sample2pop.keys()
    else:
        pop = 'dummy_pop'
        all_samples=all_samples_in_file
        pop2samples={pop:all_samples_in_file}
        sample2pop={}
        for sample in all_samples:sample2pop[sample]=pop    
    
    all_lines={}
    headers = []
    for sample in all_samples:
        all_lines[sample]=[]
    nb_sample = len(all_samples)
    max_missing =int( len(all_samples) * max_prop_missing)
    nb_sequence = 0
    count_more_than_one_allele=0
    count_indel=0
    count_too_many_missing=0
    count_non_polymorphic=0
    for vcf_record in reader:
        ref_base=vcf_record.get_reference_base()
        alt_bases=vcf_record.get_alt_bases()
        
        if len(alt_bases)>1:
            count_more_than_one_allele+=1
            continue
        if vcf_record.is_indel():
            count_indel+=1
            continue
        
        nb_missing = 0 
        all_alleles = []
        all_genotypes = set()
        for sample in all_samples:
            gt = vcf_record.get_genotype(sample)
            gq = vcf_record.get_genotype_quality(sample)
            if gt and gq>genotype_quality_threshold:
                allele1, allele2 = gt.split('/')
                alleles="%02d%02d"%(int(allele1)+1, int(allele2)+1)
                all_genotypes.add(gt)
            else:
                nb_missing+=1
                alleles ='0000'
            all_alleles.append(alleles)
        
        if len(all_genotypes)==1:
            count_non_polymorphic+=1
            continue    
        
        if nb_missing <= max_missing:
            nb_sequence+=1
            variant_name = '%s:%s'%(vcf_record.get_reference(),vcf_record.get_position())
            headers.append(variant_name)
            for i,sample in enumerate(all_samples):
                all_lines[sample].append(all_alleles[i])
        else:
            count_too_many_missing+=1
    title_line="Generated by %s from %s on %s"%(getpass.getuser(), vcf_file, time.ctime())
    if count_more_than_one_allele:
        logging.warning("%s snps remove because they had more than 2 alleles"%(count_more_than_one_allele))
    if count_indel:
        logging.warning("%s indels removed"%(count_indel))
    if count_non_polymorphic:
        logging.warning("%s snps removed because no polymorphism was found between populations"%(count_non_polymorphic))
    if count_too_many_missing:
        logging.warning("%s snps removed because >%s missing samples"%(count_too_many_missing, max_missing))
    logging.info("%s samples and %s SNPs output"%(nb_sample, nb_sequence))
    with open(output_file,'w') as open_file:
        open_file.write(title_line+'\n')
        open_file.write('%s\n'%('\n'.join(headers)))
        for pop in pop2samples:
            open_file.write("Pop\n")
            samples = pop2samples.get(pop)
            for sample in samples:
                open_file.write("%s %s, %s\n"%(sample, pop, '\t'.join(all_lines.get(sample)) ) )
        

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
    if options.debug:
        utils_logging.init_logging(logging.DEBUG)
    if options.type_output==OUTPUT_TYPE_PHYLIP:
        vcf_to_phylip(options.input_vcf_file, options.output_file, options.pop_file, options.geno_qual_threshold, options.max_prop_missing)
    elif options.type_output==OUTPUT_TYPE_STRUCTURE:
        vcf_to_structure(options.input_vcf_file, options.output_file, options.pop_file, options.geno_qual_threshold, options.max_prop_missing, options.phased)
    elif options.type_output==OUTPUT_TYPE_GENEPOP:
        vcf_2_genepop(options.input_vcf_file, options.output_file, options.pop_file, options.geno_qual_threshold, options.max_prop_missing, options.phased)


def _prepare_optparser():
    """Prepare optparser object. New options will be added in this
    function first.
    """
    usage = """usage: %prog -i <input_vcf_file> -o <output_reformated> [-t phylip -p <pop_info_file> -g 20 -x .5]"""
    description = """This script will take a vcf file and reformat the SNPs into another output. It will also filter some of the SNPs out if they do not fulfil the provided criteria"""
    
    optparser = OptionParser(description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="show this help message and exit.")
    optparser.add_option("-i","--input_vcf_file",dest="input_vcf_file",type="string",
                         help="Path to input vcf file where the SNPs are located. Default: %default")
    optparser.add_option("-o","--output_file",dest="output_file",type="string",
                         help="Path to the file that reformated output. Default: %default")
    optparser.add_option("-p","--pop_file",dest="pop_file",type="string",default=None,
                         help="Path to the file that contains population information. Default: %default")
    optparser.add_option("-g","--geno_qual_threshold",dest="geno_qual_threshold",type="int",default =20,
                         help="The genotype quality threshold above which genotypes will be used. Default: %default")
    optparser.add_option("-x","--max_prop_missing",dest="max_prop_missing",type="float",default=.5,
                         help="The maximum of missing samples across all populations. Default: %default")
    optparser.add_option("-t","--type_output",dest="type_output",type="string",default=OUTPUT_TYPE_PHYLIP,
                         help="Type of output format required. (Should be one of "+ ", ".join(OUTPUT_TYPE)+") Default: %default")
    optparser.add_option("--phased",dest="phased",action="store_true",default=False,
                         help="Use Phasing information (if available) to create larger markers. Default: %default")
    optparser.add_option("--debug",dest="debug",action="store_true",default=False,
                         help="Set the verbosity to debug mode. Default: %default")
    return optparser


def _verifyOption(options):
    """Check if the mandatory option are present in the options objects.
    @return False if any argument is wrong."""
    arg_pass=True
    
    if not options.input_vcf_file or not os.path.exists(options.input_vcf_file):
        logging.error("You must specify a valid input file.")
        arg_pass=False
    if not options.output_file or not os.path.exists(os.path.dirname(os.path.abspath(options.output_file))):
        logging.error("You must specify a valid output file.")
        arg_pass=False
    #if not options.pop_file or not os.path.exists(options.pop_file):
    #    logging.error("You must specify a valid population information file.")
    #    arg_pass=False
    if not options.type_output in OUTPUT_TYPE:
        logging.error("You must specify a valid output type (%s)."%", ".join(OUTPUT_TYPE))
    #elif (options.type_output == OUTPUT_TYPE_STRUCTURE) and \
    #not (options.pop_file and os.path.exists(options.pop_file)):
    #    logging.error("You must specify a population file with -p for "+OUTPUT_TYPE_STRUCTURE+" format output.")
    #    arg_pass=False
    return arg_pass


if __name__=="__main__":
    main()   

