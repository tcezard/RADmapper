'''
Created on Mar 9, 2011

@author: tcezard
'''
import sys
import logging
from optparse import OptionParser
from collections import defaultdict

from utils import utils_logging
from scipy import stats
import numpy


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

def generate_empty_hash_with_sample(all_samples):
    hash={}
    for sample in all_samples:
        hash[sample]='.'
    return hash

def get_normalize_coverage(coverage_file, nb_sample_required=0):
    #contig    coverage        coverage_mrk_dup        nb_sample

    file_handle = utils_logging.open_input_file(coverage_file, pipe=False)
    all_samples_to_coverage={}
    all_samples=[]
    all_markers=[]
    for line in file_handle:
        sp_line = line.strip().split("\t")
        if line.startswith("#"):
            for i in range(4, len(sp_line), 2):
                sample = sp_line[i]
                all_samples.append(sample)
                all_samples_to_coverage[sample]=[]
        elif int(sp_line[3])>=nb_sample_required:
            i=0
            all_markers.append(sp_line[0])
            j=0
            for i in range(4, len(sp_line), 2):
                data=sp_line[i]
                all_samples_to_coverage[all_samples[j]].append(int(data))
                j+=1
    all_samples_to_norm_coverage={}
    for sample in all_samples:
        coverage_info = all_samples_to_coverage.get(sample)
        total=sum(coverage_info)
        normalized_coverage=[]
        for coverage in coverage_info:
            normalized_coverage.append(float(coverage)/total*1000000)
        all_samples_to_norm_coverage[sample]=normalized_coverage
    
    return all_markers, all_samples, all_samples_to_norm_coverage

def detect_presence_absence_markers(pop_file , coverage_file, nb_sample_required=0):
    sample2pop, pop2sample = read_pop_file(pop_file)
    if len(pop2sample)!=2:
        logging.critical('Hemizygous Markers can only be search between two set of samples. Edit your population file to have two populations')
        return -1
    pop1,pop2 = pop2sample.keys()
    samples_pop1 = pop2sample.get(pop1)
    samples_pop2 = pop2sample.get(pop2)
    
    all_markers, all_samples, all_samples_to_norm_coverage = get_normalize_coverage(coverage_file, nb_sample_required)
    
    sample_errors =  set(sample2pop.keys()).difference(set(all_samples)) 
    if len(sample_errors)>0:
        logging.critical('%s samples (%s) from the population file not found in the coverage file'%(len(sample_errors), ', '.join(sample_errors)))
        return -2
    for i, marker in enumerate(all_markers):
        out=[marker]
        out_pop1=[]
        out_pop2=[]
        pop1_values = []
        pop2_values = []
        for sample in samples_pop1:
            cov = all_samples_to_norm_coverage.get(sample)[i]
            pop1_values.append(cov)
            out_pop1.append(str(cov))
        
        for sample in samples_pop2:
            cov = all_samples_to_norm_coverage.get(sample)[i]
            pop2_values.append(cov)
            out_pop2.append(str(cov))
        pop1_nvalues = numpy.array(pop1_values)
        pop1_nvalues_2 =pop1_nvalues*2
        pop2_nvalues = numpy.array(pop2_values)
        pop2_nvalues_2 = pop2_nvalues*2
        t_stat, pvalue =  stats.ttest_ind(pop1_nvalues,pop2_nvalues)
        if pvalue<.05:
            t_stat_comp1, pvalue_comp1 =  stats.ttest_ind(pop1_nvalues_2,pop2_nvalues)
            if pvalue_comp1 >0.5:
                fold=pop2_nvalues.mean()/pop1_nvalues.mean()
                if fold < 2.2 and fold >1.8: 
                    print ' '.join(out), fold, pop1_nvalues.mean(), pop2_nvalues.mean(), pvalue, pvalue_comp1

def detect_hemizygous_markers(pop_file , coverage_file, nb_sample_required=0):
    sample2pop, pop2sample = read_pop_file(pop_file)
    if len(pop2sample)!=2:
        logging.critical('Hemizygous Markers can only be search between two set of samples. Edit you population file to have two populations')
        return -1
    pop1,pop2 = pop2sample.keys()
    samples_pop1 = pop2sample.get(pop1)
    samples_pop2 = pop2sample.get(pop2)
    
    all_markers, all_samples, all_samples_to_norm_coverage = get_normalize_coverage(coverage_file, nb_sample_required)
    
    sample_errors =  set(sample2pop.keys()).difference(set(all_samples)) 
    if len(sample_errors)>0:
        logging.critical('%s samples (%s) from the population file not found in the coverage file'%(len(sample_errors), ', '.join(sample_errors)))
        return -2
    print "#contigs\tmean_male\tsd_male\tmean_female\tsd_female\t%s\t%s"%("\t".join(samples_pop1), "\t".join(samples_pop2))
    for i, marker in enumerate(all_markers):
        out=[marker]
        out_pop1=[]
        out_pop2=[]
        pop1_values = []
        pop2_values = []
        for sample in samples_pop1:
            cov = all_samples_to_norm_coverage.get(sample)[i]
            pop1_values.append(cov)
            out_pop1.append(str(cov))
        
        for sample in samples_pop2:
            cov = all_samples_to_norm_coverage.get(sample)[i]
            pop2_values.append(cov)
            out_pop2.append(str(cov))
        pop1_nvalues = numpy.array(pop1_values)
        pop1_nvalues_2 =pop1_nvalues*2
        pop2_nvalues = numpy.array(pop2_values)
        pop2_nvalues_2 = pop2_nvalues*2
        t_stat, pvalue =  stats.ttest_ind(pop1_nvalues,pop2_nvalues)
        if pvalue<.05:
            t_stat_comp1, pvalue_comp1 =  stats.ttest_ind(pop1_nvalues_2,pop2_nvalues)
            if pvalue_comp1 >0.5:
                fold=pop2_nvalues.mean()/pop1_nvalues.mean()
                if fold < 2.2 and fold >1.8: 
                    print ' '.join(out), fold, pop1_nvalues.mean(), pop2_nvalues.mean(), pvalue, pvalue_comp1
        
        #out.append(str(male_values.get_mean()))
        #out.append(str(male_values.get_std_dev()))
        #out.append(str(female_values.get_mean()))
        #out.append(str(female_values.get_std_dev()))
        
        #print "%s\t%s\t%s"%("\t".join(out),"\t".join(out_male),"\t".join(out_female)) 

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
    #detect_hemizygous_markers(options.pop_file, options.coverage_file)
    detect_presence_absence_markers(options.pop_file, options.coverage_file)

def _prepare_optparser():
    """Prepare optparser object. New options will be added in this
    function first.
    """
    usage = """usage: %prog -i <input_vcf_file> -o <output_reformated> [-t phylip -p <pop_info_file> -g 20 -x .5]"""
    description = """This script will take a vcf file and reformat the SNPs into another output. It will also filter some of the SNPs out if they do not fulfil the provided criteria"""
    
    optparser = OptionParser(description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="show this help message and exit.")
    optparser.add_option("-c","--coverage_file",dest="coverage_file",type="string",
                         help="Path to coverage file where read1 coverage is stored for all samples. Default: %default")
    optparser.add_option("-o","--output_file",dest="output_file",type="string",
                         help="Path to the file that reformated output. Default: %default")
    optparser.add_option("-p","--pop_file",dest="pop_file",type="string",
                         help="Path to the file that contains population information. Default: %default")
    optparser.add_option("--debug",dest="debug",action="store_true",default=False,
                         help="Set the verbosity to debug mode. Default: %default")
    return optparser


def _verifyOption(options):
    """Check if the mandatory option are present in the options objects.
    @return False if any argument is wrong."""
    arg_pass=True
    
    return arg_pass


if __name__=="__main__":
    main()   

