#!/usr/bin/env python
'''
Created on Mar 9, 2011

@author: tcezard
'''
__version__='0.1'
from utils import utils_logging
from IO_interface import vcfIO
import sys, os
import pprint
from collections import Counter
import logging
from optparse import OptionParser


def generate_empty_hash_with_sample(all_samples):
    hash={}
    for sample in all_samples:
        hash[sample]='.'
    return hash

def test_mandatory_samples(vcf_records, genotype_quality_threshold, minimum_depth, mandatory_list_sample):
    """Test that all the mandatory samples have high quality genotype and that more than 1 allele is represented"""
    
    keep=False
    if mandatory_list_sample:
        genotypes = vcf_records.get_valid_genotype_per_sample(genotype_quality_threshold, minimum_depth, sample_list=mandatory_list_sample)
        #print vcf_records.get_reference(),vcf_records.get_position()
        if len(genotypes)>0:
            all_gen=set()
            all_samples=[]
            for genotype in genotypes.keys():
                all_gen.update(set(genotype.split('/')))
                all_samples.extend(genotypes.get(genotype))
            if len(all_gen)>1 and len(all_samples)==len(mandatory_list_sample):
                keep=True
    else:
        keep=True
    return keep 


def test_all_samples(vcf_records, genotype_quality_threshold, minimum_depth, list_samples, min_nb_high_qual_sample):
    keep=False
    genotypes = vcf_records.get_valid_genotype_per_sample(genotype_quality_threshold, minimum_depth, sample_list=list_samples)
    if len(genotypes)>0:
        nb_of_valid_genotype=0
        all_hap_count = Counter()
        all_gen = set()
        for genotype in genotypes.keys():
            for h in genotype.split('/') : all_hap_count[h]+=len(genotypes.get(genotype))
            all_gen.update(set(genotype.split('/')))
            nb_of_valid_genotype+=len(genotypes.get(genotype))
        #good = 0
        #for count in all_hap_count.values():
        #    if count > min_nb_high_qual_sample:
        #        good += 1
        #Test if the site is variable in the samples and that the are a minimal number of valid genotype
        if len(all_gen)>1 and nb_of_valid_genotype>=min_nb_high_qual_sample :
            keep=True
    return keep




def vcf_to_simple_genotype(vcf_file, mandatory_list_sample, list_samples, min_nb_high_qual_sample=1, print_all_genotype=False):
    file_handle = utils_logging.open_input_file(vcf_file, pipe=False)
    reader  = vcfIO.VcfReader(file_handle)
    if not list_samples:
        list_samples=reader.get_sample_names()
    discarded_parents=0
    nb_markers=0
    print '#chr\tpos\tal1\tal2\t%s'%('\t'.join(list_samples))
    for vcf_records in reader:
        #First check that the parent are callable
        nb_markers+=1
        
        keep1=test_mandatory_samples(vcf_records, genotype_quality_threshold=20, minimum_depth=6, mandatory_list_sample=mandatory_list_sample)
        
        keep2 = test_all_samples(vcf_records, genotype_quality_threshold=20, minimum_depth=6, list_samples=list_samples,
                                 min_nb_high_qual_sample=min_nb_high_qual_sample)
        
        if keep1 and keep2:
            ref_base = vcf_records.get_reference_base()
            alt_bases = vcf_records.get_alt_bases()
            if len(alt_bases)>1:
                continue
            else:
                alt_base=alt_bases[0]
            if not print_all_genotype:
                genotypes_all = vcf_records.get_valid_genotype_per_sample(genotype_quality_threshold=20, minimum_depth=6)
                samples2genotype = generate_empty_hash_with_sample(list_samples)
                
                for genotype in genotypes_all:
                    sample_list = genotypes_all.get(genotype)
                    genotype_str = genotype.replace('0', ref_base).replace('1', alt_base)
                    #genotype_str = genotype
                    
                    for sample in sample_list:
                        samples2genotype[sample] = genotype_str
                out=[]
                out.append(vcf_records.get_reference())
                out.append(str(vcf_records.get_position()))
                out.append(ref_base)
                out.append(alt_base)
                for sample in list_samples:
                    out.append(samples2genotype.get(sample))
                print '\t'.join(out)
            else:
                all_genotypes=vcf_records.get_all_genotype(list_samples)
                all_genotype_quals=vcf_records.get_all_genotype_quality(list_samples)
                all_sample_depth=vcf_records.get_all_sample_depth(list_samples)
                
                out=[]
                out.append(vcf_records.get_reference())
                out.append(str(vcf_records.get_position()))
                out.append(ref_base)
                out.append(alt_base)
                for i in range(len(all_genotypes)):
                    if all_genotypes[i]:
                        genotype_str = all_genotypes[i].replace('0', ref_base).replace('1', alt_base)
                    else:
                        genotype_str='.'
                    if not all_genotype_quals[i]: all_genotype_quals[i]=0
                    if not all_sample_depth[i]: all_sample_depth[i]=0
                    out.append('%s:%s:%s'%(genotype_str,all_genotype_quals[i],all_sample_depth[i]))
                print '\t'.join(out)
        else:
            discarded_parents+=1
    sys.stderr.write("%s markers %s filtered out\n"%(nb_markers,discarded_parents))
    return True



def main():
    #initialize the logging
    utils_logging.init_logging()
    #Setup options
    optparser=_prepare_optparser()
    (options,args) = optparser.parse_args()
    #verify options
    arg_pass=_verifyOption(options)
    utils_logging.change_log_stdout_to_log_stderr()
    if not arg_pass:
        logging.warning(optparser.get_usage())
        logging.critical("Non valid arguments: exit")
        sys.exit(1)
    if options.mandatory_list_samples:
        mandatory_list_samples = options.mandatory_list_samples.split(',')
    else:
        mandatory_list_samples = None
    if options.list_samples:
        list_samples = options.list_samples.split(',')
    else:
        list_samples = None
    run_fine=vcf_to_simple_genotype(options.vcf_file, mandatory_list_samples, list_samples,
                                    options.nb_of_required_high_qual_genotype, options.print_all_genotype)
    if run_fine:
        logging.info('Run completed')
    else:
        logging.error('Run Failed')
        sys.exit(1)
        
def _prepare_optparser():
    """Prepare optparser object. New options will be added in this
    function first.
    """
    usage = """usage: %prog """
    description = """This script will take a vcf file and change the format to a tab delimited file"""
    
    prog_version=__version__
    optparser = OptionParser(version="%prog v"+prog_version,description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="show this help message and exit.")
    optparser.add_option("-i","--vcf_file",dest="vcf_file",type="string",
                         help="Path to a vcf file where the data is located. Default: %default")
    optparser.add_option("-m","--mandatory_list_samples",dest="mandatory_list_samples",type="string",default='',
                         help="Comma separated list of sample that need to be genotyped and variable for a SNPs ot be kept. Default: %default")
    optparser.add_option("-l","--list_samples",dest="list_samples",type="string",default='',
                         help="Comma separated list of sample name that will be output. Default: %default")
    optparser.add_option("-n","--nb_of_required_high_qual_genotype",dest="nb_of_required_high_qual_genotype",type="int",default=1,
                         help="Minimum number of high quality genotype to output a SNPs. Default: %default")
    optparser.add_option("-a", "--print_all_genotype",dest="print_all_genotype",action='store_true',default=False,
                         help="Output all genotypes along with genotype quality and depth information. Default: %default")
    return optparser


def _verifyOption(options):
    """Check if the mandatory option are present in the options objects.
    @return False if any argument is wrong."""
    arg_pass=True
    
    if not options.vcf_file :
        logging.error("You must specify a vcf file.")
        arg_pass=False
    elif not os.path.exists(options.vcf_file):
        logging.error("Vcf file (%s) not found: You must specify a valid file path."%(options.vcf_file))
        arg_pass=False
    return arg_pass



if __name__=="__main__":
    main()

if __name__=="1__main__":
    input_file = sys.argv[1]
    mandatory_list_sample = sys.argv[2].split(',')
    list_samples = sys.argv[3].split(',')
    if len(sys.argv)>4: all=True
    else: all=False
    vcf_to_simple_genotype(input_file, mandatory_list_sample, list_samples, all)

if __name__=="1__main__":
    input_file = sys.argv[1]
    #mandatory_list_sample = sys.argv[2].split(',')
    #list_samples = sys.argv[3].split(',')
    #if len(sys.argv)>4: all=True
    #else: all=False
    vcf_to_simple_genotype(input_file, None, None, all=False)
