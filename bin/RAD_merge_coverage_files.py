#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 31 Oct 2013
@author: tcezard
'''
import sys
import os
import logging
from optparse import OptionParser
from collections import Counter

from utils import utils_logging


header_to_record=["length","known_novel"]
header_to_ignore = ["#locus", "#contig", "coverage", "coverage_mrk_dup", "nb_allele", "nb_sample"]


def parse_coverage_file(coverage_file, all_loci_per_sample={"all_sample": Counter()},
                        all_loci_mrk_dup_per_sample={"all_sample": Counter()},
                        all_loci_characteristic={}, all_loci=set(), all_samples=set(),
                        all_characteristics=set()):
    samples_index={}
    characteristic_index={}
    with open(coverage_file) as open_file:
        for line in open_file:
            if line.startswith('#'):
                sp_line = line.strip().split('\t')
                for i, element in enumerate(sp_line):
                    if element in header_to_ignore:
                        continue
                    if element in header_to_record:
                        characteristic_index[i]=element
                        all_loci_characteristic[element]={}
                        all_characteristics.add(element)
                    else:
                        samples_index[i]=element
                        if element.endswith('mrk_dup'):
                            if not all_loci_mrk_dup_per_sample.has_key(element):
                                all_loci_mrk_dup_per_sample[element]=Counter()
                        elif not all_loci_per_sample.has_key(element):
                            all_loci_per_sample[element]=Counter()
                            all_samples.add(element)
            else:
                sp_line = line.strip().split('\t')
                all_loci.add(sp_line[0])
                for i, element in enumerate(sp_line):
                    if i in samples_index.keys():
                        sample = samples_index.get(i)
                        if sample.endswith("mrk_dup"):
                            all_loci_mrk_dup_per_sample[sample][sp_line[0]]+=int(sp_line[i])
                            all_loci_mrk_dup_per_sample["all_sample"][sp_line[0]]+=int(sp_line[i])
                        else:
                            all_loci_per_sample[sample][sp_line[0]]+=int(sp_line[i])
                            all_loci_per_sample["all_sample"][sp_line[0]]+=int(sp_line[i])
                    elif i in characteristic_index:
                        all_loci_characteristic[characteristic_index.get(i)][sp_line[0]]=sp_line[i]
    return (all_loci_per_sample, all_loci_mrk_dup_per_sample, all_loci_characteristic,
             all_loci, all_samples, all_characteristics)


def merge_coverage(coverage_files, output_file):
    all_loci_per_sample={"all_sample":Counter()}
    all_loci_mrk_dup_per_sample={"all_sample":Counter()}
    all_loci_characteristic={}
    all_loci=set()
    all_samples=set()
    all_characteristics=set()
    for file in coverage_files:
        sys.stderr.write("load %s\n"%file)
        parse_coverage_file(file, all_loci_per_sample, all_loci_mrk_dup_per_sample,
                            all_loci_characteristic, all_loci, all_samples, all_characteristics)


    #We're done reading now onto writing
    header=['#locus','coverage','coverage_mrk_dup']
    header.extend(sorted(all_characteristics))
    header.append('nb_sample')
    header.extend(['%s\t%s_mrk_dup'%(sample,sample) for sample in sorted(all_samples)])
    with open(output_file,'w') as open_output:
        open_output.write('\t'.join(header  )+'\n')
        for loci in sorted(all_loci):
            final_out=[loci]
            final_out.append(str(all_loci_per_sample["all_sample"].get(loci,0)))
            final_out.append(str(all_loci_mrk_dup_per_sample["all_sample"].get(loci,0)))
            for characteristics in sorted(all_characteristics):
                final_out.append(all_loci_characteristic[characteristics].get(loci,""))
            out=[]
            nb_sample=0
            for sample in sorted(all_samples):
                out.append(str(all_loci_per_sample[sample].get(loci,0)) )
                nb_dup=all_loci_mrk_dup_per_sample[sample+"_mrk_dup"].get(loci,0)
                if nb_dup>2:
                    nb_sample+=1
                out.append(str(nb_dup))
            final_out.append(str(nb_sample))
            final_out.extend(out)
            open_output.write('\t'.join(final_out)+'\n')


def main():
    #initialize the logging
    utils_logging.init_logging(logging.INFO)
    #Setup options
    optparser=_prepare_optparser()
    (options,args) = optparser.parse_args()
    #verify options
    arg_pass=_verifyOption(options)
    if not arg_pass:
        logging.warning(optparser.get_usage())
        logging.critical("Non valid arguments: exit")
        sys.exit(1)
    coverage_files=[options.coverage_files]
    if len(args)>0:
        coverage_files.extend(args)
    merge_coverage(coverage_files,options.output_file)
    

def _prepare_optparser():
    """Prepare optparser object. New options will be added in this
    function first.
    """
    usage = """usage: %prog <-b bam_file> [ -o output_file]"""
    description = """This script will take aligned RAD read to the consensuses and calculate per consensus coverage."""
    
    optparser = OptionParser(version="None",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="show this help message and exit.")
    optparser.add_option("-c","--coverage_files",dest="coverage_files",type="string",
                         help="Path to files that will be used to get the coverage information to merge. Default: %default")
    optparser.add_option("-o","--output_file",dest="output_file",type="string",
                         help="Path to the output file the merged results will be stored. Default: %default")
    return optparser


def _verifyOption(options):
    """Check if the mandatory option are present in the options objects.
    @return False if any argument is wrong."""
    arg_pass=True
    if not options.coverage_files or not os.path.exists(options.coverage_files):
        logging.error("Coverage file can't be found. please specify an existing coverage file with -c")
        arg_pass=False
    if not options.output_file or not os.path.exists(os.path.dirname(os.path.abspath(options.coverage_files))):
        logging.error("Output file can't be found or parent does not exist. please specify a valid output file with -o")
        arg_pass=False
    return arg_pass


if __name__=="__main__":
    main()
    
    