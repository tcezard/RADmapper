# -*- coding: utf-8 -*-
'''
Created on 25 April 2012
@author: tcezard
'''
import sys, os, re
from utils import utils_logging
import logging, threading
from optparse import OptionParser
from glob import glob
import command_runner
from overlap import merge_ranges, get_length_of_list_range
from utils.FastaFormat import FastaReader
import overlap
from collections import Counter


all_sites_headers=["coverage", "span_region","span_length", "alignment_length", "overlap_length", "number_contig", "number_contig_matching",
                   "max_length", "max_ol_length", "mismatches", "contig_coverage1", "contig_coverage2", "contig_coverage3", "contig_coverage4plus"]

def run_blast(contig_file, genome_file):
    blastn_plus_bin='/ifs/software/linux_x86_64/blast+/current/bin/blastn'
    output_file='%s.blast6out'%contig_file
    command='%s -query %s -db %s -max_target_seqs 1 -outfmt 6 -out %s'%(blastn_plus_bin, contig_file, genome_file, output_file)
    return_code = command_runner.run_command(command)
    if return_code!=0:
        return None
    return output_file


def read_site_info_file(input_file):
    open_file = open(input_file)
    all_sites={}
    for line in open_file:
        sp_line = line.strip().split()
        span_length=0
        for region in sp_line[2:]:
            start, end = region.split('-')
            span_length+=int(end)-int(start)+1
        all_sites[sp_line[0]]={"coverage":sp_line[1], "span_region":';'.join(sp_line[2:]), "span_length":span_length,
                               "contig_span":[],"alignment_length":0, "number_contig":0, "number_contig_matching":0,
                               "max_length":0, "max_ol_length":0, "mismatches":0, "overlap_length":0, "contig_coverage1":0, 
                               "contig_coverage2":0, "contig_coverage3":0,"contig_coverage4plus":0}
    return all_sites

def get_assessment(contig_file, all_sites,  genome_file):
    output_file = run_blast(contig_file, genome_file)
    if output_file:
        for contig_name, percent_identity, alignment_length, mismatches, contig_span in read_blast_fm6_output(output_file):
            match=re.match('(.+)_pair_\d+_length_\d+',contig_name)
            site_name = match.group(1)
            #all_sites[site_name]["alignment_length"]+=alignment_length
            all_sites[site_name]["contig_span"].append(contig_span)
            all_sites[site_name]["number_contig_matching"]+=1
            all_sites[site_name]["mismatches"]+=mismatches
        return all_sites
    else:
        return None
    
def get_basic_stats(contig_file,all_sites):
    open_file = open(contig_file)
    fasta_reader = FastaReader(open_file)
    for fasta_record in fasta_reader:
        header, sequence=fasta_record
        match=re.match('(.+)_pair_\d+_length_\d+',header)
        site_name = match.group(1)
        all_sites[site_name]["number_contig"]+=1
        if len(sequence)>all_sites[site_name].get("max_length"):
            all_sites[site_name]["max_length"]=len(sequence)
    return all_sites

def read_blast_fm6_output(blastout_file):
    open_file = open(blastout_file)
    for line in open_file:
        if line.startswith('#'):
            continue
        sp_line=line.strip().split()
        contig_name = sp_line[0]
        #reference_name = sp_line[1]
        percent_identity = float(sp_line[2])
        alignment_length = int(sp_line[3])
        mismatches = int(sp_line[4])
        #gap_opens = sp_line[5]
        #c_start = sp_line[6]
        #c_end = sp_line[7]
        r_start = int(sp_line[8])
        r_end = int(sp_line[9])
        if r_start>r_end:
            tmp=r_end
            r_end=r_start
            r_start=tmp
        yield(contig_name,percent_identity,alignment_length,mismatches,'%s-%s'%(r_start,r_end))
    

def output_all_sites(all_sites, output_sites):
    open_file = utils_logging.open_output_file(output_sites)
    open_file.write("sites\t%s\n"%("\t".join(all_sites_headers)))
    for site_name in all_sites.keys():
        open_file.write("%s\t%s\n"%(site_name,"\t".join([str(all_sites.get(site_name).get(key)) for key in all_sites_headers])))
    open_file.close()
    
    
def summarize(all_sites):
    for site_name in all_sites.keys():
        site_info=all_sites.get(site_name)
        contig_span_list = []
        for value in site_info.get("contig_span"):
            s,e = value.split('-')
            contig_span_list.append((int(s),int(e)))
        merged_contig_span_list=merge_ranges(contig_span_list)
        regions = site_info.get("span_region").split(';')
        overlap_ranges=[]
        pileup=Counter()
        #Find overlap between the contigs and the expected position
        for start2,end2 in contig_span_list:
            tmp_overlap_ranges=[]
            for i in range(int(start2), int(end2)+1):
                pileup[i]+=1
            for region in regions:
                if len(region)==0:
                    logging.warning("Site %s has no expected regions"%(site_name))
                    continue
                start1,end1 =region.split('-')
                
                overlap_res = overlap.get_overlap_if_exist(int(start1), int(end1), int(start2), int(end2))
                if overlap_res:
                    tmp_overlap_ranges.append(overlap_res)
                    s, e = overlap_res
            if len(tmp_overlap_ranges)==0: 
                logging.warning("%s: No overlap found between expected region %s and contig alignment %s-%s"%(site_name, site_info.get("span_region") ,start2, end2))
            else:
                max_ol_length = get_length_of_list_range(tmp_overlap_ranges)
                if site_info["max_ol_length"]<max_ol_length:
                    site_info["max_ol_length"]=max_ol_length
                overlap_ranges.extend(tmp_overlap_ranges)
        for pileup_value in pileup.values():
            if pileup_value>=4:
                site_info["contig_coverage4plus"]+=1
            else:
                site_info["contig_coverage%s"%pileup_value]+=1
        site_info["overlap_length"]=get_length_of_list_range(overlap_ranges)
        
        length=get_length_of_list_range(merged_contig_span_list)
        
        site_info["alignment_length"]=length
    return all_sites
        
def assess_read2_contigs(contig_file, genome_file, all_sites_file, output_sites):
    all_sites=read_site_info_file(all_sites_file)
    all_sites=get_basic_stats(contig_file, all_sites)
    all_sites=get_assessment(contig_file, all_sites, genome_file)
    all_sites=summarize(all_sites)
    if all_sites:
        output_all_sites(all_sites,output_sites)

        
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
    if not options.print_command:
        command_runner.set_command_to_run_localy()
        
    assess_read2_contigs(options.contig_file, options.genome_file,options.sites_file, options.output_sites)
    
    
    

def _prepare_optparser():
    """Prepare optparser object. New options will be added in this
    function first.
    """
    usage = """usage: %prog <-b bam_file> [ -o output_file]"""
    description = """This script will take aligned RAD read to the consensuses and calculate per consensus coverage."""
    
    optparser = OptionParser(version="None",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="show this help message and exit.")
    optparser.add_option("-c","--contig_file",dest="contig_file",type="string",
                         help="Path to one contig_file. Default: %default")
    optparser.add_option("-g","--genome_file",dest="genome_file",type="string",
                         help="Path to one genome_file. Default: %default")
    optparser.add_option("-s","--sites_file",dest="sites_file",type="string",
                         help="Path to the file containing the site information. Default: %default")
    optparser.add_option("-o","--output_sites",dest="output_sites",type="string",
                         help="Path to the output file . Default: %default")
    optparser.add_option("--print",dest="print_command",action='store_true',default=False,
                         help="print the commands instead of running them. Default: %default")
    return optparser


def _verifyOption(options):
    """Check if the mandatory option are present in the options objects.
    @return False if any argument is wrong."""
    arg_pass=True
    
    return arg_pass



if __name__=="__main__":
    main()
    
