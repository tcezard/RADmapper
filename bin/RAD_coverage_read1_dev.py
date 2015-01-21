#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 20 October 2011
@author: tcezard
'''
import sys
import os
import logging
import threading
from optparse import OptionParser
from collections import Counter, defaultdict

from utils import utils_param, utils_logging
from utils.parameters import Config_file_error
from utils.utils_commands import get_output_stream_from_command
from IO_interface.samIterator import Sam_record


class FuncThread(threading.Thread):
    def __init__(self, target, *args, **kwargs):
        self._target = target
        self._args = args
        self.kwargs=kwargs
        threading.Thread.__init__(self)
        self.returned_value=None
 
    def run(self):
        self.returned_value = self._target(*self._args, **self.kwargs)
        
    def get_returned_value(self):
        return self.returned_value
 

class AllContigsInfo():
    def __init__(self):
        self.per_loci_info={}
        self.sample_list=set()
        self.lock = threading.RLock()

    def add_values(self, contig, position, coverage, duplicate, alleles=None, sample=None):
        with self.lock:
            if position:
                loci = '%s:%s' % (contig, position)
            else:
                loci = '%s' % (contig)
            info = self.per_loci_info.get(loci)
            if not info:
                info = {'coverage': 0, 'coverage_mrk_dup': 0, "alleles": Counter()}
                self.per_loci_info[loci]=info
            info['coverage']+=coverage
            info['coverage_mrk_dup']+=coverage-duplicate
            if alleles:
                for allele in alleles:
                    info['alleles'][allele] += alleles.get(allele)
            if sample:
                self.sample_list.add(sample)
                if not info.has_key(sample):
                    info[sample] = {'coverage': 0, 'coverage_mrk_dup': 0, 'alleles': Counter()}
                info[sample]['coverage'] += coverage
                info[sample]['coverage_mrk_dup'] += coverage - duplicate
                if alleles:
                    for allele in alleles:
                        info[sample]['alleles'][allele] += alleles.get(allele)


    def get_all_loci(self):
        return self.per_loci_info.keys()
    
    def get_all_samples(self):
        return list(self.sample_list)


    def get_loci_coverage_value(self, loci, sample=None):
        # loci='%s\t%s'%(contig, position)
        info = self.per_loci_info.get(loci)
        if info:
            if sample:
                tmp = info.get(sample)
                if tmp: 
                    return tmp.get('coverage')
            else:
                return info.get('coverage')

    def get_locus_coverage_mrk_dup_value(self, loci, sample=None):
        #loci='%s\t%s'%(contig, position)
        info = self.per_loci_info.get(loci)
        if info:
            if sample:
                tmp = info.get(sample)
                if tmp: 
                    return tmp.get('coverage_mrk_dup')
            else:
                return info.get('coverage_mrk_dup')

    def get_nb_alleles(self, loci, sample=None):
        info = self.per_loci_info.get(loci)
        alleles = None
        if info:
            if sample:
                tmp = info.get(sample)
                if tmp:
                    alleles = tmp.get('alleles')
            else:
                alleles = info.get('alleles')
        return detect_alleles(alleles)

    def is_known_length(self, locus):
        return False


class AllSitesInfo():
    def __init__(self, known_sites=None):
        self.per_loci_info={}
        self.sample_list=set()
        self.lock = threading.RLock()
        if known_sites:
            for site in known_sites:
                self.add_values(site, coverage=0, duplicate=0, known_length=known_sites.get(site))

    def add_values(self, loci, coverage, duplicate, alleles=None, sample=None, known_length=False):
        with self.lock:
            info = self.per_loci_info.get(loci)
            if not info:
                info = {'coverage': 0, 'coverage_mrk_dup': 0, 'alleles': {}}
                self.per_loci_info[loci]=info
            info['coverage']+=coverage
            info['coverage_mrk_dup']+=coverage-duplicate
            if alleles:
                for allele in alleles:
                    info['alleles'] += allele
            if known_length:
                info['known_length']=known_length
            if sample:
                self.sample_list.add(sample)
                if not info.has_key(sample):
                    info[sample] = {'coverage': 0, 'coverage_mrk_dup': 0, 'alleles': {}}
                info[sample]['coverage'] += coverage
                info[sample]['coverage_mrk_dup'] += coverage - duplicate
                if alleles:
                    for allele in alleles:
                        info['alleles'] += allele


    def get_all_loci(self):
        return self.per_loci_info.keys()

    def get_all_samples(self):
        return list(self.sample_list)


    def get_loci_coverage_value(self, loci, sample=None):
        info = self.per_loci_info.get(loci)
        if info:
            if sample:
                tmp = info.get(sample)
                if tmp:
                    return tmp.get('coverage')
            else:
                return info.get('coverage')

    def get_locus_coverage_mrk_dup_value(self, loci, sample=None):
        info = self.per_loci_info.get(loci)
        if info:
            if sample:
                tmp = info.get(sample)
                if tmp:
                    return tmp.get('coverage_mrk_dup')
            else:
                return info.get('coverage_mrk_dup')

    def get_nb_alleles(self, loci, sample=None):
        info = self.per_loci_info.get(loci)
        alleles = None
        if info:
            if sample:
                tmp = info.get(sample)
                if tmp:
                    alleles = tmp.get('alleles')
            else:
                alleles = info.get('alleles')
        return detect_alleles(alleles)


    def is_known_length(self, loci):
        info = self.per_loci_info.get(loci)
        return info.get("known_length", None)


def detect_alleles(alleles):
    nb_allele = 0
    if alleles:
        for allele in alleles:
            if alleles.get(allele) > 10:
                nb_allele += 1
    return nb_allele


def process_single_samtools_run(bam_file, all_contigs_info, samtools_bin):
    command="%s view -F 132 %s"%(samtools_bin, bam_file)
    open_stream, process=get_output_stream_from_command(command)
    current_contig=None
    coverage=0
    duplicate=0
    sample_name, ext = os.path.splitext(bam_file)
    for line in open_stream:
        sp_line=line.strip().split()
        if current_contig!=sp_line[2] and current_contig != None:
            all_contigs_info.add_values(current_contig, coverage, duplicate, sample=sample_name)
            coverage=0
            duplicate=0
        current_contig=sp_line[2]
        if int(sp_line[3])==1:
            if int(sp_line[1]) & 1024 == 1024:
                duplicate+=1
            coverage+=1
            
    open_stream.close()


def process_single_samtools_run_with_read_group(bam_file,all_contigs_info,samtools_bin):
    command="%s view -h -F 132 %s"%(samtools_bin, bam_file)
    open_stream, process = get_output_stream_from_command(command)
    current_contig=None
    sample_name, ext = os.path.splitext(bam_file)
    read_groups={}
    try:
        for line in open_stream:
            if not line.startswith("@"):
                break
            if line.startswith("@RG"):
                sp_line = line.strip().split()
                rg_id=rg_sample=rg_library=None
                for value in sp_line:
                    if value.startswith("ID"):
                        rg_id=value[3:]
                    elif value.startswith("SM"):
                        rg_sample=value[3:]
                    elif value.startswith("LB"):
                        rg_library=value[3:]
                if rg_id:
                    if rg_sample:
                        read_groups[rg_id]=rg_sample
                    elif rg_library:
                        read_groups[rg_id]=rg_library
                    else:
                        read_groups[rg_id]=rg_id
        all_sample_coverage={}
        all_sample_coverage_reads = {}
        all_sample_duplicate={}
        for sample in read_groups.values():
            all_sample_coverage[sample]=Counter()
            all_sample_duplicate[sample]=Counter()
            all_sample_coverage_reads[sample] = defaultdict(Counter)
        #process the first read
        sam_record = Sam_record(line.strip())
        current_contig = sam_record.get_reference_name()
        if not sam_record.is_unmapped():
            rg_id = sam_record.get_tag("RG")
            read_sequence = sam_record.get_query_sequence()
            loci = get_loci_from_read(sam_record)
            if sam_record.is_duplicate_read():
                all_sample_duplicate[read_groups.get(rg_id)][str(loci)]+=1
            all_sample_coverage[read_groups.get(rg_id)][str(loci)]+=1
            all_sample_coverage_reads[read_groups.get(rg_id)][str(loci)][read_sequence] +=1
        i=1
        #process all the others
        for line in open_stream:
            i+=1
            if i%1000000==0:
                print i
            sam_record = Sam_record(line.strip())
            if current_contig != sam_record.get_reference_name() and current_contig != None:
                for sample in read_groups.values():
                    for loci in all_sample_coverage.get(sample):
                        alleles = all_sample_coverage_reads[sample].get(loci)
                        all_contigs_info.add_values(current_contig, loci, all_sample_coverage.get(sample).get(loci, 0),
                                                    all_sample_duplicate.get(sample).get(loci, 0), alleles=alleles,
                                                    sample=sample)

                    all_sample_coverage[sample]=Counter()
                    all_sample_duplicate[sample]=Counter()
                    all_sample_coverage_reads[sample] = defaultdict(Counter)
            current_contig = sam_record.get_reference_name()
            
            if not sam_record.is_unmapped():
                rg_id = sam_record.get_tag("RG")
                loci = get_loci_from_read(sam_record)
                read_sequence = sam_record.get_query_sequence()
                if sam_record.is_duplicate_read():
                    all_sample_duplicate[read_groups.get(rg_id)][str(loci)]+=1
                all_sample_coverage[read_groups.get(rg_id)][str(loci)]+=1
                all_sample_coverage_reads[read_groups.get(rg_id)][str(loci)][read_sequence] +=1
        if current_contig != None:
            for sample in read_groups.values():
                for loci in all_sample_coverage.get(sample):
                    alleles = all_sample_coverage_reads[sample].get(loci)
                    all_contigs_info.add_values(current_contig, loci, all_sample_coverage.get(sample).get(loci, 0),
                                                all_sample_duplicate.get(sample).get(loci, 0), alleles=alleles,
                                                sample=sample)
                all_sample_coverage[sample]=Counter()
                all_sample_duplicate[sample]=Counter()
                all_sample_coverage_reads[sample] = defaultdict(Counter)
    finally:
        open_stream.close()

def load_from_sites_generator(stream):
    all_unmatched_read1={}
    all_unmatched_read2={}
    count_line=0
    for line in stream:
        count_line+=1
        if count_line%10000==0:
            sys.stderr.write('%s %s %s\n'%(count_line, len(all_unmatched_read1), len(all_unmatched_read2)))
        sam_record = Sam_record(line)
        if sam_record.is_first_read():
            sam_record_r1 = sam_record
            sam_record_r2 = all_unmatched_read2.pop(sam_record.get_query_name(),None)
            if not sam_record_r2:
               all_unmatched_read1[sam_record.get_query_name()]=sam_record
        else:
            sam_record_r2 = sam_record
            sam_record_r1 = all_unmatched_read1.pop(sam_record.get_query_name(),None)
            if not sam_record_r1:
                all_unmatched_read2[sam_record.get_query_name()]=sam_record

        if sam_record_r1 and sam_record_r2:
            yield  ((sam_record_r1,sam_record_r2))

def load_known_ddRAD_restriction_site_file(restriction_file):
    """Load known ddRAD sites from a file
    Excpected format is:
    chr1  885347  885975  +  631
    chr1  892323  892573  +  253"""
    all_known_sites={}
    with open(restriction_file) as open_file:
        for line in open_file:
            sp_line = line.strip().split()
            all_known_sites["_".join(sp_line[0:4])]=sp_line[4]
    return all_known_sites


def get_restriction_site_from_sam_record(sam_record):
    """Return the position of the alignment and orientation."""
    if not sam_record.is_unmapped():
        ref=sam_record.get_reference_name()
        pos=int(sam_record.get_alignment_start())
        align_length=sam_record.get_alignment_length()
        sequence = sam_record.get_query_sequence()
        if sam_record.is_reverse_strand():
            site_position = pos + align_length
            strand='-'
            sequence = DNA_tools.rev_complements(sequence)
        else:
            site_position = pos-1
            strand='+'
        return (ref, int(site_position), strand)
    return None

def process_double_digest_rad_run(bam_file,all_sites_info,samtools_bin):
    command="%s view -h %s"%(samtools_bin, bam_file)
    open_stream, process = get_output_stream_from_command(command)
    sample_name, ext = os.path.splitext(bam_file)
    read_groups={}
    try:
        for line in open_stream:
            if not line.startswith("@"):
                break
            if line.startswith("@RG"):
                sp_line = line.strip().split()
                rg_id=rg_sample=rg_library=None
                for value in sp_line:
                    if value.startswith("ID"):
                        rg_id=value[3:]
                    elif value.startswith("SM"):
                        rg_sample=value[3:]
                    elif value.startswith("LB"):
                        rg_library=value[3:]
                if rg_id:
                    if rg_sample:
                        read_groups[rg_id]=rg_sample
                    elif rg_library:
                        read_groups[rg_id]=rg_library
                    else:
                        read_groups[rg_id]=rg_id

        all_sample_coverage={}
        for sample in read_groups.values():
            all_sample_coverage[sample]=Counter()

        i=0
        for sam_record_r1,sam_record_r2 in load_from_sites_generator(open_stream):
            duplicate=0
            i+=1
            if i%1000000==0:
                print i
            if not sam_record_r1.is_unmapped() and not sam_record_r2.is_unmapped():
                loci = get_dd_RAD_loci_from_read_pair(sam_record_r1,sam_record_r2)
                if sam_record_r1.is_duplicate_read():
                    duplicate=1
                all_sites_info.add_values(loci, coverage=1, duplicate=duplicate,
                                          sample=read_groups.get(sam_record_r1.get_tag("RG")))
    finally:
        open_stream.close()


def get_dd_RAD_loci_from_read_pair(sam_record_r1, sam_record_r2):
    #The strand is reverese compare to the known site (I'm not sure why)
    if sam_record_r1.is_reverse_strand():
        position1 = sam_record_r1.get_position() + sam_record_r1.get_alignment_length()
        strand1 = "-"
    else:
        position1 = sam_record_r1.get_position()
        strand1 = "+"

    if sam_record_r2.is_reverse_strand():
        position2 = sam_record_r2.get_position() + sam_record_r2.get_alignment_length()
        strand2 = "-"
    else:
        position2 = sam_record_r2.get_position()
        strand2="+"

    if sam_record_r1.get_reference_name() == sam_record_r2.get_reference_name():
        ref=sam_record_r1.get_reference_name()
    else:
        ref=sam_record_r1.get_reference_name()+"_"+sam_record_r2.get_reference_name()
    if strand1 == '-':
        return  "%s_%s_%s_%s"%(ref,position2,position1,strand1)
    else:
        return  "%s_%s_%s_%s"%(ref,position1,position2,strand1)

def get_loci_from_read(sam_record):
    if sam_record.is_reverse_strand():
        return '%s_%s'%(sam_record.get_position() + sam_record.get_alignment_length() , '-')
    else:
        return '%s_%s'%(sam_record.get_position(), '+')


def RAD_median_coverage(bam_files,output_file, with_rg=False, dd_rad=False, restriction_file=None):
    try:
        pipeline_param=utils_param.get_pipeline_parameters()
        samtools_dir=pipeline_param.get_samtools_dir()
    except Config_file_error, e:
        #logging.exception('Config_file_error:')
        logging.warning("You'll need to have samtools in your path")
        samtools_dir=''
    samtools_bin=os.path.join(samtools_dir,"samtools")
    all_threads=[]
    if dd_rad:
        known_sites=load_known_ddRAD_restriction_site_file(restriction_file)
        all_contigs_info=AllSitesInfo(known_sites=known_sites)
    else:
        all_contigs_info=AllContigsInfo()
    
    if with_rg and dd_rad:
        function = process_double_digest_rad_run
    elif with_rg:
        function = process_single_samtools_run_with_read_group
    else:
        function = process_single_samtools_run

    for bam_file in bam_files:
        sample_name, ext = os.path.splitext(bam_file)
        t = FuncThread(function, bam_file, all_contigs_info, samtools_bin)
        t.start()
        all_threads.append(t)
    for t in all_threads:
        t.join()

    all_loci = all_contigs_info.get_all_loci()
    all_loci.sort()
    sample_names=all_contigs_info.get_all_samples()
    sample_names.sort()
    open_output=utils_logging.open_output_file(output_file)
    if restriction_file:
        open_output.write("#locus\tcoverage\tcoverage_mrk_dup\tnb_allele\tlength\tknown_novel\tnb_sample")
    else:
        open_output.write("#locus\tcoverage\tcoverage_mrk_dup\tnb_allele\tnb_sample")
    for sample in sample_names:
        open_output.write("\t%s\t%s_mrk_dup\t%s_nb_allele" % (sample, sample, sample))
    open_output.write("\n")
    for locus in all_loci:
        out=[locus]
        nb_sample=0
        coverage=all_contigs_info.get_loci_coverage_value(locus)
        coverage_mrk_dup=all_contigs_info.get_locus_coverage_mrk_dup_value(locus)
        nb_allele = all_contigs_info.get_nb_alleles(locus)
        out.append("%s"%(coverage))
        out.append("%s"%(coverage_mrk_dup))
        out.append("%s" % (nb_allele))
        if restriction_file:
            length = all_contigs_info.is_known_length(locus)
            if length:
                out.append("%s"%length)
                out.append("known")
            else:
                out.append("0")
                out.append("novel")

        for sample in sample_names:
            coverage_mrk_dup=all_contigs_info.get_locus_coverage_mrk_dup_value(locus, sample)
            if coverage_mrk_dup and coverage_mrk_dup>2:
                nb_sample+=1
        out.append("%s"%(nb_sample))
        
        for sample in sample_names:
            coverage = all_contigs_info.get_loci_coverage_value(locus, sample)
            coverage_mrk_dup=all_contigs_info.get_locus_coverage_mrk_dup_value(locus, sample)
            nb_allele = all_contigs_info.get_nb_alleles(locus, sample)
            if not coverage:
                coverage=0
            if not coverage_mrk_dup:
                coverage_mrk_dup=0
            if not nb_allele:
                nb_allele = 0
            out.append("%s\t%s\t%s" % (coverage, coverage_mrk_dup, nb_allele))
            
        open_output.write("%s\n"%('\t'.join(out)))
        
    
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
    bam_files=[options.bam_files]
    for file in args:
        if os.path.exists(file):
            bam_files.append(file)
    RAD_median_coverage(bam_files,options.output_file, with_rg=options.with_read_group, dd_rad=options.ddRAD,
                        restriction_file=options.known_sites_file)


def _prepare_optparser():
    """Prepare optparser object. New options will be added in this
    function first.
    """
    usage = """usage: %prog <-b bam_file> [ -o output_file]"""
    description = """This script will take aligned RAD read to an assembled genome and calculate per locus coverage."""
    
    optparser = OptionParser(version="None",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="show this help message and exit.")
    optparser.add_option("-b","--bam_files",dest="bam_files",type="string",
                         help="Path to one or several bam files from which the coverage should be extracted. Default: %default")
    optparser.add_option("-o","--output_file",dest="output_file",type="string",
                         help="The path to the file where the data should be output. If not set, the results will be print to stdout. Default: %default")
    optparser.add_option("-k","--known_sites_file",dest="known_sites_file",type="string",
                         help="The path to the file where That contains the known sites. Default: %default")
    optparser.add_option("-r","--with_read_group",dest="with_read_group",action="store_true",default=False,
                         help="Make the script use the read group of the bam files to determine samples instead of the file name. Default: %default")
    optparser.add_option("-d","--ddRAD",dest="ddRAD",action="store_true",default=False,
                         help="Process the data as double digest RAD assuming a restriction site on both side. Default: %default")
    return optparser


def _verifyOption(options):
    """Check if the mandatory option are present in the options objects.
    @return False if any argument is wrong."""
    arg_pass=True
    
    if not options.bam_files :
        logging.error("You must specify a bam file.")
        arg_pass=False
    return arg_pass



if __name__=="__main__":
    main()
    
