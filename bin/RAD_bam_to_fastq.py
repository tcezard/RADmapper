#!/usr/bin/env python
import sys
import utils
from IO_interface import samIterator
from utils.GenomeLoader import GenomeLoader
from utils import DNA_tools, utils_logging
import os, pysam

import logging, threading
from optparse import OptionParser
from overlap import merge_ranges
import random
from collections import OrderedDict
import command_runner
from utils.utils_commands import get_output_stream_from_command
from multiprocessing import Pool, Value, Lock, Manager
import time
import RAD_assemble_read2, RAD_merge_results, RAD_smalt_align_reads_back_to_consensus

def read_white_list(white_list_file):
    all_consensus=[]
    with open(white_list_file) as open_file:
        for line in open_file:
            sp_line = line.strip().split()
            if sp_line and len(sp_line)>0:
                all_consensus.append(sp_line[0])
    return all_consensus

class FuncThread(threading.Thread):
    def __init__(self, target, *args):
        self._target = target
        self._args = args
        threading.Thread.__init__(self)
        self.returned_value=None
 
    def run(self):
        self.returned_value = self._target(*self._args)
        
    def get_returned_value(self):
        return self.returned_value

_all_open_bam={}

def _get_open_bam_file(bam_file):
    if _all_open_bam.has_key(bam_file):
        return _all_open_bam.get(bam_file)
    
    logging.info('open %s'%bam_file)
    open_bam = pysam.Samfile( bam_file, "rb" )
    _all_open_bam[bam_file]=open_bam
    return open_bam
        

def get_pysam_iterator(bam_file, options, consensus):
    sam_file = _get_open_bam_file(bam_file)
    return sam_file.fetch(consensus)


def get_read_group_from_bam_files(bam_files):
    all_read_groups=[]
    for bam_file in bam_files:
        open_bam_file = _get_open_bam_file(bam_file)
        read_groups = open_bam_file.header.get('RG')
        out = ['@RG']
        for read_group in read_groups:
            out.append('ID:%s'%read_group.pop('ID'))
            for key in read_group.keys():
                out.append('%s:%s'%(key,read_group.get(key)))
        all_read_groups.append('\t'.join(out))
    return all_read_groups


def load_from_sites_generator(bam_file, options='', consensus=''):
    """This function return a generator that iterates over read pair where both pairs are mapping.
    @return a tuple containing read1 read2 and the restriction site position."""
    stream, process = utils.get_sam_stream(bam_file, options=options, chomosome_and_position=consensus)
    all_unmatched_read1={}
    all_unmatched_read2={}
    count_line=0
    for line in stream:
        count_line+=1
        sam_record = samIterator.Sam_record(line)
        if sam_record.is_first_read():
            sam_record_r2 = all_unmatched_read2.pop(sam_record.get_query_name(),None)
            if sam_record_r2:
                yield((sam_record,sam_record_r2))
            else:
                all_unmatched_read1[sam_record.get_query_name()]=sam_record
        else:
            sam_record_r1 = all_unmatched_read1.pop(sam_record.get_query_name(),None)
            if sam_record_r1:
                yield((sam_record_r1, sam_record))
            else:
                all_unmatched_read2[sam_record.get_query_name()]=sam_record
    #return_code = process.wait()


def load_from_sites_generator2(bam_file, options='', consensus=''):
    """This function return a generator that iterates over read pair where both pairs are mapping.
    @return a tuple containing read1 read2 and the restriction site position."""
    iter_sam = get_pysam_iterator(bam_file, options=options, consensus=consensus)
    all_unmatched_read1={}
    all_unmatched_read2={}
    count_line=0
    for align_read in iter_sam:
        count_line+=1
        if align_read.is_read1:
            align_read_r2 = all_unmatched_read2.pop(align_read.qname,None)
            if align_read_r2:
                yield((align_read, align_read_r2))
            else:
                all_unmatched_read1[align_read.qname] = align_read
        else:
            sam_record_r1 = all_unmatched_read1.pop(align_read.qname,None)
            if sam_record_r1:
                yield((sam_record_r1, align_read))
            else:
                all_unmatched_read2[align_read.qname]=align_read
    
class File_factory():
    def __init__(self, max_nb_file=50):
        self.all_dir_paths={}
        self.in_dir_counter=0
        self.dir_counter=0
        self.max_nb_file=max_nb_file
        
    def get_dir_name_from_site(self, reference):
        dir_path = self.all_dir_paths.get('%s'%(reference))
        if not dir_path:
            self.in_dir_counter+=1
            if self.in_dir_counter>self.max_nb_file or self.dir_counter==0:
                self.in_dir_counter=1
                self.dir_counter+=1
                if not os.path.exists('%s_dir'%self.dir_counter):
                    os.mkdir('%s_dir'%self.dir_counter)
                    logging.debug("create %s_dir"%self.dir_counter)
            dir_path='%s_dir/%s_dir'%(self.dir_counter,reference)
            self.all_dir_paths['%s'%(reference)]=dir_path
            logging.debug("create %s -- %s"%(dir_path,self.dir_counter))
            print "create %s -- %s"%(dir_path,self.dir_counter)
        return dir_path
    
    def close(self):
        self.list_file_open=[]
        
    def __del__(self):
        self.close()


def get_dir_name_from_site(reference, in_dir_counter, dir_counter, lock, nb_consensus_per_dir=50):
    lock.acquire()
    in_dir_counter.value+=1
    if in_dir_counter.value>nb_consensus_per_dir or dir_counter.value==0:
        in_dir_counter.value=1
        dir_counter.value+=1
        if not os.path.exists('%s_dir'%dir_counter.value):
            os.mkdir('%s_dir'%dir_counter.value)
            logging.debug("create %s_dir"%dir_counter.value)
    dir_path='%s_dir/%s_dir'%(dir_counter.value,reference)
    lock.release()
    return dir_path


def process_one_bam_file_one_consensus(bam_file, consensus_name):
    sam_pair_generator = load_from_sites_generator2(bam_file, consensus=consensus_name)
    all_first_reads_for_consensus=[]
    all_second_reads_for_consensus=[]
    for aligned_read_r1,aligned_read_r2, in sam_pair_generator:
        rgid=aligned_read_r1.opt('RG')
        all_first_reads_for_consensus.append(aligned_read_to_fastq(aligned_read_r1,rgid=rgid))
        all_second_reads_for_consensus.append(aligned_read_to_fastq(aligned_read_r2,rgid=rgid))
        aligned_read_r1=None
        aligned_read_r2=None
    return all_first_reads_for_consensus, all_second_reads_for_consensus


def extract_reads_from_one_consensus(bam_files, output_dir, consensus_name, consensus_sequence):
    all_read1_for_that_consensus=[]
    all_read2_for_that_consensus=[]
    for bam_file in bam_files:
        all_first_reads_for_consensus, all_second_reads_for_consensus = process_one_bam_file_one_consensus(bam_file, consensus_name)
        all_read1_for_that_consensus.extend(all_first_reads_for_consensus)
        all_read2_for_that_consensus.extend(all_second_reads_for_consensus)
    
    consensus_directory = os.path.join(output_dir, consensus_name+'_dir')
    if not os.path.exists(consensus_directory):
        os.mkdir(consensus_directory)
    read1_file=os.path.join(consensus_directory,consensus_name+"_1.fastq")
    read2_file=os.path.join(consensus_directory,consensus_name+"_2.fastq")
    read1_consensus=os.path.join(consensus_directory,consensus_name+"_1.fa")
    open_file = open(read1_consensus,'w')
    open_file.write('>%s\n%s\n'%(consensus_name,consensus_sequence))
    open_file.close()
    
    open_file = open(read1_file,'w')
    open_file.write('\n'.join(all_read1_for_that_consensus))
    open_file.close()
    
    open_file = open(read2_file,'w')
    open_file.write('\n'.join(all_read2_for_that_consensus))
    open_file.close()
    
    return read1_consensus, read1_file, read2_file

def process_one_chomosome(in_dir_counter, dir_counter, lock,
                          bam_files, consensus_name, consensus_sequence, min_nb_read=20, 
                          nb_consensus_per_dir=50):
    print "process %s "%consensus_name
    all_read1_for_that_consensus=[]
    all_read2_for_that_consensus=[]
    for bam_file in bam_files:
        all_first_reads_for_consensus, all_second_reads_for_consensus = process_one_bam_file_one_consensus(bam_file, consensus_name)
        print "process %s --> %s reads from %s"%(consensus_name, len(all_first_reads_for_consensus), bam_file)
        all_read1_for_that_consensus.extend(all_first_reads_for_consensus)
        all_read2_for_that_consensus.extend(all_second_reads_for_consensus)
    
    if len(all_read1_for_that_consensus) > min_nb_read:
        output_directory = get_dir_name_from_site(consensus_name, in_dir_counter, dir_counter, lock, nb_consensus_per_dir)
        if not os.path.exists(output_directory):
            os.mkdir(output_directory)
        read1_file=os.path.join(output_directory,consensus_name+"_1.fastq")
        read2_file=os.path.join(output_directory,consensus_name+"_2.fastq")
        read1_consensus=os.path.join(output_directory,consensus_name+"_1.fa")
        open_file = open(read1_consensus,'w')
        open_file.write('>%s\n%s\n'%(consensus_name,consensus_sequence))
        open_file.close()
        
        open_file = open(read1_file,'w')
        open_file.write('\n'.join(all_read1_for_that_consensus))
        open_file.close()
        
        open_file = open(read2_file,'w')
        open_file.write('\n'.join(all_read2_for_that_consensus))
        open_file.close()

def get_all_consensus_names(bam_file):
    all_consensus_names=[]
    command = "samtools view -H %s | grep 'SQ'"%(bam_file)
    
    stream, process = get_output_stream_from_command(command)
    
    for line in stream:
        sp_line = line.split("\t")
        if len(sp_line)>1 and sp_line[0]=="@SQ":
            for element in sp_line[1:]:
                tag, value = element.split(":")
                if tag=="SN":
                    all_consensus_names.append(value)
                    break
    return all_consensus_names



def process_all_bam_files_with_pool(bam_files, all_read1_consensus_file, min_nb_read, nb_consensus_per_dir=50, number_of_process=10):
    genome_loader = GenomeLoader(all_read1_consensus_file)
    pool = Pool(processes=number_of_process)
    manager = Manager()
    lock = manager.Lock()
    in_dir_counter = manager.Value('i', 0)
    dir_counter = manager.Value('i', 0)
    #consensus_name, consensus_sequence = genome_loader.next()
    #pool.apply(process_one_chomosome, args=[in_dir_counter, dir_counter, lock, bam_files, consensus_name, consensus_sequence, min_nb_read])
    for consensus_name, consensus_sequence in genome_loader:
        pool.apply_async(process_one_chomosome, args=[in_dir_counter, dir_counter, lock, bam_files, consensus_name, consensus_sequence, min_nb_read, nb_consensus_per_dir])
    pool.close()
    pool.join() 

def extract_reads_from_all_bam_files_all_consensus(bam_files, all_read1_consensus_file, white_list_file,
                                                   nb_consensus_per_dir=50, number_of_process=10):
    list_consensus = read_white_list(white_list_file)
    list_of_list_consensus=[]
    temp_list = [] 
    i=0
    for consensus in list_consensus:
        i+=1
        temp_list.append(consensus)
        if i%nb_consensus_per_dir == 0:
            list_of_list_consensus.append(temp_list)
            temp_list = []
    if len(temp_list)>0:
        list_of_list_consensus.append(temp_list)
    print "split into %s jobs of %s consensus"%(len(list_of_list_consensus),nb_consensus_per_dir)
    
    genome_loader = GenomeLoader(all_read1_consensus_file)
    list_of_list_consensus=list_of_list_consensus[:10]
    
    #all_read_groups = get_read_group_from_bam_files(bam_files)
    
    #Generate the assembly function list from the assembler to try
    #all_assembler_to_try=['velvet']
    #assembly_function_list=[]
    #for assembler_name in all_assembler_to_try:
    #    assembly_function=RAD_assemble_read2.get_assembly_function(assembler_name)
    #    assembly_function_list.append(assembly_function)
    
    
    for i, list_consensus in enumerate(list_of_list_consensus):
        output_dir = os.path.join(os.path.curdir, '%s_dir'%(i+1))
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
        extract_reads_from_all_bam_files_set_of_consensus(bam_files, genome_loader, list_consensus, output_dir)
        

def split_whitelist(white_list_file,nb_consensus_per_dir):
    list_consensus = read_white_list(white_list_file)
    list_of_list_consensus=[]
    temp_list = [] 
    i=0
    for consensus in list_consensus:
        i+=1
        temp_list.append(consensus)
        if i%nb_consensus_per_dir == 0:
            list_of_list_consensus.append(temp_list)
            temp_list = []
    if len(temp_list)>0:
        list_of_list_consensus.append(temp_list)
    #print "split into %s jobs of %s consensus"%(len(list_of_list_consensus),nb_consensus_per_dir)
    for i, list_of_consensus in enumerate(list_of_list_consensus):
        output_dir = os.path.join(os.path.curdir, '%s_dir'%(i+1))
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
        output_file=os.path.join(output_dir, 'whitelist.txt')
        with open(output_file, 'w') as open_file:
            open_file.write('\n'.join(list_of_consensus))
        yield (output_dir, output_file)
        


def extract_reads_from_all_bam_files_set_of_consensus(bam_files, list_consensus, output_dir, genome_loader=None, all_read1_consensus_file=None):
    if genome_loader is None:
        genome_loader = GenomeLoader(all_read1_consensus_file)
    for consensus_name in list_consensus:
        logging.info("Extract reads from %s "%consensus_name)
        consensus_name, consensus_sequence = genome_loader.get_chr(consensus_name)
        extract_reads_from_one_consensus(bam_files, output_dir, consensus_name, consensus_sequence)
    
def sam_record_to_fastq(sam_record, rgid=None):
    out=[]
    if sam_record.is_first_read():
        suffix='/1'
    elif sam_record.is_second_read():
        suffix='/2'
    else:
        suffix=''
    if rgid:
        rgid_str='RGID:%s'%(rgid)
    else:
        rgid_str=''
    out.append('@%s%s%s'%(sam_record.get_query_name(), rgid_str, suffix))
    if sam_record.is_second_read():
        out.append(DNA_tools.rev_complements(sam_record.get_query_sequence()))
        out.append('+')
        out.append(sam_record.get_query_quality()[::-1])
    else:
        out.append(sam_record.get_query_sequence())
        out.append('+')
        out.append(sam_record.get_query_quality())
    return '\n'.join(out)

def aligned_read_to_fastq(aligned_read, rgid=None):
    out=[]
    if aligned_read.is_read1:
        suffix='/1'
    elif aligned_read.is_read2:
        suffix='/2'
    else:
        suffix=''
    if rgid:
        rgid_str='RGID:%s'%(rgid)
    else:
        rgid_str=''
    out.append('@%s%s%s'%(aligned_read.qname, rgid_str, suffix))
    if aligned_read.is_read2:
        out.append(DNA_tools.rev_complements(aligned_read.seq))
        out.append('+')
        out.append(aligned_read.qual[::-1])
    else:
        out.append(aligned_read.seq)
        out.append('+')
        out.append(aligned_read.qual)
    return '\n'.join(out)


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
    if options.debug:
        utils_logging.init_logging(logging.DEBUG)
    bam_files=[options.bam_file]
    if len(args)>0:
        bam_files.extend(args)
    command_runner.set_command_to_run_localy()
    if options.output_dir:
        list_consensus = read_white_list(options.white_list_file)
        extract_reads_from_all_bam_files_set_of_consensus(bam_files, list_consensus, options.output_dir,
                                                          genome_loader=None, all_read1_consensus_file=options.read1_consensus_file)
    else:
        iter_dir_and_file = split_whitelist(white_list_file=options.white_list_file, nb_consensus_per_dir=options.nb_consensus_per_dir)
        for output_dir, sub_whitelist in iter_dir_and_file:
            command = 'python %s -g %s -w %s -o %s -b %s'%(sys.argv[0], options.read1_consensus_file, sub_whitelist, output_dir, ' '.join(bam_files))
            print command
            


def _prepare_optparser():
    """Prepare optparser object. New options will be added in this
    function first.
    """
    usage = """usage: %prog <-b bam_file> <-g genome_file> <-k known_sites>"""
    description = """This script extract reads from an aligned bam file and create the corresponding fastq files."""
    
    optparser = OptionParser(description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="show this help message and exit.")
    optparser.add_option("-b","--bam_file",dest="bam_file",type="string",
                         help="The bam file from which the reads should be extracted.")
    optparser.add_option("-g","--read1_consensus_file",dest="read1_consensus_file",type="string",
                         help="The fasta file containing the read1 consensus file corresponding to the bam file.")
    optparser.add_option("-w","--white_list_file",dest="white_list_file",type="string",
                         help="The file containing the name of the consensus to assemble")
    optparser.add_option("-o","--output_dir",dest="output_dir",type="string",
                         help="This act as a flag to bypass the nb_consensus_per_dir. This means that all consensus in the whitelist will be extracted in the output_dir.")
    optparser.add_option("-p","--number_of_process",dest="number_of_process",type="int",default=10,
                         help="The number of process used to extract the information from the bam files.")
    optparser.add_option("-n","--nb_consensus_per_dir",dest="nb_consensus_per_dir",type="int",default=50,
                         help="The number of that should be held in each directory.")
    optparser.add_option("--debug",dest="debug",action='store_true',default=False,
                         help="Output debug statment. Default: %default")
    
    return optparser


def _verifyOption(options):
    """Check if the mandatory option are present in the options objects.
    @return False if any argument is wrong."""
    arg_pass=True
    
    if not options.bam_file:
        logging.error("You must specify a bam file -b.")
        arg_pass=False
    if not options.read1_consensus_file:
        logging.error("You must specify a fasta file with the read1 consensus with -g.")
        arg_pass=False
    elif not os.path.exists(options.read1_consensus_file):
        logging.error("You must specify a an existing fasta file with -g.")
        arg_pass=False
    if not options.white_list_file:
        logging.error("You must specify a list of consensus to use with -w.")
        arg_pass=False
    elif not os.path.exists(options.white_list_file):
        logging.error("You must specify a an existing list of consensus to use with -w.")
        arg_pass=False
    return arg_pass



if __name__=="__main__":
    main()        


if __name__=="1__main__":
    bam_files=sys.argv[1:]
    for bam_file in bam_files:
        open_bam = _get_open_bam_file(bam_file)
    