import sys
import utils
from IO_interface import samIterator
from utils.GenomeLoader import GenomeLoader
from utils import DNA_tools, utils_logging
import os
import logging, threading
from optparse import OptionParser
from overlap import merge_ranges
import random
from collections import OrderedDict
import command_runner
from utils.utils_commands import get_output_stream_from_command
from multiprocessing import Pool, Value, Lock, Manager

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
 
        
def load_from_sites_generator(bam_file, options='', consensus=''):
    """This function return a generator that iterates over read pair where both pairs are mapping.
    @return a tuple containing read1 read2 and the restriction site position."""
    stream = utils.get_sam_stream(bam_file, options=options, chomosome_and_position=consensus)
    all_unmatched_read1={}
    all_unmatched_read2={}
    count_line=0
    for line in stream:
        count_line+=1
        #if count_line%10000==0:
        #    print count_line, len(all_unmatched_read1), len(all_unmatched_read2)
        sam_record = samIterator.Sam_record(line)
        if sam_record.is_first_read():
            sam_record_r2 = all_unmatched_read2.pop(sam_record.get_query_name(),None)
            info=None
            if not sam_record.is_unmapped():
                info = get_RAD_site_from_sam_record(sam_record)
            if sam_record_r2:
                if info:
                    yield((sam_record,sam_record_r2, info))
            else:
                all_unmatched_read1[sam_record.get_query_name()]=sam_record
        else:
            sam_record_r1 = all_unmatched_read1.pop(sam_record.get_query_name(),None)
            if sam_record_r1:
                info = get_RAD_site_from_sam_record(sam_record_r1)
                if info:
                    yield((sam_record_r1, sam_record, info))
            else:
                all_unmatched_read2[sam_record.get_query_name()]=sam_record


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


def get_dir_name_from_site(reference, in_dir_counter, dir_counter, lock):
    lock.acquire()
    in_dir_counter.value+=1
    if in_dir_counter.value>50 or dir_counter.value==0:
        in_dir_counter.value=1
        dir_counter.value+=1
        if not os.path.exists('%s_dir'%dir_counter.value):
            os.mkdir('%s_dir'%dir_counter.value)
            logging.debug("create %s_dir"%dir_counter.value)
    dir_path='%s_dir/%s_dir'%(dir_counter.value,reference)
    lock.release()
    return dir_path


def process_one_bam_file_one_consensus(bam_file, consensus_name):
    sam_pair_generator = load_from_sites_generator(bam_file, consensus=consensus_name)
    all_first_reads_for_consensus=[]
    all_second_reads_for_consensus=[]
    for sam_record_r1,sam_record_r2, site_info in sam_pair_generator:
        rgid=sam_record_r1.get_tag('RG')
        all_first_reads_for_consensus.append(sam_record_to_fastq(sam_record_r1,rgid=rgid))
        all_second_reads_for_consensus.append(sam_record_to_fastq(sam_record_r2,rgid=rgid))
        
    return all_first_reads_for_consensus, all_second_reads_for_consensus



def process_one_chomosome(in_dir_counter, dir_counter, lock,
                          bam_files, consensus_name, consensus_sequence, min_nb_read=20):
    print "process %s "%consensus_name
    all_read1_for_that_consensus=[]
    all_read2_for_that_consensus=[]
    for bam_file in bam_files:
        all_first_reads_for_consensus, all_second_reads_for_consensus = process_one_bam_file_one_consensus(bam_file, consensus_name)
        all_read1_for_that_consensus.extend(all_first_reads_for_consensus)
        all_read2_for_that_consensus.extend(all_second_reads_for_consensus)
    
    if len(all_read1_for_that_consensus) > min_nb_read:
        output_directory = get_dir_name_from_site(consensus_name, in_dir_counter, dir_counter, lock)
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



def process_all_bam_files_with_pool(bam_files, all_read1_consensus_file, min_nb_read):
    genome_loader = GenomeLoader(all_read1_consensus_file)
    pool = Pool(processes=10)
    manager = Manager()
    lock = manager.Lock()
    in_dir_counter = manager.Value('i', 0)
    dir_counter = manager.Value('i', 0)
    #consensus_name, consensus_sequence = genome_loader.next()
    #pool.apply(process_one_chomosome, args=[in_dir_counter, dir_counter, lock, bam_files, consensus_name, consensus_sequence, min_nb_read])
    for consensus_name, consensus_sequence in genome_loader:
        pool.apply_async(process_one_chomosome, args=[in_dir_counter, dir_counter, lock, bam_files, consensus_name, consensus_sequence, min_nb_read])
    pool.close()
    pool.join() 


def get_RAD_site_from_sam_record(sam_record):
    if sam_record.is_first_read() and not sam_record.is_unmapped():
        ref=sam_record.get_reference_name()
        pos=int(sam_record.get_alignment_start())
        return (ref, pos)
    return None    
            
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
    bam_files=[options.bam_file]
    if len(args)>0:
        bam_files.extend(args)
    
    #process_all_bam_files(bam_files=bam_files, 
    #                      min_nb_read=options.min_coverage)
    process_all_bam_files_with_pool(bam_files=bam_files,
                                    all_read1_consensus_file=options.read1_consensus_file,
                                    min_nb_read=options.min_coverage)
    
    


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
                         help="The fasta file contanins the read1 consensus file corresponding to the bam file.")
    optparser.add_option("-c","--min_coverage",dest="min_coverage",type="int",default=20,
                         help="The minimum coverage to take create a file with read2 reads.")
    return optparser


def _verifyOption(options):
    """Check if the mandatory option are present in the options objects.
    @return False if any argument is wrong."""
    arg_pass=True
    
    if not options.bam_file:
        logging.error("You must specify a bam file -b.")
        arg_pass=False
    return arg_pass



if __name__=="__main__":
    main()        

