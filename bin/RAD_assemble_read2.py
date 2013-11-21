#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 19 March  2012
@author: tcezard
'''
import sys, os, re
from utils import utils_logging, utils_commands
import logging
from optparse import OptionParser
from glob import glob
import command_runner
from utils.FastaFormat import FastaReader
import time
from RAD_merge_read1_and_read2 import merge_2_contigs

VELVETOPT_ALIAS=["velvetopt","velvetoptimiser"]
VELVET_ALIAS=["velvet"]
VELVETK_ALIAS=["velvetk", "velvetk.pl", "velvetK"]
VELVETKOPT_ALIAS=["velvetkopt", "velvetKopt"]
VELVETKCLC_ALIAS=["velvetkclc"]
CLC_ALIAS=["clc","clcbio"]
ABYSS_ALIAS=["abyss"]
LOCASOPT_ALIAS=["locasopt","locasoptimiser"]
LOCAS_ALIAS=["locas"]
SOAPDENOVO_ALIAS=["soap","soapdenovo"]
IDBA_ALIAS=["idba","idba_ud"]
NO_ASSEMBLER=["None", "no", ""]


def run_velvetOptimiser(fastq_file_name, low_k=19, high_k=99, outputdir='velvetopt', **kwarg):
    command='rm -rf %s'%outputdir
    if os.path.exists(outputdir):
        return_code = command_runner.run_command(command)
    log_file='%s.log'%outputdir
    command_tmp="VelvetOptimiser.pl -f '-fastq -short %s' --s %s --e %s --k max --c max --d %s 2>&1 >%s"
    command=command_tmp%(fastq_file_name,low_k,high_k,outputdir, log_file)
    return_code = command_runner.run_command(command)
    
    contig_files = glob('%s/contigs.fa'%outputdir)
    # If only one contig file exists, as it should if VelvetOptimiser runs
    # successfully, write out the assembled contig(s)
    contig_file_name=None
    if len(contig_files) == 1:
        contig_file_name = contig_files[0]
    return contig_file_name;


def run_velvet(fastq_file_name, kmer_length=29, output_dir= 'velvet', **kwarg):
    log_file='%s.log'%(output_dir)
    command = "velveth %s %s -fastq -short %s 2>&1 >%s"%(output_dir, kmer_length, fastq_file_name, log_file)
    return_code = command_runner.run_command(command)
    command = "velvetg %s 2>&1 >%s"%(output_dir, log_file)
    return_code = command_runner.run_command(command)

    contig_files = glob('%s/contigs.fa'%output_dir)
    contig_file_name=None
    if len(contig_files) == 1:
        contig_file_name = contig_files[0]
    
    return contig_file_name;


def run_velvetk(fastq_file_name, estimated_size=600, **kwarg):
    command = "/home/tcezard/velvetk.pl --size %s --best %s 2> /dev/null"%(estimated_size, fastq_file_name)
    logging.info(command)
    stream,process = utils_commands.get_output_stream_from_command(command)
    kmer_length=29
    for line in stream:
        if line.strip().isdigit():
            kmer_length = int(line.strip())
    if kmer_length<19: kmer_length=19
    elif kmer_length>99: kmer_length=99
    logging.info("velvetk kmer: %s"%kmer_length)
    return kmer_length

def run_velvetk_plus_velvet(fastq_file_name, estimated_size=600, **kwarg):
    kmer_length = run_velvetk(fastq_file_name, estimated_size)
    output_dir='velvetk_plus_velvet'
    return run_velvet(fastq_file_name, kmer_length, output_dir)

def run_velvetk_plus_velvetopt(fastq_file_name, estimated_size=600, **kwarg):
    half_range_size=10
    kmer_length = run_velvetk(fastq_file_name, estimated_size)
    output_dir='velvetk_plus_velvetopt'
    if kmer_length-half_range_size<=19: 
        min_k=19
        max_k=39
    elif kmer_length+half_range_size>=99: 
        min_k=79
        max_k=99
    else:
        min_k=kmer_length-half_range_size
        max_k=kmer_length+half_range_size
    logging.info("velvetk kmer range: %s-%s"%(min_k,max_k))
    return run_velvetOptimiser(fastq_file_name, min_k, max_k, output_dir)


def run_velvetk_plus_clc(fastq_file_name, estimated_size=600, **kwarg):
    kmer_length = run_velvetk(fastq_file_name, estimated_size)
    if kmer_length<12: kmer_length=12
    elif kmer_length>64: kmer_length=64
    output_dir='velvetk_plus_clc'
    return run_clc_assemble(fastq_file_name, kmer_length, output_dir)


def run_clc_assemble(fastq_file_name, word_size=None, output_dir='clc_bio', **kwarg):
    log_file='%s.log'%output_dir
    command='mkdir %s'%output_dir
    if not os.path.exists(output_dir):
        return_code = command_runner.run_command(command)
    if word_size:
        command = "clc_novo_assemble -v -w %s -q %s -o clc_bio/contigs.fa -m 100 2>&1 >%s "%(word_size, fastq_file_name, log_file)
    else:
        command = "clc_novo_assemble -v -q %s -o clc_bio/contigs.fa -m 100 2>&1 >%s "%(fastq_file_name, log_file)
    return_code = command_runner.run_command(command)
    contig_files = glob('%s/contigs.fa'%output_dir)
    contig_file_name=None
    if len(contig_files) == 1:
        contig_file_name = contig_files[0]
    return contig_file_name;


def run_abyss(fastq_file_name,tag_read_num=1, **kwarg):
    log_file='abyss.log'
    command='mkdir abyss'
    if not os.path.exists('abyss'):
        return_code = command_runner.run_command(command)
    kmer_length= 29

    command = "ABYSS -k %s %s -o abyss/contigs.fa --standard-quality -c 0 -e 0 --no-chastity 2>&1 >%s"%(kmer_length, fastq_file_name, log_file)
    return_code = command_runner.run_command(command)

    contig_files = glob('abyss/contigs.fa')
    contig_file_name=None
    if len(contig_files) == 1:
        contig_file_name = contig_files[0]
    
    return contig_file_name


def run_locasopt(fastq_file_name, **kwarg):
    log_file='locasopt.log'

    command = "LOCASopt -I %s -F fastq -O locasopt --noRevCompl -P kmer 11 15 2 -L 21 35 3 -ST 0.02 0.04 0.01 2>&1 >%s"%(fastq_file_name,log_file)
    return_code = command_runner.run_command(command)

    contig_files = glob('locasopt/contigs.fasta')
    contig_file_name=None
    if len(contig_files) == 1:
        contig_file_name = contig_files[0]
    
    return contig_file_name


def run_locas(fastq_file_name, **kwarg):
    log_file='locas.log'
    command = "locas -I %s -F fastq -O locas -C 101 1 2>&1 >%s"%(fastq_file_name,log_file)
    return_code = command_runner.run_command(command)
    contig_files = glob('locas/contigs.fasta')
    contig_file_name=None
    if len(contig_files) == 1:
        contig_file_name = contig_files[0]
    
    return contig_file_name


def run_soapdenovo(fastq_file_name, max_read_len=101, **kwarg):
    log_file='soapdenovo.log'
    command='mkdir soapdenovo'
    if not os.path.exists('soapdenovo'):
        return_code = command_runner.run_command(command)
    config_file='soapdenovo/config_file'
    open_file=open(config_file,'w')
    open_file.write("max_rd_len=%s\n[LIB]\nq=%s\n"%(max_read_len,fastq_file_name))
    open_file.close()
    command='SOAPdenovo31mer pregraph -K 29 -s %s -o soapdenovo/graph -p 1 2>&1 >%s'%(config_file,log_file)
    return_code = command_runner.run_command(command)
    command='SOAPdenovo31mer contig -g soapdenovo/graph 2>&1 >>%s'%log_file
    return_code = command_runner.run_command(command)
    contig_files = glob('soapdenovo/graph.contig')
    contig_file_name=None
    if len(contig_files) == 1:
        contig_file_name = contig_files[0]
    return contig_file_name


def run_idba(fastq_file_name, max_read_len=101, **kwarg):
    log_file='idba_ud.log'
    command="/ifs/software/linux_x86_64/idba_ud/idba_ud-1.0.9/bin/idba_ud -r %s -o idba_ud --min_contig %s --num_threads 1 2>&1 >%s"%(fastq_file_name, max_read_len,log_file)
    return_code = command_runner.run_command(command)
    contig_files = glob('idba_ud/contig.fa')
    contig_file_name=None
    if len(contig_files) == 1:
        contig_file_name = contig_files[0]
    return contig_file_name

def pass_function(*args, **kwarg):
    pass

def correct_contig_file(contig_file, site_name, min_contig_len=101):
    """This function will change the name of the contigs and remove the small one"""
    dir_name=os.path.dirname(contig_file)
    corrected_file=os.path.join(dir_name,'contigs_corrected.fa')
    open_file = open(contig_file)
    open_corrected=open(corrected_file,'w')
    fasta_reader = FastaReader(open_file)
    nb_seq=0
    max_len=0
    for fasta_record in fasta_reader:
        header, sequence=fasta_record
        if len(sequence)>min_contig_len:
            nb_seq+=1
            if len(sequence)>max_len:
                max_len=len(sequence)
            header='%s_pair_%s_length_%s'%(site_name,nb_seq,len(sequence))
            open_corrected.write('>%s\n%s\n'%(header,sequence))
    open_corrected.close()
    open_file.close()
    return (corrected_file, nb_seq, max_len)

def run_assembly(assembly_function, fastq_file, output_dir=None, estimated_size=600, name=None):
    if name is None:
        name,ext =os.path.splitext(os.path.basename(fastq_file))
    current_dir=None
    if output_dir and os.path.exists(output_dir):
        logging.debug('change directory to %s'%output_dir)
        current_dir=os.getcwd()
        os.chdir(os.path.abspath(output_dir))
    contig_file = assembly_function(fastq_file, estimated_size=estimated_size)
    if contig_file:
        contig_file = os.path.abspath(contig_file)
        merged_consensus = os.path.join(os.path.dirname(contig_file),'merged_consensus.fa')
        if os.path.exists(merged_consensus):
            logging.debug('remove the merged_consensus.fa that already exists before assembling')
            command = 'rm -f %s'%(merged_consensus)
            command_runner.run_command(command)
    
    if current_dir:
        logging.debug('change directory back to %s'%current_dir)
        os.chdir(current_dir)
    nb_seq=max_len=0
    corrected_contig_file=None
    if contig_file:
        corrected_contig_file, nb_seq, max_len = correct_contig_file(contig_file, name)
    return (corrected_contig_file,nb_seq, max_len)

def get_best_assembly(assembly_dir):
    contigs_file_to_compare = glob(os.path.join(assembly_dir,"*/contigs_corrected.fa"))
    contigs_file_to_compare.sort(cmp=compare_contig_file)
    return contigs_file_to_compare[0]

def get_best_assembly_merged(assembly_dir, read1_fasta, name, force_merge=False):
    contigs_file_to_compare = glob(os.path.join(assembly_dir,"*/contigs_corrected.fa"))
    contigs_file_to_compare.sort(cmp=compare_contig_file)
    
    best_assembly=None
    for contig_file2 in contigs_file_to_compare:
        output_dir=os.path.dirname(contig_file2)
        merged_contig_file = merge_read1_and_read2_contigs(name, read1_fasta, contig_file2)
        if merged_contig_file:
            best_assembly=merged_contig_file
            best_assembler_name=os.path.basename(os.path.dirname(contig_file2))
            logging.info("Best assembly with %s: Merged"%best_assembler_name)
            break
    
    if best_assembly is None :
        best_assembler_name=os.path.basename(os.path.dirname(contigs_file_to_compare[0]))
        logging.info("Best assembly with %s: Concatenated"%best_assembler_name)
        best_assembly=os.path.join(assembly_dir,"merged_consensus.fa")
        if force_merge:
            logging.info("Best assembly with %s: Force merged"%best_assembler_name)
            force_merge_consensus(read1_fasta, contigs_file_to_compare[0], best_assembly)
        else:
            concatenate_consensus([read1_fasta, contigs_file_to_compare[0]], best_assembly)
    return best_assembler_name,best_assembly

def compare_fasta_length(fasta_rec1,fasta_rec2):
    h1, s1 = fasta_rec1
    h2, s2 = fasta_rec2
    return len(s2)-len(s1)

def merge_read1_and_read2_contigs(name, read1_contig, read2_contigs):
    output_dir = os.path.dirname(read2_contigs)
    output_file = os.path.join(output_dir,'merged_consensus.fa')
    if os.path.exists(output_file):
        logging.warning('%s already exist no need to merge'%output_file)
        return output_file

    all_fasta2_entries=[]    
    with open(read2_contigs) as open_read2:
        read2_reader = FastaReader(open_read2)
        for header, sequence in read2_reader:
            all_fasta2_entries.append((header,sequence))
    if len(all_fasta2_entries)==1:
        results = merge_2_contigs(name, read1_contig, read2_contigs, output_dir)
        if results:
            with open(output_file,'w') as open_output:
                open_output.write('>%s\n%s\n'%results)
            return output_file
    else:
        all_fasta2_entries.sort(cmp=compare_fasta_length)
        merged_consensus = None
        remaining=[]
        for header,sequence in all_fasta2_entries:
            cur_pair=os.path.join(output_dir,header+".fa")
            with open(cur_pair,'w') as open_pair:
                open_pair.write(">%s\n%s\n"%(header,sequence))
            
            if not merged_consensus: 
                #If we have not successfully merged a read2 contig
                results = merge_2_contigs(name, read1_contig, cur_pair, output_dir)
                if not results:
                    remaining.append(cur_pair)
                else:
                    tmp_output_file=os.path.join(output_dir,'tmp_merged_consensus.fa')
                    with open(tmp_output_file,'w') as open_output:
                        open_output.write('>%s\n%s\n'%results)
                    merged_consensus = tmp_output_file                    
            else:
                results = merge_2_contigs(name+"add", read1_contig, cur_pair, output_dir)
                #TODO: fix this as it doesn't seems to trim and the output file is the same as above     but in the mean time disable
                if False and results:
                    additional_merged_pair=os.path.join(output_dir,'tmp_merged_consensus.fa')
                    with open(additional_merged_pair,'w') as open_output:
                        open_output.write('>%s\n%s\n'%results)
                    #trim this contig
                    trim_additional_merged_contigs(cur_pair, additional_merged_pair)
                remaining.append(cur_pair)
            
        if merged_consensus:
            tmp = [merged_consensus]
            tmp.extend(remaining)
            concatenate_consensus(tmp, output_file)
            logging.debug('Merging successfully in %s'%output_file)
            #TODO:delete the tmp file from the assembler directory
            return output_file
        else:
            logging.debug('Merging Failed in %s: %s'%(output_file, os.path.exists(output_file)))
            return None

def trim_additional_merged_contigs(original_contig, merged_contig):
    open_merged=open(merged_contig)
    fasta_reader=FastaReader(open_merged)
    header, sequence = fasta_reader.next()
    sp_header = header.split('_')
    trim=int(sp_header[-2].strip('os'))
    cigar=sp_header[-1]
    all_cigar = re.findall('(\d+)([MXDI])', cigar)
    for count, type  in all_cigar:
        if type=="M" or type=="X":
            trim+=int(count)
    
    open_merged.close()
    if trim>50:
        logging.info("trim %s of %s"%(trim, original_contig))
        open_contig=open(original_contig)
        fasta_reader=FastaReader(open_contig)
        header, sequence,fasta_reader.next()
        header=header+"trim_%s"%trim
        sequence=sequence[trim:]
        open_contig.close()
        open_contig=open(original_contig, 'w')
        if len(sequence)>100:
            open_contig.write(">%s\n%s\n"%(header,sequence))
        open_contig.close()


def concatenate_consensus(all_fasta_files,output_merge_file):
    if len(all_fasta_files) > 50:
        counter=1
        files_to_concatenate=[]
        for i in range(0,len(all_fasta_files), 50):
            output_merge_file_tmp="%s.%s"%(output_merge_file,counter)
            concatenate_consensus(all_fasta_files[i:i+50],output_merge_file_tmp)
            files_to_concatenate.append(output_merge_file_tmp)
            counter+=1
        concatenate_consensus(files_to_concatenate, output_merge_file)
        command='rm -f %s'%(' '.join(files_to_concatenate))
        command_runner.run_command(command)
    else:
        command = "cat %s > %s "%(' '.join(all_fasta_files), output_merge_file)
        command_runner.run_command(command)


def force_merge_consensus(read1_consensus, read2_consensus, output_merge_file):
    open_output = open(output_merge_file,'w')
    open_read1 = open(read1_consensus)
    open_read2 = open(read2_consensus)
    fasta_reader1 = FastaReader(open_read1)
    read1_name, read1_sequence = fasta_reader1.next()
    open_read1.close()
    name="%s_forced_merged"%read1_name
    array=[read1_sequence]
    fasta_reader2 = FastaReader(open_read2)
    print read2_consensus
    for read2_name, read2_sequence in fasta_reader2:
        print read2_sequence
        array.append("N"*100)
        array.append(read2_sequence)
    
    open_output.write(">%s\n%s\n"%(name, ''.join(array)))
    open_read2.close()
    open_output.close()


def get_list_of_length(contig_file):
    open_file = open(contig_file)
    list_length=[]
    nb_contig=0
    max_length=0
    fasta_reader = FastaReader(open_file)
    for fasta_record in fasta_reader:
        header, sequence=fasta_record
        nb_contig+=1
        if len(sequence)>max_length: max_length=len(sequence)
    return nb_contig, max_length


def compare_contig_file(contig_file1,contig_file2):
    nb_contig1, max_length1 = get_list_of_length(contig_file1)
    nb_contig2, max_length2 = get_list_of_length(contig_file2)
    if nb_contig1==nb_contig2:
        return max_length2-max_length1
    else:
        return nb_contig1-nb_contig2


def run_one_fastq_file(fastq_file, output_dir, assembly_function_list, estimated_size=600, read1_fasta=None, name=None, force_merge=False):
    fastq_file=os.path.abspath(fastq_file)
    #output_dir='%s_dir'%fastq_file
    if not os.path.exists(output_dir):
        command='mkdir %s'%(output_dir)
        return_code = command_runner.run_command(command)
    for assembly_function in assembly_function_list:
        #Assemble with provided assembler
        (contig_file, nb_seq, max_len) = run_assembly(assembly_function,fastq_file,output_dir,estimated_size=estimated_size, name=name)
        #Merge read one and read2 contig
        if contig_file:
            #TODO: This function gets run twice need to change that as the second run is not useful
            merge_read1_and_read2_contigs(name, read1_contig=read1_fasta, read2_contigs=contig_file)
        
    best_assembler_name, best_assembly_file = get_best_assembly_merged(output_dir, read1_fasta, name, force_merge)
    command="cp %s %s"%(best_assembly_file, os.path.join(output_dir, "best_assembly.fa"))
    return_code = command_runner.run_command(command)
    return os.path.join(output_dir, "best_assembly.fa")


def run_all_fastq_files(directory,assembly_function_list,estimated_size, force_merge=False):
    directory=os.path.abspath(directory)
    all_fastqs = glob(os.path.join(directory,'*','*_2.fastq'))

    all_contig_list=[]
    for fastq_file in all_fastqs:
        output_dir=os.path.dirname(fastq_file)
        name=os.path.basename(fastq_file)[:-len("_2.fastq")]
        read1_fasta=os.path.join(output_dir,name+"_1.fa")
        
        contig_file = run_one_fastq_file(fastq_file, output_dir, assembly_function_list, estimated_size=estimated_size, read1_fasta=read1_fasta, name=name)
        if contig_file:
            all_contig_list.append(contig_file)
        logging.info("\n")


def get_assembly_function(assembler_name):
    if assembler_name in VELVETOPT_ALIAS:
        return run_velvetOptimiser
    if assembler_name in VELVET_ALIAS:
        return run_velvet
    if assembler_name in VELVETK_ALIAS:
        return run_velvetk_plus_velvet
    if assembler_name in VELVETKOPT_ALIAS:
        return run_velvetk_plus_velvetopt
    if assembler_name in VELVETKCLC_ALIAS:
        return run_velvetk_plus_clc
    if assembler_name in CLC_ALIAS:
        return run_clc_assemble
    if assembler_name in ABYSS_ALIAS:
        return run_abyss
    if assembler_name in LOCASOPT_ALIAS:
        return run_locasopt
    if assembler_name in LOCAS_ALIAS:
        return run_locas
    if assembler_name in SOAPDENOVO_ALIAS:
        return run_soapdenovo
    if assembler_name in IDBA_ALIAS:
        return run_idba
    if assembler_name in NO_ASSEMBLER:
        return pass_function
    logging.error("Unknown assembler: %s"%assembler_name)
    return None

    
def main():
    #initialize the logging
    utils_logging.init_logging(logging.INFO)
    #utils_logging.init_logging(logging.CRITICAL)
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
    if not options.print_command:
        command_runner.set_command_to_run_localy()
    start_time = time.time()
    all_assembler_to_try = options.assembler_name.split(',')
    assembly_function_list=[]
    for assembler_name in all_assembler_to_try:
        assembly_function=get_assembly_function(assembler_name)
        assembly_function_list.append(assembly_function)
    if options.fastq_dir:
        run_all_fastq_files(options.fastq_dir,assembly_function_list, options.estimated_size, options.force_merge)
    elif options.fastq_file:
        name=os.path.basename(options.fastq_file)[:-len("_2.fastq")]
        output_dir=os.path.dirname(options.fastq_file)
        read1_fasta=os.path.join(output_dir,name+"_1.fa")
        contig_file = run_one_fastq_file(options.fastq_file, output_dir, assembly_function_list, estimated_size=options.estimated_size, read1_fasta=read1_fasta, name=name, force_merge=options.force_merge)
    logging.info("Elapsed time:%.1f seconds"%(time.time()-start_time))

def _prepare_optparser():
    """Prepare optparser object. New options will be added in this
    function first.
    """
    usage = """usage: %prog <-b bam_file> [ -o output_file]"""
    description = """This script will take aligned RAD read to the consensuses and calculate per consensus coverage."""
    
    optparser = OptionParser(version="None",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="show this help message and exit.")
    optparser.add_option("-f","--fastq_file",dest="fastq_file",type="string",
                         help="Path to one fastq_file. Default: %default")
    optparser.add_option("-d","--fastq_dir",dest="fastq_dir",type="string",
                         help="Path to a directory containing fastq file (only extension .fastq will be processed). Default: %default")
    optparser.add_option("-a","--assembler_name",dest="assembler_name",type="string",default=None,
                         help="The name of the assembler that will be used on the fastq files. Default: %default")
    optparser.add_option("-s","--estimated_size",dest="estimated_size",type="string",default="600",
                         help="The estimated size of the contig to assemble. It is used by velvetk to estimate the best kmer. Default: %default")
    optparser.add_option("--force_merge",dest="force_merge",action='store_true',default=False,
                         help="Force merged the consensus: add a run of 100 N in between each sequence. Default: %default")
    optparser.add_option("--print",dest="print_command",action='store_true',default=False,
                         help="print the commands instead of running them. Default: %default")
    optparser.add_option("--debug",dest="debug",action='store_true',default=False,
                         help="Output debug statment. Default: %default")
    return optparser


def _verifyOption(options):
    """Check if the mandatory option are present in the options objects.
    @return False if any argument is wrong."""
    arg_pass=True
    
    return arg_pass



if __name__=="__main__":
    main()
    
