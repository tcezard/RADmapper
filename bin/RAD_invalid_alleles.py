'''
Created on Mar 9, 2011

@author: tcezard
'''
from utils import utils_logging, utils_commands
from IO_interface import vcfIO
import sys, os
import pprint
import utils
from IO_interface.samIterator import Sam_record


def generate_empty_hash_with_sample(all_samples):
    hash={}
    for sample in all_samples:
        hash[sample]={}
    hash['all']={}
    return hash

def count_with_hash(hashtable, key):
    if hashtable.has_key(key):
        hashtable[key]+=1
    else:
        hashtable[key]=1

def snps_to_allele(vcf_file, bam_file):
    
    file_handle = utils_logging.open_input_file(vcf_file, pipe=False)
    reader  = vcfIO.VcfReader(file_handle)
    sample_names = reader.get_sample_names()
    #all_samples_in_file=reader.get_sample_names()
    vcf_record_in_one_contig={}
    curr_reference=None
    for vcf_records in reader:
        #First check that the parent are callable
        if curr_reference!=vcf_records.get_reference():
            if curr_reference:
                process_alleles(vcf_record_in_one_contig, sample_names, curr_reference)
            vcf_record_in_one_contig={}
        vcf_record_in_one_contig[vcf_records.get_position()]=vcf_records
        curr_reference = vcf_records.get_reference()
    if curr_reference:
        process_alleles(vcf_record_in_one_contig, sample_names, curr_reference)
    

def process_alleles(vcf_record_in_one_contig, sample_names,curr_reference):
    sample_to_allele = generate_empty_hash_with_sample(sample_names)
    command = "samtools view -F 1028 %s %s"%(bam_file,curr_reference)
    stream,process=utils_commands.get_output_stream_from_command(command)
    for line in stream:
        sam_record=Sam_record(line)
        allele_array=[]
        sequence = sam_record.get_query_sequence()
        sample = sam_record.get_tag("RG")
        for position in vcf_record_in_one_contig.keys():
            #if vcf_record_in_one_contig.get(position).get_genotype_quality(sample)>20:
            allele_array.append(sequence[position-1])
            #else:
            #    allele_array.append('.')
        count_with_hash(sample_to_allele[sample], ''.join(allele_array))
        count_with_hash(sample_to_allele['all'], ''.join(allele_array))
    process.wait()
    pprint.pprint(sample_to_allele)
    filter_alleles(sample_to_allele)
    pprint.pprint(sample_to_allele)
    all_alleles=set()
    valid=True
    for sample in sample_to_allele.keys():
        alleles = sample_to_allele.get(sample)
        all_alleles.update(set(alleles.keys()))
        if len(alleles)>2:
            valid=False
    if len(all_alleles)>4:
        valid=False
    if not valid:
        print curr_reference

def filter_alleles(sample_to_allele):
    all_alleles_count = sample_to_allele.pop('all')
    total=0.0
    valid_alleles=[]
    for alleles in all_alleles_count:
        total+=all_alleles_count.get(alleles)
    for allele in all_alleles_count:
        if all_alleles_count.get(allele)>1 and all_alleles_count.get(allele)/total>0.05:
            valid_alleles.append(allele)
    for sample in sample_to_allele.keys():
        test_alleles = sample_to_allele.get(sample)
        total=0.0
        for alleles in test_alleles:
            total+=test_alleles.get(alleles)
        for allele in test_alleles.keys():
            if not allele in valid_alleles or test_alleles.get(allele)<3 or test_alleles.get(allele)/total<=0.05 :
                test_alleles.pop(allele)



if __name__=="__main__":
    """This script will output a list of contig that have more than 4 alleles on the whole set and more than 2 per individuals"""
    vcf_file = sys.argv[1]
    bam_file = sys.argv[2]
    
    snps_to_allele(vcf_file,bam_file)
