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
from RAD_genotype_to_pattern import Phased_vcf_reader, load_sex_info


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

def allele_presence_abscence(vcf_file,parent_name1,parent_name2):
    
    file_handle = utils_logging.open_input_file(vcf_file, pipe=False)
    reader  = vcfIO.VcfReader(file_handle)
    sample_names = reader.get_sample_names()
    #all_samples_in_file=reader.get_sample_names()
    vcf_record_in_one_contig={}
    curr_reference=None
    qual_threshold=10
    cov_threshold=0
    
    children_samples=[]
    for sample in sample_names:
        if sample != parent_name1 and sample != parent_name2:
            children_samples.append(sample) 
    for vcf_records in reader:
        #First check that the parent are callable
        vcf_record_in_one_contig[vcf_records.get_position()]=vcf_records
        curr_reference = vcf_records.get_reference()
        genotype_p1 = vcf_records.get_valid_genotype_per_sample(genotype_quality_threshold=qual_threshold, minimum_depth=cov_threshold,sample_list=[parent_name1])
        genotype_p2 = vcf_records.get_valid_genotype_per_sample(genotype_quality_threshold=qual_threshold, minimum_depth=cov_threshold,sample_list=[parent_name2])
        
        #if the parent are callable
        
        if len(genotype_p1)==1 and len(genotype_p2)==1:
            g1=genotype_p1.keys()[0]
            g2=genotype_p2.keys()[0]
            all_alleles=set()
            all_alleles.update(set(g1.split('/')))
            all_alleles.update(set(g2.split('/')))
            if g1!=g2 and len(all_alleles)==2 and ( len(set(g1.split('/')))>1 or len(set(g2.split('/')))>1 ):
                allele1=set(all_alleles.pop())
                allele2=set(all_alleles.pop())
                allele1_seg_pattern=[]
                allele2_seg_pattern=[]
                
                
                for sample in children_samples:
                    remaining_genotype = vcf_records.get_valid_genotype_per_sample(genotype_quality_threshold=qual_threshold,
                                                                                   minimum_depth=cov_threshold,sample_list=[sample])
                    if len(remaining_genotype)>0:
                        haplotypes = set(remaining_genotype.keys()[0].split('/'))
                        if len(allele1.intersection(haplotypes)):
                            allele1_seg_pattern.append('1')
                        else:
                            allele1_seg_pattern.append('0')
                        if len(allele2.intersection(haplotypes)):
                            allele2_seg_pattern.append('1')
                        else:
                            allele2_seg_pattern.append('0')
                    else:
                        allele1_seg_pattern.append('-')
                        allele2_seg_pattern.append('-')
                if len(set(g1.split('/')))>1:
                    print "%s\t%s\t10\t%s"%(curr_reference, vcf_records.get_position(),''.join(allele1_seg_pattern))
                elif len(set(g2.split('/')))>1:
                    print "%s\t%s\t01\t%s"%(curr_reference, vcf_records.get_position(),''.join(allele2_seg_pattern))
                    
                        
                                 
                            
                            

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
            if not allele in valid_alleles or test_alleles.get(allele)==1 or test_alleles.get(allele)/total<=0.05 :
                test_alleles.pop(allele)

def sex_specific_markers(vcf_file, mother, father, offsprings_file):
    sample_to_sex,sex_to_sample, ordered_sample = load_sex_info(offsprings_file)
    file_handle = utils_logging.open_input_file(vcf_file, pipe=False)
    reader  = vcfIO.VcfReader(file_handle)
    sex_to_sample.get("M").remove(father)
    sex_to_sample.get("F").remove(mother)
    
    for vcf_record in reader:
        gt_mother = vcf_record.get_genotype(mother)
        sp_mother = vcf_record.get_sample_depth(mother)
        gt_father = vcf_record.get_genotype(father)
        sp_father = vcf_record.get_sample_depth(father)
        sample_female=[]
        sample_male=[]
        if (sp_mother is None or int(sp_mother)<4) and sp_father is not None and int(sp_father)>10:
            valid=True
            for sample in sex_to_sample.get('F'):
                gt_of_fem = vcf_record.get_genotype(sample)
                sp_of_fem = vcf_record.get_sample_depth(sample)
                sample_female.append("%s:%s"%(gt_of_fem,sp_of_fem))
                if sp_of_fem is not None and int(sp_of_fem)>2:
                    valid=False
                    break
            nb_male_offspring=0
            for sample in sex_to_sample.get('M'):
                gt_of_mal = vcf_record.get_genotype(sample)
                sp_of_mal = vcf_record.get_sample_depth(sample)
                sample_male.append("%s:%s"%(gt_of_mal,sp_of_mal))
                if sp_of_mal is not None and int(sp_of_mal)>5:
                    nb_male_offspring+=1
            if nb_male_offspring < len(sex_to_sample.get('M'))-4:
                valid=False
            if valid :
                print vcf_record.get_reference(), vcf_record.get_position(), "%s:%s"%(gt_father,sp_father), '  '.join(sample_male), "\t\t", "%s:%s"%(gt_mother,sp_mother), '  '.join(sample_female)

def sex_specific_markers_female(vcf_file, mother, father, offsprings_file):
    sample_to_sex,sex_to_sample, ordered_sample = load_sex_info(offsprings_file)
    file_handle = utils_logging.open_input_file(vcf_file, pipe=False)
    reader  = vcfIO.VcfReader(file_handle)
    sex_to_sample.get("M").remove(father)
    sex_to_sample.get("F").remove(mother)

    for vcf_record in reader:
        gt_mother = vcf_record.get_genotype(mother)
        sp_mother = vcf_record.get_sample_depth(mother)
        gt_father = vcf_record.get_genotype(father)
        sp_father = vcf_record.get_sample_depth(father)
        sample_female=[]
        sample_male=[]

        if (sp_mother is None or int(sp_mother)<4) and sp_father is not None and int(sp_father)>10:
            valid=True
            for sample in sex_to_sample.get('F'):
                gt_of_fem = vcf_record.get_genotype(sample)
                sp_of_fem = vcf_record.get_sample_depth(sample)
                sample_female.append("%s:%s"%(gt_of_fem,sp_of_fem))
                if sp_of_fem is not None and int(sp_of_fem)>2:
                    valid=False
                    break
            nb_male_offspring=0
            for sample in sex_to_sample.get('M'):
                gt_of_mal = vcf_record.get_genotype(sample)
                sp_of_mal = vcf_record.get_sample_depth(sample)
                sample_male.append("%s:%s"%(gt_of_mal,sp_of_mal))
                if sp_of_mal is not None and int(sp_of_mal)>5:
                    nb_male_offspring+=1
            if nb_male_offspring < len(sex_to_sample.get('M'))-4:
                valid=False
            if valid :
                print vcf_record.get_reference(), vcf_record.get_position(), "%s:%s"%(gt_father,sp_father), '  '.join(sample_male), "\t\t", "%s:%s"%(gt_mother,sp_mother), '  '.j$



if __name__=="1__main__":
    """This script will output a list of contig that have more than 4 alleles on the whole set and more than 2 per individuals"""
    vcf_file = sys.argv[1]
    parent_1 = sys.argv[2]
    parent_2 = sys.argv[3]
    
    allele_presence_abscence(vcf_file,parent_1,parent_2)
    
if __name__=="__main__":
    """This script will output a list of contig that have more than 4 alleles on the whole set and more than 2 per individuals"""
    vcf_file = sys.argv[1]
    mother = sys.argv[2]
    father = sys.argv[3]
    offsprings_file = sys.argv[4]
    sex_specific_markers(vcf_file, mother, father, offsprings_file)
