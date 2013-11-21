from utils import utils_logging
import logging
import sys
from collections import Counter, defaultdict
import re
import copy
GENOTYPE='GT'
DEPTH='DP'
GENOTYPE_QUALITY='GQ'
GENOTYPE_SUPPORT='AD'

class NoSuchTagException(StandardError):
    pass

class NoSuchSampleException(StandardError):
    pass

class VcfRecord():
    def __init__(self, vcf_line, sample_names):
        """Contructor to read a vcf file header and content."""
        self.vcf_line_split=vcf_line.strip().split('\t')
        if len(self.vcf_line_split)>9:
            self.genotype_format = self.vcf_line_split[8]
            self.genotype_format_number = {}
            for pos,format in enumerate(self.genotype_format.split(':')):
                self.genotype_format_number[format]=pos
        else:
            self.genotype_format=None
        self.genotypes={}
        for pos, sample in enumerate(sample_names):
            self.genotypes[sample]=self.vcf_line_split[9+pos]
    
    def __str__(self):
        return '\t'.join(self.vcf_line_split)
        
    def get_reference(self):
        return self.vcf_line_split[0]
    
    def get_position(self):
        return int(self.vcf_line_split[1])
    
    def get_id(self):
        return self.vcf_line_split[2]
    
    def get_reference_base(self):
        return self.vcf_line_split[3]
    
    def get_alt_bases(self):
        return self.vcf_line_split[4].split(',')
    
    def get_qual(self):
        return float(self.vcf_line_split[5])
    
    def is_indel(self):
        """Check if this variant is an insertion of deletion"""
        if self.vcf_line_split[7].startswith('INDEL') or len(self.vcf_line_split[3])>1 or max(map(len, self.get_alt_bases()))>1:
            return True
        else:
            return False
    
    def get_info_tag_value(self,info_flag):
        """Extract the values of a specific tag from the info string."""
        #DP=16;AF1=0.4576;CI95=0.1667,0.6667;G3=0.118,0.8797,0.002267;DP4=3,1,2,0;MQ=23;FQ=9.46;PV4=1,0.48,1,1
        info=self.vcf_line_split[7]
        if info is None or info=='.':
            return None
        for tag_value in info.split(';'):
            if len(tag_value.split('='))==1:
                tag=tag_value
                value=True
            else:
                tag, value = tag_value.split('=')
            if tag == info_flag:
                return value
        return None
    
    def get_biais_pvalue(self):
        key='PV4'
        strand_pv=None
        baseQ_pv=None
        mapQ_pv=None
        tail_distance_pv=None
        
        value = self.get_info_tag_value(key)
        if value:
            strand_pv, baseQ_pv, mapQ_pv, tail_distance_pv=value.split(',')
        return (strand_pv, baseQ_pv, mapQ_pv, tail_distance_pv)
    
    def get_dp4_values(self):
        key='DP4'
        ref_forward=None
        ref_reverse=None
        alt_forward=None
        alt_reverse=None
        
        value = self.get_info_tag_value(key)
        if value:
            ref_forward, ref_reverse, alt_forward, alt_reverse=value.split(',')
        return (ref_forward, ref_reverse, alt_forward, alt_reverse)
        
    def get_genotype_quality(self,sample):
        """Extract the values of the genotype quality from the genotype string for a sample."""
        try :
            value = self.get_genotype_tag_value(GENOTYPE_QUALITY, sample)
            if value:
                return float(value)
            else:
                return 0
        except NoSuchTagException:
            return 0
    
    def get_sample_depth(self,sample):
        """Extract the values of the depth from the depth string for a sample."""
        try :
            return self.get_genotype_tag_value(DEPTH, sample)
        except NoSuchTagException:
            return None
        
    def get_all_genotype_quality(self,sample_list=None):
        """Extract the values of the genotype quality from the genotype string
        The order of the genotype qualities can be  different from the order in the file."""
        return self.get_all_genotype_tag_value(GENOTYPE_QUALITY, sample_list=sample_list)
    
    def get_all_sample_depth(self,sample_list=None):
        """Extract the depth of sequence covering this sample from the genotype string
        The order of the depth values can be  different from the order in the file."""
        return self.get_all_genotype_tag_value(DEPTH, sample_list=sample_list)
    
    def get_genotype(self,sample):
        """Extract the values of the genotype from the genotype string for a sample."""
        return self.get_genotype_tag_value(GENOTYPE, sample)
        
    def get_all_genotype(self,sample_list=None):
        """Extract the values of the genotype from the genotype string
        The order of the genotype can be  different from the order in the file."""
        return self.get_all_genotype_tag_value(GENOTYPE, sample_list=sample_list)

    def get_genotype_tag_value(self, genotype_tag, sample):
        """Extract the values of a specific tag from the genotype string for a sample."""
        genotype_string = self.genotypes.get(sample)
        if genotype_string is None:
            raise NoSuchSampleException("Sample %s not in vcf header"%sample)
        if genotype_string=='.' or genotype_string=='./.':
            return None
        genotype_elements = genotype_string.split(':')
        if not self.genotype_format_number.has_key(genotype_tag):
            raise NoSuchTagException("Tag %s not in genotype format"%(genotype_tag))
        if self.genotype_format_number[genotype_tag]<len(genotype_elements):
            return genotype_elements[self.genotype_format_number[genotype_tag]]
        else:
            return None
    
    def get_all_genotype_tag_value(self, genotype_tag, sample_list=None):
        """Extract the values of a specific tag from the genotype string for all samples.
        The order of the values can be  different from the order in the file."""
        all_values=[]
        if not self.genotype_format_number.has_key(genotype_tag):
            raise NoSuchTagException
        if sample_list:
            genotype_list=[]
            for sample in sample_list:
                genotype_list.append(self.genotypes.get(sample))
        else:
            genotype_list=self.genotypes.values()
        for genotype in genotype_list:
            if genotype=='.' or genotype=='./.':
                all_values.append(None)
            else:
                all_values.append(genotype.split(':')[self.genotype_format_number[genotype_tag]])
        return all_values
    
    def get_valid_genotype_per_sample(self, genotype_quality_threshold=0, minimum_depth=0, sample_list=None):
        """return a list of valid genotype and their associated sample name."""
        all_genotypes = {}
        if sample_list is None:
            sample_list=self.genotypes.keys()
        
        for sample in sample_list:
            genotype_quality = self.get_genotype_quality(sample)
            depth = self.get_sample_depth(sample)
            if genotype_quality is not None and float(genotype_quality) >= genotype_quality_threshold and\
            depth is not None and int(depth)>=minimum_depth:
                genotype = self.get_genotype(sample)
                if not all_genotypes.has_key(genotype):
                    all_genotypes[genotype]=[]
                all_genotypes[genotype].append(sample)
        return all_genotypes
    
    def get_genotype_support_for_ref_and_alt_allele(self,sample):
        """Extract the number of bases seen supporting reference and alternate allele for a sample."""
        try :
            sample_support =  self.get_genotype_tag_value(GENOTYPE_SUPPORT, sample)
            if sample_support and sample_support!='.':
                ref_support, alt_support = sample_support.split(',')
                return (int(ref_support), int(alt_support))
            else:
                return (0,0)
        except NoSuchTagException:
            return (0,0)
    
    def get_bases_support_for_ref_and_alt_allele(self, sample_list=None):
        """return a list of valid genotype and their associated sample name."""
        all_samples_support = []
        if sample_list is None:
            sample_list=self.genotypes.keys()
        
        for sample in sample_list:
            ref_support, alt_support = self.get_genotype_support_for_ref_and_alt_allele(sample)
            all_samples_support.append((ref_support, alt_support))
        return sample_list, all_samples_support
    
    def is_phased_with_previous(self, sample):
        geno = self.get_genotype(sample)
        if geno: return len(self.get_genotype(sample).split('|'))>1
        else: False
            
    def set_id(self,id):
        self.vcf_line_split[2]=id
    
        
class VcfReader():
    def __init__(self, file_handle):
        """Contructor to read a vcf file header and content."""
        self.file_handle=file_handle
        self._parse_header()
        
        #TODO: Can specify what type of iterator we want based on options passed to the constructor
        self.iterator=None
        
    def __iter__(self):
        """Allow the object to be used as an iterator"""
        if not self.iterator:
            self.iterator=self._get_vcfrecord()
        return self.iterator
    
    def __del__(self):
        pass
        
    def _parse_header(self):
        """Parse the header of a vcf file to get the metadata associated"""
        #TODO need some more information about the header
        self.header_lines=[]
        for line in self.file_handle:
            self.header_lines.append(line.strip())
            if '#CHROM' in line:
                self.all_fields = line.strip().split('\t')
                self.all_samples = self.all_fields[9::]
                break
        #self.all_samples = get_vcfsamples(self.file_handle)
    
    def _get_vcfdataline(self):
        """return all entries line by line."""
        for line in self.file_handle:
            if not line.startswith('#') : yield line.strip()
            
    def _get_vcfrecord(self):
        """return all entries record by record."""
        for line in self.file_handle:
            if not line.startswith('#') : yield VcfRecord(line.strip(),self.all_samples)
    
    def get_sample_names(self):
        """return the sample name from the header"""
        return self.all_samples

    def get_header_lines(self):
        """return the header lines in one string"""
        return '\n'.join(self.header_lines)
        

class PhasedVcfRecord():
    def __init__(self, vcf_records, all_samples, sample_references=None):
        """Contructor to create one phased vcf record out of several vcf Records."""
        self.vcf_records=vcf_records
        self.all_samples=all_samples
        self.vcf_sample_references=sample_references
        self.all_alternates=[]
        self._construct_phased_genotype()
        
    def __iter__(self):
        """Allow the object to be used as an iterator"""
        if not self.iterator:
            self.iterator=iter(self.vcf_records)
            return self.iterator

    def _construct_phased_genotype(self):
        """This method will create a list of genotype one for each sample along with minimum genotype quality and depth encountered."""
        all_samples_genotype=[]
        for sample in self.all_samples:
            all_samples_genotype.append([])
        haplotype1_per_sample=defaultdict(list)
        haplotype2_per_sample=defaultdict(list)
        genotypes_quality_per_sample=defaultdict(list)
        genotypes_depth_per_sample=defaultdict(list)
        
        for i, vcf_record in enumerate(self.vcf_records):
            for sample in self.all_samples:
                genotype = vcf_record.get_genotype(sample)
                genotype_quality = vcf_record.get_genotype_quality(sample)
                genotype_depth = vcf_record.get_sample_depth(sample)
                #Check that this genotype is phased with the previous one unless we're looking at the first one
                if i!=0 and not vcf_record.is_phased_with_previous(sample):
                    genotype_quality=0
                    genotype_depth=0
                if genotype:
                    h1,h2 = re.split(r'[|/]', genotype)
                else:
                    #if genotype is not defined the haplotype won't be either
                    h1,h2,genotype_quality,genotype_depth=('.','.',0,0)
                haplotype1_per_sample[sample].append(h1)
                haplotype2_per_sample[sample].append(h2)
                genotypes_quality_per_sample[sample].append(genotype_quality)
                genotypes_depth_per_sample[sample].append(genotype_depth)
        all_haplotypes={}
        self.non_collapse_genotypes={}
        self.genotypes={}
        self.genotypes_depth={}
        self.genotypes_quality={}
        for sample in self.all_samples:
            h1=''.join(haplotype1_per_sample.get(sample))
            if '.' in h1:
                haplotype_code1='.'
            elif not h1 in all_haplotypes:
                haplotype_code1 = len(all_haplotypes)
                all_haplotypes[h1]=haplotype_code1
                alternate_bases = self._haplotype_to_alternates(haplotype1_per_sample.get(sample)) 
                self.all_alternates.append('_'.join(alternate_bases))
            else:
                haplotype_code1 = all_haplotypes.get(h1)
            
            h2=''.join(haplotype2_per_sample.get(sample))
            if '.' in h2:
                haplotype_code2='.'
            elif not h2 in all_haplotypes:
                haplotype_code2 = len(all_haplotypes)
                all_haplotypes[h2]=haplotype_code2
                alternate_bases = self._haplotype_to_alternates(haplotype2_per_sample.get(sample))
                self.all_alternates.append('_'.join(alternate_bases))
            else:
                haplotype_code2 = all_haplotypes.get(h2)
            
            ref_base = self.get_reference_base()
            if ref_base in self.all_alternates:
                self.all_alternates.remove(ref_base)
                
            self.non_collapse_genotypes[sample]='%s/%s'%(h1,h2)
            self.genotypes[sample]='%s/%s'%(haplotype_code1,haplotype_code2)
            self.genotypes_depth[sample]=min(genotypes_depth_per_sample.get(sample))
            self.genotypes_quality[sample]=min(genotypes_quality_per_sample.get(sample))
            
    def _haplotype_to_alternates(self,haplotype_array):
        alternates = []
        for i, haplotype in enumerate(haplotype_array):
            if haplotype == '0' :
                base = self.vcf_records[i].get_reference_base()
            else:
                base = self.vcf_records[i].get_alt_bases()[int(haplotype)-1]
            alternates.append(base)
        return alternates
            
    
    def __str__(self):
        """String representation of Phased vcf records"""
        out=[]
        out.append(self.get_reference())
        out.append(str(self.get_position()))
        out.append('.')
        out.append(str(self.get_reference_base()))
        out.append(','.join(self.get_alt_bases()))
        out.append(str(self.get_qual()))
        out.append('.')
        out.append('.')
        out.append('GT:OG:DP:GQ')
        for sample in self.all_samples:
            out.append('%s:%s:%s:%s'%(self.genotypes.get(sample), self.non_collapse_genotypes.get(sample),
                                      self.genotypes_depth.get(sample), self.genotypes_quality.get(sample) ) )
        return '\t'.join(out)

    def get_reference(self):
        return self.vcf_records[0].get_reference()
    
    def get_position(self):
        return min([v_rec.get_position() for v_rec in self.vcf_records])
    
    def get_id(self):
        raise NotImplementedError("Method get_id is not implemented for VcfPhasedRecord")
    
    def get_reference_base(self):
        return '_'.join(rec.get_reference_base() for rec in self.vcf_records)
    
    def get_alt_bases(self):
        return copy.copy(self.all_alternates)
    
    def get_qual(self):
        return min([v_rec.get_qual() for v_rec in self.vcf_records])
    
    def is_indel(self):
        ret=False
        for rec in self.vcf_records:
            if rec.is_indel(): 
                ret=True
                break
        return ret

    def get_info_tag_value(self,info_flag):
        raise NotImplementedError("Method get_info_tag_value is not implemented for VcfPhasedRecord")
    
    def get_biais_pvalue(self):
        raise NotImplementedError("Method get_biais_pvalue is not implemented for VcfPhasedRecord")
    
    def get_dp4_values(self):
        raise NotImplementedError("Method get_dp4_values is not implemented for VcfPhasedRecord")
    
    def get_genotype(self,sample):
        """Returns the phased genotype of the marker."""
        return self.genotypes.get(sample)
    
    
    def get_genotype_quality(self,sample):
        """Returns the phased genotype quality of the marker."""
        return self.genotypes_quality.get(sample)
        
    
    def get_sample_depth(self,sample):
        """Returns the depth of coverage of the marker"""
        return self.genotypes_depth.get(sample)
    
    def get_all_genotype(self,sample_list=None):
        """Return all the genotype quality for a list of sample"""
        genotypes=[]
        for sample in sample_list:genotypes.append(self.get_genotype(sample))
        return genotypes
    
    def get_all_genotype_quality(self,sample_list=None):
        """Return all the genotype quality for a list of sample"""
        genotype_qualities=[]
        for sample in sample_list:genotype_qualities.append(self.get_genotype_quality(sample))
        return genotype_qualities
    
    def get_all_sample_depth(self,sample_list=None):
        """Return all the depth information for a list of sample"""
        depth_info=[]
        for sample in sample_list:depth_info.append(self.get_depth(sample))
        return depth_info
        
    
    def get_genotype_tag_value(self, genotype_tag, sample):
        raise NotImplementedError("Method get_genotype_tag_value is not implemented for VcfPhasedRecord")
    
    
    def get_all_genotype_tag_value(self, genotype_tag, sample_list=None):
        raise NotImplementedError("Method get_all_genotype_tag_value is not implemented for VcfPhasedRecord")
    
    def get_valid_genotype_per_sample(self, genotype_quality_threshold=0, minimum_depth=0, sample_list=None):
        """return a list of valid genotype and their associated sample name."""
        all_genotypes = {}
        if sample_list is None:
            sample_list=self.genotypes.keys()
        
        for sample in sample_list:
            genotype_quality = self.get_genotype_quality(sample)
            depth = self.get_sample_depth(sample)
            if genotype_quality is not None and float(genotype_quality) >= genotype_quality_threshold and\
            depth is not None and int(depth)>=minimum_depth:
                genotype = self.get_genotype(sample)
                if not all_genotypes.has_key(genotype):
                    all_genotypes[genotype]=[]
                all_genotypes[genotype].append(sample)
        return all_genotypes

def get_list_of_SNPs_name(list_vcf_record):
    ref=None
    positions=[]
    for vcf_record in list_vcf_record:
        if ref and ref!=vcf_record.get_reference():
            logging.error("%s and %s are different: can only get a name to SNPs of the same reference"%(ref,vcf_record.get_reference()))
            return None
        else:
            ref=vcf_record.get_reference()
        positions.append(str(vcf_record.get_position()))
    return '%s:%s'%(ref,'_'.join(positions))

class PhasedVcfReader():
    def __init__(self, file_handle, reference_samples=None, number_mandatory_sample=None, geno_qual_threshold=20):
        
        self.file_handle = file_handle
        self.reader  = VcfReader(self.file_handle)
        self.reader.header_lines.append('##FORMAT=<ID=OG,Number=1,Type=String,Description="Original Genotypes">')
        self.geno_qual_threshold=geno_qual_threshold
        self.all_samples=self.reader.get_sample_names()
        
        if reference_samples:
            logging.debug('Use %s (%s) samples as reference'%(len(reference_samples), ', '.join(reference_samples) ))
            self.reference_samples=reference_samples
            self.number_mandatory_sample=len(reference_samples)
        elif number_mandatory_sample:
            logging.debug('Use %s samples out of %s as reference'%(number_mandatory_sample, len(self.all_samples) ))
            self.reference_samples=self.all_samples
            if number_mandatory_sample>len(self.all_samples):self.number_mandatory_sample=len(self.all_samples)
            else:self.number_mandatory_sample=number_mandatory_sample
        else:
            logging.debug('Use all (%s) samples as reference'%(len(self.all_samples) ))
            self.reference_samples=self.all_samples
            self.number_mandatory_sample=len(self.all_samples)
            
        self.nb_snps=0
        self.marker_created=0
        self.discard_quality_reference_samples=0
        self.remaining_sample=self.reader.get_sample_names()
        #for sample in self.reference_samples:
        #    if sample in self.remaining_sample:
        #        self.remaining_sample.remove(sample)
        #    else:
        #        raise NoSuchSampleException("Sample %s not in vcf header"%sample)
        self.iterator=None
        
    def __iter__(self):
        """Allow the object to be used as an iterator"""
        if not self.iterator:
            self.iterator=self._get_phased_vcfrecord()
            return self.iterator
            
     
    def _get_phased_vcfrecord(self):
        end_of_phase=False
        phased_vcf_record=[]
        previous_reference=None
        start_session_reference=Counter()
        
        for vcf_records in self.reader:
            
            logging.debug("Try phasing %s:%s with %s previous SNPs "%(vcf_records.get_reference(),vcf_records.get_position(), len(phased_vcf_record)))
            #First check that enough reference_samples are callable
            self.nb_snps+=1
            quality_failed_samples=[]
            in_phase_samples=[]
            for sample in self.reference_samples:
                #geno_sample = vcf_records.get_genotype(sample)
                gq_sample = vcf_records.get_genotype_quality(sample)
                if gq_sample < self.geno_qual_threshold:
                    quality_failed_samples.append(sample)
                
                if vcf_records.is_phased_with_previous(sample=sample):
                    in_phase_samples.append(sample)
            
            #Test if enough reference samples are good quality for that SNP
            if len(quality_failed_samples) > len(self.reference_samples) - self.number_mandatory_sample:
                self.discard_quality_reference_samples+=1
                #break the Phase if we the quality is bad
                #yield this SNP on its own
                add_SNP_before_yield=False
                add_SNP_after_yield=True
                end_of_phase=True
                logging.debug("stop phasing because %s:%s is low quality"%(vcf_records.get_reference(),vcf_records.get_position()))
            elif len(in_phase_samples) >= self.number_mandatory_sample:
                #Keep the Phase if enough samples are keeping the phase
                #Add this SNP to the phased set
                add_SNP_before_yield=True
                add_SNP_after_yield=False
                end_of_phase=False
                
            else:
                #Break the Phase if too many samples break the phase
                #Add this SNP to a new phased set
                add_SNP_before_yield=False
                add_SNP_after_yield=True
                #to correctly phase homozygous location at the beginning of contigs records if an het has occured  
                for sample in self.reference_samples:
                    if not vcf_records.is_phased_with_previous(sample=sample):
                        start_session_reference[sample]=1
                if sum(start_session_reference.values()) >= self.number_mandatory_sample:
                    logging.debug("stop phasing because not phased with previous %s:%s"%(vcf_records.get_reference(), vcf_records.get_position()))
                    end_of_phase=True
                else:
                    logging.debug("Continue phasing because only %s sessions started so far %s:%s"%(sum(start_session_reference.values()), vcf_records.get_reference(), vcf_records.get_position()))
                    end_of_phase=False
                
            
            
            if previous_reference and previous_reference!=vcf_records.get_reference():
                #break the Phase if we change reference
                add_SNP_before_yield=False
                add_SNP_after_yield=True
                end_of_phase=True
                logging.debug("stop phasing because changing reference: now %s --> %s"%(previous_reference, vcf_records.get_reference()))
                for sample in self.reference_samples: start_session_reference[sample]=False

            
            if vcf_records.get_info_tag_value('PhasingInconsistent'):
                #break the phase if GATK says the phasing is inconsistent
                logging.debug("stop phasing because Found a PhasingInconsistent mark: %s:%s"%(vcf_records.get_reference(), vcf_records.get_position()))
                end_of_phase=True
                
            if add_SNP_before_yield:
                #phased_genotype_mother.append(re.split(r'[|/]',geno_mother))
                #phased_genotype_father.append(re.split(r'[|/]',geno_father))
                phased_vcf_record.append(vcf_records)
            #print previous_reference, vcf_records.get_reference() , vcf_records.get_position()
            #print 'in phase for mother=%s, in phase for father=%s, end_of_phase=%s, add_SNP_before_yield=%s, add_SNP_after_yield=%s'%(
            #                                                         vcf_records.is_phased_with_previous(sample=self.mother),
            #                                                         vcf_records.is_phased_with_previous(sample=self.father),
            #                                                         end_of_phase, add_SNP_before_yield,add_SNP_after_yield)
            #print 'start_session_mother=%s, start_session_father=%s'%(start_session_mother,start_session_father)
            #print 'nb_of SNPs in phased set = %s'%len(phased_vcf_record)
            #if  len(phased_vcf_record)>0:
            #    print ' || '.join(['%s--%s'%(tmp.get_reference(),tmp.get_position()) for tmp in phased_vcf_record])
            if end_of_phase and len(phased_vcf_record)>0:
                self.marker_created+=1
                logging.debug("Create marker: %s"%(get_list_of_SNPs_name(phased_vcf_record)))
                yield PhasedVcfRecord(phased_vcf_record,self.all_samples)
                phased_vcf_record=[]    
                
            if add_SNP_after_yield:
                phased_vcf_record.append(vcf_records)
            previous_reference = vcf_records.get_reference()
        if len(phased_vcf_record)>0:
            self.marker_created+=1
            yield PhasedVcfRecord(phased_vcf_record,self.all_samples)

    def get_sample_names(self):
        """return the sample name from the header"""
        return self.reader.get_sample_names()

    def get_header_lines(self):
        """return the header lines in one string"""
        return self.reader.get_header_lines()

    def get_nb_marker_created(self):
        return self.marker_created
    
    def get_discard_quality_parent(self):
        return self.discard_quality_parent
    
    def get_nb_snps(self):
        return self.nb_snps
 



if __name__=="__main__":
    utils_logging.init_logging(logging.DEBUG)
    input_vcf = sys.argv[1]
    reader = PhasedVcfReader(input_vcf, reference_samples=None, number_mandatory_sample=None, geno_qual_threshold=20)
    for val in reader:
        print val
        print 
