import sys
import re
import math
from utils import utils_logging, benchmarck_timer
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
        return self.vcf_line_split[7].startswith('INDEL')

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
        return len(self.get_genotype(sample).split('|'))>1
    
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
        


def compareSampleLists(list1, list2):
    """ given two list of samplenames, compare element by element; return 1 if they are all the same 0 otherwise """
    if len(list1) != len(list2): return 0

    for i in range(0, len(list1) ):
        if list1[i] != list2[i]: return 0

    return 1

def get_vcfdataline(fh):
    """yield datalines """
    for line in fh:
        if '#' not in line: yield line.strip()
        

def yield_bedinterval(fh):
    """ yield tuple of (chr, start, end) zero-based half-open interval of vcf dataline """
    for line in fh:
        if '#' not in line:
            fields=split_vcfdataline( line.strip () )
            (chr,start)=fields[0:2]
            yield (chr, int(start)-1, int(start) )


def split_vcfdataline(line):
    """ split the fields of a vcf dataline and return in a list """
    return line.strip().split('\t')

def get_vcfdataline_passfilter(fh):
    """ yield list that has not been filtered """
    for line in fh:
        if '#' not in line:
            fields = line.strip().split('\t')
            if fields[6] == '.' or fields[6]=='PASS':
                yield line.strip()

def get_vcfdatafields(fh):
    """ yield data fields """
    for line in fh:
        if '#' not in line:
            fields = line.strip().split('\t')
            yield fields[0:9]

def get_vcfgenotypes(fh):
    """ yield genotypes from vcf dataline """
    for line in fh:
        if '#' not in line:
            fields = line.strip().split('\t')
            yield fields[9::]


def get_vcftuples_keep(fh, keepList):
    """yield  (chrom, pos, ref, alt, [ (sample,gt), ...., (sample,gt) ] ) but only for those samples in the keepList"""
    samples = get_vcfsamples(fh)
    for line in fh:
        if '#' not in line:
            fields = line.strip().split('\t')
            tuple = zip(samples,fields[9::])
            filtered_tuple= [x for x in tuple if x[0] in keepList]
            yield(fields[0], fields[1], fields[3], fields[4], filtered_tuple)




def get_vcftuple_passfilter(fh):
    """yield tuple with (chrom, pos, ref, alt, [ (sample,gt), ...., (sample,gt) ] )  of datalines that were not filtered"""
    samples = get_vcfsamples(fh)
    for line in fh:
        if '#' not in line:
            fields = line.strip().split('\t')
            if fields[6] == '.' or fields[6] == 'PASS':
                t=zip(samples,fields[9::])
                yield (fields[0],fields[1],fields[3], fields[4], t)




def get_vcftuples(fh):
    """yield tuple with (chrom, pos, ref, alt, [ (sample,gt), ..., (sample,gt) ] ) of datalines regarless of filter field"""
    samples = get_vcfsamples(fh)
    
    for line in fh:
        if '#' not in line:
            fields = line.strip().split('\t')
            t=zip(samples,fields[9::])
            yield (fields[0], fields[1], fields[3], fields[4], t)


def get_vcfsamples_keep(fh, keepList):
    """ yield a list of samples that are in the keepList """
    for line in fh:
        line.strip()
        if '#CHROM' in line:
            fields = line.strip().split('\t')
            samples = fields[9::]
            kept=[s for s in samples if s in keepList]
            return kept
     

def get_vcfsamples(fh):
    """return list of samples """
    for line in fh:
        if '#CHROM' in line:
            fields = line.strip().split('\t')
            samples = fields[9::]
            return samples
    
def getFormatfield(gt_string, index):
    """ return the ith field of a genotype format string specified by index"""
    gtformatfields=gt_string.split(':')
    return gtformatfields[index]

def stripGT(gt_string):
    """ strip a genotype string to only its GT field GT: => GT  """

    if gt_string == '.':
        return './.'
    else:
        gtformatfields=gt_string.split(':')
        return gtformatfields[0]

def yieldGQ (g1, formatstr):
    """ given a list of [ (sample, gt_string) ] parse out the GQ and yield the genotype quality for the genotype """
    formatfields=formatstr.split(':')
    if 'GQ' in formatfields:
        gq_index = formatfields.index('GQ')
    else:
        sys.stderr.write("GQ not in genotype format string, cannot return genotype quality!\n")
        exit(1)
    for (sample, gtstring) in g1:
        (allele1, allele2) = returnAlleles( stripGT(gtstring)  )
        if (allele1 == '.' or allele2 == '.'): continue
        gq=getFormatfield(gtstring, gq_index)
        yield gq
    


def returnAlleles(gt):
    """ return tuple (allele1, allele2) of stripped genotype GT field"""

    if ':' in gt:
        sys.stderr.write("strip string to only its GT field first!")
        return None
    elif '/' in gt:
        (a1,a2)= gt.split('/')
        return (a1,a2)
    elif '|' in gt:
        (a1, a2) = gt.split('|')
        return (a1,a2)
    elif '.' in gt:
        return ('.', '.')
    else:
        return None
        



def returnAlleles_unphased(gt):
    """ return tuple of (allele1, allele2) of unphased and stripped  GT string e.g. 0/0 returns (0,0) """
    if ':' in gt:
        sys.stderr.write("strip string to only its GT field first!")
        return None
    if '/' in gt:
        (a1,a2)= gt.split('/')
        return (a1,a2)
    else:
        return None

def returnAlleles_phased(gt):
    """ return tuple of (allele1, allele2) of phased and stripped  GT string e.g. 0|0 returns (0,0) """
    if ':' in gt:
        sys.stderr.write("strip string to only its GT field first!")
        return None
    if '|' in gt:
        (a1,a2)= gt.split('|')
        return (a1,a2)
    else:
        return None

def getCalledGenotypeCount (g):
    """ return the numebr of called genotpes given a list  [ (sample, gt_string) ... ] return the number of called genotypes """
    total=0
    for (sample, genotype) in g:
        (p1,p2) = returnAlleles ( stripGT( genotype ) )

        if typeofGenotype(p1,p2) != 3:
            total+=1

    return total


def calMaf (g):
    """ given a list in the form [ (sample, gt_string) ] determine the minor allele freq= alt/2*called genotypes """
    total=0
    alt_count=0
    for (sample, genotype) in g:
        (p1,p2) = returnAlleles ( stripGT( genotype ) )
        
        alleletype=typeofGenotype(p1,p2)

        if alleletype == 3:
            continue
        total+=1
        if alleletype ==1 : 
            alt_count+=1
        elif alleletype == 2:
            alt_count +=2
        else:
            pass
    maf = float(alt_count) / float(2*total)
    #print float(2*total)
    #print p1, p2, sample, maf, alt_count, total
    return  ( round(maf,3) )


def return_expected_genotype_counts ( g1 ):
    """ given a list in the form [ (sample, gt_string), ... ] return expected counts of 3 genotypes """

    theta = calMaf(g1)

    pAA =  float(1-theta) * float(1-theta) * float( getCalledGenotypeCount(g1) )
    pAB =  float( 2*theta*(1-theta) ) * float ( getCalledGenotypeCount(g1) )
    pBB =  float( theta) * float (theta) * float ( getCalledGenotypeCount(g1) )

    return (round(pAA,3), round(pAB,3), round(pBB,3) )


def return_observed_genotype_counts( g1 ):
    """ given a list in the form [ (sample, gt_string), ... ] return the counts of the 3 types of genotypes """
    homoz_ref=0
    het=0
    homoz_nonref=0
    for i in range (0, len( g1 ) ):
        (p1,p2) = returnAlleles ( stripGT( g1[i][1] ) )
        if typeofGenotype(p1,p2) == 3:
            continue
        else:
            alleletype=typeofGenotype(p1,p2)
            #print g1[i][0], alleletype, g1[i][1]
            if alleletype==0:
                homoz_ref+=1
            if alleletype==1:
                het+=1
            if alleletype==2:
                homoz_nonref+=1
    return( homoz_ref,  het, homoz_nonref )

def hweLRT ( g1 ):
    """ return the LRT for HWE """

    #2 * sum (obs * log(obs/exp)) = LRT

    lrt=0

    obs = (pAAobs, pABobs, pBBobs ) = return_observed_genotype_counts ( g1 )
    exp = (pAAexp, pABexp, pBBexp) =  return_expected_genotype_counts ( g1 )

    for i in range(0, len(obs)):
        if obs[i] == 0:
            lrt+=0
        else:
            lrt+= float(obs[i]) * math.log(obs[i]/exp[i])

    lrt *= 2
    return  lrt


def gq_calibration ( g1, truth, formatstr):
    """ given two lists in the form [ (sample, gt_string), ... ] and the second list are 'truth' genotypes return a list of [ (gq, 0|1) ] were 1
    indicates the genotype in g1 matched the truth and 0 indicates it did not. Note, 'GQ' needs to be in format string inorder to calibrate """

    formatfields = formatstr.split(':')
    if 'GQ' in formatstr:
        gq_index = formatfields.index('GQ')
        
    else:
        sys.stderr.write("genotype format string doesn't contain GQ, cannot perform genotype quality calibration!\n")
        exit(1)

    gq_calibration = []
    
    if len(g1) != len(truth):
        sys.stderr.write("gq_calibration: cannot compare  genotypes; unequal size of lists!\n")
        exit(1)
    else:
        for i in range(0, len(g1) ):
            if g1[i][0] != truth[i][0]:
                sys.stderr.write("gq_calibration:samples don't match in comparison\n")
                exit(1)
            (p1,p2) = returnAlleles ( stripGT( g1[i][1] ) )
            if typeofGenotype(p1,p2) == 3: continue
            
            gq=getFormatfield(g1[i][1], gq_index)
            
            (u1, u2) = returnAlleles ( stripGT( truth[i][1]) )
            if typeofGenotype(u1,u2) == 3: continue

            correctly_genotyped=0

            g1_alleletype=typeofGenotype(p1,p2)
            g2_alleletype=typeofGenotype(u1,u2)

            if g1_alleletype != 3 and g2_alleletype !=3 and ( g1_alleletype != g2_alleletype ):
                correctly_genotyped=0
            elif g1_alleletype != 3 and g2_alleletype !=3 and ( g1_alleletype == g2_alleletype ):
                correctly_genotyped=1
                #print g1_maxprob, g2[i][1], correctly_imputed
            else:
                pass

            gq_calibration.append( (gq, correctly_genotyped) )
            
    return gq_calibration

def posterior_imputed_gprob_calibration ( g1, g2, formatstr):
    """ given two lists of the form [ (sample, gt_string), ... ] and the second one is the 'truth' genotypes
        return a list [ (max_gprob, 0|1), ... ] where 1 is the imputed genotype matched the truth 0 is it did not
        note the list returned will only be as long as how many imputed gentypes there were for that marker  """

    formatfields=formatstr.split(':')
    if 'GPROB' in formatstr and 'OG' in formatstr:
        gprobs_index=formatfields.index('GPROB')
        og_index=formatfields.index('OG')
    else:
        sys.stderr.write("genotype format doesn't contain GPROB and/or OG , cannot perform posterior genotype calibration!\n")
        exit(1)
    imputed_genotype_calibration=[]

    if len(g1) != len(g2):
        sys.stderr.write("posterior_imputed_gprob_calibration: cannot compare  genotypes; unequal size of lists!\n")
    else:
        for i in range( 0, len( g1 ) ):

            if g1[i][0] != g2[i][0]:
                sys.stderr.write("samples don't match in genotypes comparison!")
                exit(1)

            #get the original genotype and make sure it was originally missing data (only want to compare imputed genotypes)
            (og1, og2) = returnAlleles ( getFormatfield( g1[i][1], og_index ) )
            if (typeofGenotype(og1,og2) != 3):
                continue

            #get the genotyep proba of the imputed genotype
            g1_maxprob= max ( getFormatfield( g1[i][1], gprobs_index ).split(';') )
            
            
            (p1,p2) = returnAlleles ( stripGT( g1[i][1] ) )
            (u1, u2) = returnAlleles ( stripGT( g2[i][1]) )

            

            g1_alleletype=typeofGenotype(p1,p2)
            g2_alleletype=typeofGenotype(u1,u2)
            correctly_imputed=0

            if g1_alleletype != 3 and g2_alleletype !=3 and ( g1_alleletype != g2_alleletype ):
                correctly_imputed=0
                #print  g1_maxprob, g1[i][1], g2[i][1], correctly_imputed
                
            elif g1_alleletype != 3 and g2_alleletype !=3 and ( g1_alleletype == g2_alleletype ):
                correctly_imputed=1
                #print g1_maxprob, g2[i][1], correctly_imputed
            else:
                pass

            imputed_genotype_calibration.append( (g1_maxprob, correctly_imputed) )
            

    return imputed_genotype_calibration


def compare_nonimputed_genotypes( g1, g2, formatstr, thresh):
    """ given two lists of the form [ (sample, gt_string), ... ] iterate and compare *non-imputed* genotypes in g1 to genotypes in g2  """
    """ to find out which genotypes in g1 are imputed the FORMAT string must contain the OG tag and the OG must be a nocall """
    """ return list of [ ( g1_alleletype, g2_alleletype, samplename) ] where alleletype is [0, homref; 1, het; 2 hom_nonref; 3, nocall  """

    formatfields=formatstr.split(':')
    if 'GPROB' in formatstr and 'OG' in formatstr:
        gprobs_index=formatfields.index('GPROB')
        og_index=formatfields.index('OG')
    else:
        sys.stderr.write("genotype format doesn't contain GPROB/OG, cannot determine if genotype was imputed!")
        exit(1)
    compared_imputed_genotypes=[]

    if len(g1) != len(g2):
        sys.stderr.write("cannot compare phased unphased genotypes; unequal size of lists!")
    else:
        for i in range( 0, len( g1 ) ):

            if g1[i][0] != g2[i][0]:
                sys.stderr.write("samples don't match in genotypes comparison!")
                exit(1)

            #check if the original genotype in g1 was imputed; if not pass on the comparison!
            (og1, og2) = returnAlleles ( getFormatfield( g1[i][1], og_index ) )
            if (typeofGenotype(og1,og2) == 3):
                continue
            g1_maxprob= max ( getFormatfield( g1[i][1], gprobs_index ).split(';') )
            g1_maxprob=float(g1_maxprob)
            #posterior prob of imputed genotype did not meet threshold!
            if g1_maxprob <=float(thresh):
                continue
            (p1,p2) = returnAlleles ( stripGT( g1[i][1] ) )
            (u1, u2) = returnAlleles ( stripGT( g2[i][1]) )

            g1_alleletype=typeofGenotype(p1,p2)
            g2_alleletype=typeofGenotype(u1,u2)
            #print g1[i], g1_alleletype
            #print g2[i], g2_alleletype
            compared_imputed_genotypes.append( (g1_alleletype, g2_alleletype, g1[i][0])  )
    return  compared_imputed_genotypes




def compare_imputed_genotypes( g1, g2, formatstr, thresh):
    """ given two lists of the form [ (sample, gt_string), ... ] iterate and compare *imputed* genotypes in g1 to genotypes in g2  """
    """ to find out which genotypes in g1 are imputed the FORMAT string must contain the OG tag and the OG must be a nocall """
    """ return list of [ ( g1_alleletype, g2_alleletype, samplename) ] where alleletype is [0, homref; 1, het; 2 hom_nonref; 3, nocall  """
    
    formatfields=formatstr.split(':')
    if 'GPROB' in formatstr and 'OG' in formatstr:
        gprobs_index=formatfields.index('GPROB')
        og_index=formatfields.index('OG')
    else:
        sys.stderr.write("genotype format doesn't contain GPROB/OG, cannot determine if genotype was imputed!")
        exit(1)
    compared_imputed_genotypes=[]

    if len(g1) != len(g2):
        sys.stderr.write("cannot compare phased unphased genotypes; unequal size of lists!")
    else:
        for i in range( 0, len( g1 ) ):

            if g1[i][0] != g2[i][0]:
                sys.stderr.write("samples don't match in genotypes comparison!")
                exit(1)

            #check if the original genotype in g1 was imputed; if not pass on the comparison!
            (og1, og2) = returnAlleles ( getFormatfield( g1[i][1], og_index ) )
            if (typeofGenotype(og1,og2) != 3):
                continue
            g1_maxprob= max ( getFormatfield( g1[i][1], gprobs_index ).split(';') )
            g1_maxprob=float(g1_maxprob)
            #posterior prob of imputed genotype did not meet threshold!
            if g1_maxprob <=float(thresh):
                continue
            (p1,p2) = returnAlleles ( stripGT( g1[i][1] ) )
            (u1, u2) = returnAlleles ( stripGT( g2[i][1]) )

            g1_alleletype=typeofGenotype(p1,p2)
            g2_alleletype=typeofGenotype(u1,u2)
            
            compared_imputed_genotypes.append( (g1_alleletype, g2_alleletype, g1[i][0])  )

    
    return  compared_imputed_genotypes

def compare_genotypes(g1, g2):
    """ given two lists of the form [ (sample, gt_string), ....  ], iterate thru and compare genotypes per individual    """
    """ return list of [ ( g1_alleletype, g2_alleletype, samplename) ] where alleletype is [0, homref; 1, het; 2 hom_nonref; 3, nocall """


    
    compared_genotypes=[] # [ (g1_type, g2_type, sample) ... ]

    if len(g1) != len(g2):
        sys.stderr.write("cannot compare phased unphased genotypes; unequal size of lists!")
    else:
        for i in range( 0, len( g1 ) ):
            if g1[i][0] != g2[i][0]: 
                sys.stderr.write("samples don't match in genotypes comparison!")
                exit(1)
            (p1,p2) = returnAlleles ( stripGT(g1[i][1] ) )
            (u1, u2) = returnAlleles ( stripGT(g2[i][1]) )

            g1_alleletype=typeofGenotype(p1,p2)
            g2_alleletype=typeofGenotype(u1,u2)
            #print g1[i], g1_alleletype
            #print g2[i], g2_alleletype
            compared_genotypes.append( (g1_alleletype, g2_alleletype, g1[i][0])  )
    return  compared_genotypes

def compare_phased_to_unphased(phased, unphased):
    """ give two lists of the form  [ (sample, gt_string), ....  ] in which one is unphased and theother phased iterate thru and compare the genotypes per individual  """

    unmatched_genotypes=[]
    if len(phased) != len(unphased):
        sys.stderr.write("cannot compare phased unphased genotypes; unequal size of lists!")
    else:
        for i in range( 0, len( phased ) ):
            if stripGT(unphased[i][1]) != '.':
                (p1,p2) = returnAlleles_phased ( stripGT(phased[i][1] ) )
                (u1, u2) = returnAlleles_unphased ( stripGT(unphased[i][1]) )
                if getNonRefDosage(p1,p2) != getNonRefDosage(u1,u2):
                    unmatched_genotypes.append( (phased[i][1], unphased[i][1], unphased[i][0] )    )
    return  unmatched_genotypes

    
def getNonRefDosage(allele1,allele2):
    """ return the number of non-ref alleles given the two alleles from a genotype """

    dosage=0
    if '1' in allele1: dosage+=1
    if '1' in allele2: dosage+=1

    return dosage

def typeofGenotype(allele1, allele2):
    """ I really should be a python version of a typedef here, but dont know how 
        hom_ref =1 het =2 hom_nonref=3 no_call=4 """

    if allele1 == '0' and allele2 == '0': return 0

    if  allele1 == '0' and allele2== '1': return 1
    if allele1 =='1' and allele2 == '0': return 1

    if allele1== '1' and allele2== '1': return 2

    if allele1== '.' or allele2 == '.': return 3
    if allele1 == '.' and allele2 == '.': return 3


def doPatternSearch(pattern, string):
    """ given a pattern object that represents a regular expression search for the pattern in string
        return matched string. if there is no match return None """

    match=pattern.search(string)

    if match.group() == None:
        return None
    else:
        return match.group()
    
    
def get_genotype_specific_snps(vcf_file,output_file_16,output_file_3,output_file_het):
    file_handle = utils_logging.open_input_file(vcf_file, pipe=False)
    output_fh_16 = utils_logging.open_output_file(output_file_16, pipe=False)
    output_fh_3 = utils_logging.open_output_file(output_file_3, pipe=False)
    output_fh_het = utils_logging.open_output_file(output_file_het, pipe=False)
    
    reader  = VcfReader(file_handle)
    all_samples = reader.get_sample_names()
    sample_16=[]
    sample_3=[]
    sample_4=[]
    for sample in all_samples:
        if sample.startswith('16'):
            sample_16.append(sample)
        if sample.startswith('4'):
            sample_4.append(sample)
        if sample.startswith('3'):
            sample_3.append(sample)
    number_discordant_16=0
    number_discordant_3=0
    number_discordant_4=0
    number_discordant=0
    nb_het_for_4=0
    nb_hom16_for_4=0
    nb_hom3_for_4=0
    nb_no_info_for_4=0
    nb_others_weird=0
    nb_no_interest=0
    for line in reader:
        stop=False
        #test if all the samples of the same type say the same things
        record = VcfRecord(line, reader.get_sample_names())
        all_genotypes_16 = record.get_valid_genotype_per_sample(genotype_quality_threshold=20, sample_list=sample_16)
        all_genotypes_4 = record.get_valid_genotype_per_sample(genotype_quality_threshold=20, sample_list=sample_4)
        all_genotypes_3 = record.get_valid_genotype_per_sample(genotype_quality_threshold=20, sample_list=sample_3)
        if len(all_genotypes_16)>1:
            number_discordant_16+=1
            stop =True
        if len(all_genotypes_4)>1:
            number_discordant_4+=1
            stop =True
        if len(all_genotypes_3)>1:
            number_discordant_3+=1
            stop = True
        if stop:
            number_discordant+=1
            continue
        if len(all_genotypes_16)>0 and len(all_genotypes_3)>0 and all_genotypes_16.keys()[0] != all_genotypes_3.keys()[0]:
            al16_1, al16_2 = returnAlleles(all_genotypes_16.keys()[0])
            al3_1, al3_2 = returnAlleles(all_genotypes_3.keys()[0])
            if al16_1==al16_2 and al3_1==al3_2:
                if len(all_genotypes_4)>0:
                    al4_1, al4_2 = returnAlleles(all_genotypes_4.keys()[0]) 
                    if al4_1==al4_2:
                        if al4_1==al16_1:
                            nb_hom16_for_4+=1
                            output_fh_16.write(line+'\n')
                        elif al4_1==al3_2:
                            output_fh_3.write(line+'\n')
                            nb_hom3_for_4+=1
                        else:
                            print line.strip()
                            nb_others_weird+=1
                    else:
                        output_fh_het.write(line+'\n')
                        nb_het_for_4+=1
                else:
                    nb_no_info_for_4+=1
                #print all_genotypes_16
                #print all_genotypes_3
                #print all_genotypes_4
                #print record
        else:
            nb_no_interest+=1
                
    output_fh_16.close()
    output_fh_3.close()
    output_fh_het.close()
    
    print 'number_discordant_16=%s'%number_discordant_16
    print 'number_discordant_3=%s'%number_discordant_3
    print 'number_discordant_4=%s'%number_discordant_4
    print 'number_discordant=%s'%number_discordant_4
    print 'nb_het_for_4=%s'%nb_het_for_4
    print 'nb_hom16_for_4=%s'%nb_hom16_for_4
    print 'nb_hom3_for_4=%s'%nb_hom3_for_4
    print 'nb_no_info_for_4=%s'%nb_no_info_for_4
    print 'nb_others_weird=%s'%nb_others_weird
    print 'nb_no_interest=%s'%nb_no_interest

def get_distance_between_samples(vcf_file):
    import random
    file_handle = utils_logging.open_input_file(vcf_file, pipe=False)
    reader  = VcfReader(file_handle)
    all_samples = reader.get_sample_names()
    sample_16=[]
    sample_3=[]
    sample_4=[]
    for sample in all_samples:
        if sample.startswith('16'):
            sample_16.append(sample)
        if sample.startswith('4'):
            sample_4.append(sample)
        if sample.startswith('3'):
            sample_3.append(sample)
    comp_16_3=0
    comp_16_4=0
    comp_16_R=0
    comp_3_4=0
    comp_3_R=0
    comp_4_R=0
    number_discordant_16=0
    number_discordant_3=0
    number_discordant_4=0
    number_discordant=0
    nb_interest=0
    nb_no_interest=0
    for line in reader:
        stop=False
        #test if all the samples of the same type say the same things
        record = VcfRecord(line, reader.get_sample_names())
        all_genotypes_16 = record.get_valid_genotype_per_sample(genotype_quality_threshold=20, sample_list=sample_16)
        all_genotypes_4 = record.get_valid_genotype_per_sample(genotype_quality_threshold=20, sample_list=sample_4)
        all_genotypes_3 = record.get_valid_genotype_per_sample(genotype_quality_threshold=20, sample_list=sample_3)
        if len(all_genotypes_16)>1:
            number_discordant_16+=1
            stop =True
        if len(all_genotypes_4)>1:
            number_discordant_4+=1
            stop =True
        if len(all_genotypes_3)>1:
            number_discordant_3+=1
            stop = True
        if stop:
            number_discordant+=1
            continue
        #compare 16 and 3
        if len(all_genotypes_16)>0 and len(all_genotypes_3)>0 and len(all_genotypes_4)>0:
            nb_interest+=1
            al16_1, al16_2 = returnAlleles(all_genotypes_16.keys()[0])
            al3_1, al3_2 = returnAlleles(all_genotypes_3.keys()[0])
            al4_1, al4_2 = returnAlleles(all_genotypes_4.keys()[0])
            
            if al16_1!=al16_2:
                al16 = random.choice([al16_1,al16_2])
            else:
                al16=al16_1
            if al3_1!=al3_2:
                al3 = random.choice([al3_1,al3_2])
            else:
                al3=al3_1
            if al4_1!=al4_2:
                al4 = random.choice([al4_1,al4_2])
            else:
                al4=al4_1
            if al16=='0':
                comp_16_R+=1
            if al3=='0':
                comp_3_R+=1
            if al4=='0':
                comp_4_R+=1
            if al16==al3:
                comp_16_3+=1
            if al16==al4:
                comp_16_4+=1
            if al3==al4:
                comp_3_4+=1
        else:
            nb_no_interest+=1
                
            
    print 'number_discordant_16=%s'%number_discordant_16
    print 'number_discordant_3=%s'%number_discordant_3
    print 'number_discordant_4=%s'%number_discordant_4
    print 'number_discordant=%s'%number_discordant
    print 'comp_16_R=%s'%comp_16_R
    print 'comp_3_R=%s'%comp_3_R
    print 'comp_4_R=%s'%comp_4_R
    print 'comp_16_3=%s'%comp_16_3
    print 'comp_16_4=%s'%comp_16_4
    print 'comp_3_4=%s'%comp_3_4
    print 'nb_interest=%s'%nb_interest
    print 'nb_no_interest=%s'%nb_no_interest
    comp_16_3=float(nb_interest)-float(comp_16_3)
    comp_16_4=float(nb_interest)-float(comp_16_4)
    comp_16_R=float(nb_interest)-float(comp_16_R)
    comp_3_4=float(nb_interest)-float(comp_3_4)
    comp_3_R=float(nb_interest)-float(comp_3_R)
    comp_4_R=float(nb_interest)-float(comp_4_R)
    sys.stdout.write('16        0.0 %.0f %.0f %.0f\n'%(comp_16_3,comp_16_4,comp_16_R))
    sys.stdout.write('3         %.0f 0.0 %.0f %.0f\n'%(comp_16_3,comp_3_4,comp_3_R))
    sys.stdout.write('4         %.0f %.0f 0.0 %.0f\n'%(comp_16_4,comp_3_4,comp_4_R))
    sys.stdout.write('R         %.0f %.0f %.0f 0.0\n'%(comp_16_R,comp_3_R,comp_4_R))

def get_treatment_specific_snps(vcf_file):
    file_handle = utils_logging.open_input_file(vcf_file, pipe=False)
    output_fh_16 = utils_logging.open_output_file(output_file_16, pipe=False)
    output_fh_3 = utils_logging.open_output_file(output_file_3, pipe=False)
    output_fh_het = utils_logging.open_output_file(output_file_het, pipe=False)
    
    reader  = VcfReader(file_handle)
    all_samples = reader.get_sample_names()
    sample_treatment_8=['168R1','168R2','168R3']
    sample_treatment_1=['161R1','161R2','161R3']
    number_discordant_treatment_8=0
    number_discordant_treatment_1=0
    number_discordant=0
    nb_het_for_4=0
    nb_hom16_for_4=0
    nb_hom3_for_4=0
    nb_no_info_for_4=0
    nb_others_weird=0
    nb_no_interest=0
    for line in reader:
        stop=False
        #test if all the samples of the same type say the same things
        record = VcfRecord(line, reader.get_sample_names())
        all_genotypes_treatment_8 = record.get_valid_genotype_per_sample(genotype_quality_threshold=20, sample_list=sample_treatment_8)
        all_genotypes_treatment_1 = record.get_valid_genotype_per_sample(genotype_quality_threshold=20, sample_list=sample_treatment_1)
        if len(all_genotypes_treatment_8)>1:
            number_discordant_treatment_8+=1
            stop =True
        if len(all_genotypes_treatment_1)>1:
            number_discordant_treatment_1+=1
            stop =True
        
        if stop:
            number_discordant+=1
            continue
    
        if len(all_genotypes_treatment_8)>0 and len(all_genotypes_treatment_1)>0 and all_genotypes_treatment_8.keys()[0] != all_genotypes_treatment_1.keys()[0]:
            al8_1, al8_2 = returnAlleles(all_genotypes_treatment_8.keys()[0])
            al1_1, al1_2 = returnAlleles(all_genotypes_treatment_1.keys()[0])
            if al16_1==al16_2 and al3_1==al3_2:
                if len(all_genotypes_4)>0:
                    al4_1, al4_2 = returnAlleles(all_genotypes_4.keys()[0]) 
                    if al4_1==al4_2:
                        if al4_1==al16_1:
                            nb_hom16_for_4+=1
                            output_fh_16.write(line+'\n')
                        elif al4_1==al3_2:
                            output_fh_3.write(line+'\n')
                            nb_hom3_for_4+=1
                        else:
                            print line.strip()
                            nb_others_weird+=1
                    else:
                        output_fh_het.write(line+'\n')
                        nb_het_for_4+=1
                else:
                    nb_no_info_for_4+=1
                #print all_genotypes_16
                #print all_genotypes_3
                #print all_genotypes_4
                #print record
        else:
            nb_no_interest+=1
                
    output_fh_16.close()
    output_fh_3.close()
    output_fh_het.close()
    
    print 'number_discordant_16=%s'%number_discordant_16
    print 'number_discordant_3=%s'%number_discordant_3
    print 'number_discordant_4=%s'%number_discordant_4
    print 'number_discordant=%s'%number_discordant_4
    print 'nb_het_for_4=%s'%nb_het_for_4
    print 'nb_hom16_for_4=%s'%nb_hom16_for_4
    print 'nb_hom3_for_4=%s'%nb_hom3_for_4
    print 'nb_no_info_for_4=%s'%nb_no_info_for_4
    print 'nb_others_weird=%s'%nb_others_weird
    print 'nb_no_interest=%s'%nb_no_interest

if __name__=="__main__":
    vcf_file='/home/tcezard/projects/2010053_Tom_Little_RNAseq_2/seanna_bam_files/final_bam_file/161_realigned_resorted_samtools_snp.vcf.gz'
    vcf_file='/home/tcezard/projects/2010053_Tom_Little_RNAseq_2/seanna_bam_files/final_bam_file/test.vcf.gz'
    vcf_file='/home/tcezard/projects/2010053_Tom_Little_RNAseq_2/seanna_bam_files/final_bam_file/all_samples_samtools_snp.vcf.gz'
    vcf_file='/home/tcezard/projects/2010053_Tom_Little_RNAseq_2/seanna_bam_files/final_bam_file/all_samples_GATK_snp.vcf_2.gz'
    output_file_16='/home/tcezard/projects/2010053_Tom_Little_RNAseq_2/seanna_bam_files/final_bam_file/all_samples_GATK_snp_16_diff_3_4_eq_16.vcf'
    output_file_3='/home/tcezard/projects/2010053_Tom_Little_RNAseq_2/seanna_bam_files/final_bam_file/all_samples_GATK_snp_16_diff_3_4_eq_3.vcf'
    output_file_het='/home/tcezard/projects/2010053_Tom_Little_RNAseq_2/seanna_bam_files/final_bam_file/all_samples_GATK_snp_16_diff_3_4_eq_het.vcf'
    get_genotype_specific_snps(vcf_file, output_file_16, output_file_3, output_file_het)
    
    
