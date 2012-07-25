'''
Created on 5 Mar 2010

@author: tcezard
'''
import logging
from utils import DNA_tools, utils_logging
import sys
from pprint import pprint



#############################################
#       Genes and transcript objects        #
#############################################

_all_genes={}
_all_transcripts={}
_all_exons={}

class Gene(object):
    """A gene is a structure that contains Transcripts."""
    def __init__(self, gene_id, reference, transcript_list=[]):
        self.gene_id=gene_id
        self.reference=reference
        self.transcript_list=transcript_list
        
    def __str__(self):
        return '%s\t%s\n'%(self.reference,self.gene_id)+ '\n'.join([str(t) for t in self.transcript_list])
    
    def to_string(self):
        out=['gene %s on %s with %d transcripts'%(self.gene_id,str(self.reference),len(self.transcript_list))]
        out.append('\n'.join([t.to_string() for t in self.transcript_list]))
        return '\n'.join(out)
    
    def _add_transcript(self,transcript):
        self.transcript_list.append(transcript)
    
class Transcript(object):
    """A Transcript is a structure that contains Exons."""
    def __init__(self, gene, transcript_id, strand, transcript_start=None, transcript_end=None, cds_start=None, cds_end=None, exon_list=[]):
        self.gene=gene
        self.transcript_id=transcript_id
        self.transcript_start=transcript_start
        self.transcript_end=transcript_end
        self.cds_start=cds_start
        self.cds_end=cds_end
        if DNA_tools.strand_is_positive(strand):
            self.strand='+'
            self._set_strand_specific_functions()
        else:
            self.strand='-'
            self._set_strand_specific_functions()
        self.exon_list=[]
        self._add_exons(exon_list)
        
    def __str__(self):
        return '\t'.join([str(self.transcript_id), str(self.transcript_start),
                   str(self.transcript_end), str(self.cds_start), str(self.cds_end), str(self.strand)])
    
    def __cmp__(self, other):
        return cmp(self.transcript_id,other.transcript_id)
    
    def to_string(self):
        out=['transcript %s with %d exons'%(self.transcript_id,len(self.exon_list))]
        #out.append('\n'.join([t.to_string() for t in self.transcript_list]))
        return '\n'.join(out)
    
    
    def get_number_of_exons(self):
        return len(self.exon_list)
    
    def get_exon_number(self,exon):
        try:
            return self.exon_list.index(exon)+1
        except ValueError:
            return None
    
    def get_position_in_mrna(self,position):
        pass
    
    def has_cds(self):
        if self.cds_start and self.cds_end:
            return True
        else:
            return False
    def get_position_in_cdna(self,position):
        pass

    def get_mrna_length(self):
        if not hasattr(self, "mrna_length"):
            self.mrna_length=self.get_mrna_length()
        return self.mrna_length
    
        
    def get_cdna_sequence(self, chr_sequence, complete_sequence=True, trim_sequence=False):
        """retrieve the cDNA sequence from a transcript and a chromosome sequence"""
        if hasattr(self, "cdna_sequence") and self.cdna_sequence:
            return self.cdna_sequence
        else:
            return self._get_cdna_sequence(chr_sequence,complete_sequence,trim_sequence)

    def _add_exon(self,exon):
        self.exon_list.append(exon)
        if DNA_tools.strand_is_positive(self.strand):
            self.exon_list.sort()
        else:
            self.exon_list.sort(reverse=True)
            
    def _add_exons(self,exon_list):
        self.exon_list.extend(exon_list)
        if DNA_tools.strand_is_positive(self.strand):
            self.exon_list.sort()
        else:
            self.exon_list.sort(reverse=True)

    def _get_mrna_length(self):
        length=0
        for exon in self.exon_list:
            length+=exon.exon_end-exon.exon_start+1
        return length
    
    def _set_strand_specific_functions(self):
        if self.strand=='+':
            self.get_position_in_mrna=self._get_position_in_mrna_pos_strand
            self.get_position_in_cdna=self._get_position_in_cdna_pos_strand
        else:
            self.get_position_in_mrna=self._get_position_in_mrna_neg_strand
            self.get_position_in_cdna=self._get_position_in_cdna_neg_strand

    def _get_position_in_mrna_pos_strand(self,position):
        """transform a position in the genome on a position in this mrna if the mrna is on the pos strand"""
        i=0
        over=False
        pos_in_mRNA=0
        while not over and i<len(self.exon_list):
            exon=self.exon_list[i]
            i+=1
            if position > exon.exon_end:
                pos_in_mRNA+=exon.exon_end-exon.exon_start+1
            elif position >= exon.exon_start:
                pos_in_mRNA+=(position-exon.exon_start)+1
                over=True
            else:
                logging.error("missuse of get_position_in_mrna function with position %s"%(position))
                logging.error("start=%s end=%s"%(exon.exon_start, exon.exon_end))
                raise Exception()
                pos_in_mRNA=0
        if not over: 
            logging.warning("position %s exons: %s"%(position, ', '.join(['%s-%s'%(exon.exon_start, exon.exon_end) for exon in self.exon_list])))
            return None
        else:
            return pos_in_mRNA
        
    def _get_position_in_mrna_neg_strand(self,position):
        """transform a position in the genome on a position in this mrna if the mrna is on the neg strand"""
        i=0
        over=False
        pos_in_mRNA=0
        while not over and i<len(self.exon_list):
            exon=self.exon_list[i]
            i+=1
            if position <= exon.exon_start:
                pos_in_mRNA+=exon.exon_end-exon.exon_start+1
            elif position < exon.exon_end:
                pos_in_mRNA+=(position-exon.exon_start)+1
                over=True
            else:
                logging.error("missuse of get_position_in_mrna function with position %s"%(position))
                logging.error("start=%s end=%s"%(exon.exon_start, exon.exon_end))
                raise Exception()
                pos_in_mRNA=0
        if not over: 
            logging.warning("position %s exons: %s"%(position, ', '.join(['%s-%s'%(exon.exon_start, exon.exon_end) for exon in self.exon_list])))
            return None
        else:
            return pos_in_mRNA
           
    def _get_position_in_cdna_pos_strand(self,position):
        """transform a position in the genome on a position in this cDNA if the cDNA is on the pos strand"""
        i=0
        over=False
        pos_in_cDNA=0
        while not over and i<len(self.exon_list):
            exon=self.exon_list[i]
            i+=1
            if exon.exon_cds_start is None:
                continue
            
            if position > exon.exon_cds_end:
                pos_in_cDNA+=exon.exon_cds_end-exon.exon_cds_start+1
            elif position >= exon.exon_cds_start:
                pos_in_cDNA+=(position-exon.exon_cds_start)+1
                over=True
            else:
                logging.error("missuse of get_position_in_cdna function with position %s"%(position))
                logging.error("%s"%(', '.join(['%s-%s'%(exon.exon_cds_start,exon.exon_cds_end) for exon in self.exon_list])))
                raise Exception()
        if not over: 
            logging.warning("position %s not  in cds of %s exons: %s"%(position,self.transcript_id,
                                                                       ', '.join(['%s-%s'%(exon.exon_cds_start, exon.exon_cds_end) for exon in self.exon_list])))
            return None
        else:
            return pos_in_cDNA
    
    def _get_position_in_cdna_neg_strand(self,position):
        """transform a position in the genome on a position in this cDNA if the cDNA is on the neg strand"""
        i=0
        over=False
        pos_in_cDNA=0
        while not over and i<len(self.exon_list):
            exon=self.exon_list[i]
            i+=1
            if exon.exon_cds_start is None:
                continue
            if position < exon.exon_cds_start:
                pos_in_cDNA+=exon.exon_cds_end-exon.exon_cds_start+1
            elif position <= exon.exon_cds_end:
                pos_in_cDNA+=(exon.exon_cds_end-position)+1
                over=True
            else:
                logging.error("missuse of get_position_in_cdna function with position %s"%(position))
                logging.error("%s"%(', '.join(['%s-%s'%(exon.exon_cds_start,exon.exon_cds_end) for exon in self.exon_list])))
                raise Exception()
        if not over: 
            logging.warning("position %s not  in cds of %s exons: %s"%(position,self.transcript_id,
                                                                       ', '.join(['%s-%s'%(exon.exon_cds_start, exon.exon_cds_end) for exon in self.exon_list])))
            return None
        else:
            return pos_in_cDNA
        
        
    def _get_cdna_sequence(self, chr_sequence, complete_sequence=True, trim_sequence=False):
        """retrieve the cDNA sequence from a transcript and a chromosome sequence"""
        cdna_seq=[]
        exons=self.exon_list[:]
        exons.sort()
        for exon in exons:
            if exon.exon_cds_start and exon.exon_cds_end:
                if int(exon.exon_cds_start)-1<len(chr_sequence) and int(exon.exon_cds_end)<=len(chr_sequence):
                    cdna_seq.append(chr_sequence[int(exon.exon_cds_start)-1:int(exon.exon_cds_end)])
                else:
                    logging.error('coding sequence %s-%s outside of chromosome boundaries length=%s'%(int(exon.exon_cds_start),int(exon.exon_cds_end),len(chr_sequence)))
                    return None
        full_sequence = ''.join(cdna_seq)
        if len(full_sequence)==0:
            return 
        
        #Test the sequence to see if its length is a multiple of 3
        if complete_sequence or trim_sequence:
            remaining=len(full_sequence)%3
            if not remaining==0:
                if complete_sequence:
                    to_add=3-remaining
                    if DNA_tools.strand_is_positive(self.strand):
                        #add remaining bases at the end of the sequence
                        if self.cds_end+to_add<=len(chr_sequence):
                            full_sequence=full_sequence+chr_sequence[self.cds_end-1:self.cds_end+to_add-1]
                        else:
                            full_sequence=full_sequence+('N'*to_add)
                    else:
                        #add remaining bases at the beginning of the sequence
                        if self.cds_start-to_add>0:
                            full_sequence=chr_sequence[self.cds_start-to_add-1:self.cds_start-1]+full_sequence
                        else:
                            full_sequence=('N'*to_add)+full_sequence
                    logging.warning('%s will have %d bases added cause its length=%s'%(self.transcript_id, to_add,len(full_sequence)-to_add))
                if trim_sequence:
                    if DNA_tools.strand_is_positive(self.strand):
                        full_sequence=full_sequence[:-remaining]
                    else:
                        full_sequence=full_sequence[remaining:]
                    logging.warning('%s will have %d bases trimmed of the end cause its length=%s now is %s'%(self.transcript_id, remaining,len(''.join(cdna_seq)),len(full_sequence)))
        if DNA_tools.strand_is_positive(self.strand):
            self.cdna_sequence=full_sequence.upper()
        else:
            self.cdna_sequence=DNA_tools.rev_complements(full_sequence).upper()
            
        # Test the sequence to see if it starts with an M
        AA, rest = translate_sequence(self.cdna_sequence[0:3])
        if AA != 'M':
            logging.warning('%s starts with %s instead of an M'%(self.transcript_id,AA))
        # Test the sequence to see if it ends with a STOP
        AA, rest = translate_sequence(self.cdna_sequence[-3:])
        if AA != '*':
            logging.warning('%s ends with %s instead of a STOP'%(self.transcript_id, AA))
        return self.cdna_sequence
    
    def _get_transcript_start(self):
        """calculate the transcript start from the exons"""
        if len(self.exon_list)>0:
            if DNA_tools.strand_is_positive(self.strand):
                return self.exon_list[0].exon_start
            else:
                return self.exon_list[-1].exon_start
        else:
            return None
        
    def _get_transcript_end(self):
        """calculate the transcript end from the exons"""
        if len(self.exon_list)>0:
            if DNA_tools.strand_is_positive(self.strand):
                return self.exon_list[-1].exon_end
            else:
                return self.exon_list[0].exon_end
        else:
            return None
    
    def _get_cds_start(self):
        """calculate the cds start from the exons"""
        if len(self.exon_list)>0:
            if DNA_tools.strand_is_positive(self.strand):
                for exon in self.exon_list:
                    if not exon.exon_cds_start is None:
                        return exon.exon_cds_start
            else:
                for exon in self.exon_list[::-1]:
                    if not exon.exon_cds_start is None:
                        return exon.exon_cds_start
        else:
            return None

    def _get_cds_end(self):
        """calculate the cds end from the exons"""
        if len(self.exon_list)>0:
            if DNA_tools.strand_is_positive(self.strand):
                for exon in self.exon_list[::-1]:
                    if not exon.exon_cds_end is None:
                        return exon.exon_cds_end
            else:
                for exon in self.exon_list:
                    if not exon.exon_cds_end is None:
                        return exon.exon_cds_end
        else:
            return None
    
    def _update_from_exons(self):
        if not self.transcript_start:
            self.transcript_start=self._get_transcript_start()
        if not self.transcript_end:
            self.transcript_end=self._get_transcript_end()
        if self.cds_start and self.cds_end:
            for exon in self.exon_list:
                exon._update_cds_info(self.cds_start, self.cds_end)
        if not self.cds_start:
            self.cds_start=self._get_cds_start()
        if not self.cds_end:
            self.cds_end=self._get_cds_end()
              
class Exon(object):
    """"""
    def __init__(self, transcript, exon_start, exon_end, exon_cds_start=None, exon_cds_end=None):
        self.transcript=transcript
        self.exon_start=exon_start
        self.exon_end=exon_end
        self.exon_cds_start=exon_cds_start
        self.exon_cds_end=exon_cds_end
        
    def __str__(self):
        return '\t'.join([str(self.transcript.transcript_id), str(self.get_exon_number()),
                          str(self.exon_start), str(self.exon_end), str(self.exon_cds_start), str(self.exon_cds_end) ])
    
    def __cmp__(self, other):
        return self.exon_start-other.exon_start
    
    def get_exon_number(self):
        return self.transcript.get_exon_number(self)
    
    def get_length(self):
        return self.exon_end-self.exon_start+1
    
    def _update_cds_info(self,cds_start,cds_end):
        if self.exon_start>cds_end or self.exon_end<cds_start:
            return
        if self.exon_start>cds_start:
            self.exon_cds_start=self.exon_start
        else:
            self.exon_cds_start=cds_start
        if self.exon_end<cds_end:
            self.exon_cds_end=self.exon_end
        else:
            self.exon_cds_end=cds_end

def get_gene(gene_id,reference,transcript_list=[]):
    """This function generate a gene from a list of transcript"""
    key='%s\t%s'%(reference,gene_id)
    if not _all_genes.has_key(key):
        gene = Gene(gene_id, reference, transcript_list)
        _all_genes[key]=gene
    else:
        gene = _all_genes.get(key)
    return gene

def get_transcript(reference, gene_id, transcript_id, strand, transcript_start=None, transcript_end=None, cds_start=None, cds_end=None, exon_list=[]):
    """This function generate a transcript with its list of exons"""
    key='%s\t%s'%(reference,transcript_id)
    if not _all_transcripts.has_key(key):
        gene = get_gene(gene_id, reference,transcript_list=[])
        transcript = Transcript(gene, transcript_id, strand, transcript_start, transcript_end, cds_start, cds_end, exon_list)
        gene._add_transcript(transcript)
        _all_transcripts[key]=transcript
    else:
        transcript = _all_transcripts.get(key)
    return transcript

def get_exon(reference, gene_id, transcript_id, strand, exon_start, exon_end, exon_cds_start=None, exon_cds_end=None, 
             transcript_start=None, transcript_end=None, cds_start=None, cds_end=None):
    """ """
    key = '%s\t%s\t%s\t%s'%(reference,transcript_id,exon_start, exon_end)
    if not _all_exons.has_key(key):
        exon_list=[]
        transcript  = get_transcript(reference=reference, gene_id=gene_id, transcript_id=transcript_id, 
                                     strand=strand, transcript_start=transcript_start, transcript_end=transcript_end,
                                     cds_start=cds_start, cds_end=cds_end, exon_list=exon_list)
        exon = Exon(transcript, exon_start, exon_end, exon_cds_start, exon_cds_end)
        transcript._add_exon(exon)
        _all_exons[key]=exon
    else:
        exon = _all_exons.get(key)
    return exon

def compare_exon_start(x, y):
    exon_start1=x[0]
    exon_start2=y[0]
    if int(exon_start1)>int(exon_start2): return 1
    elif int(exon_start1)<int(exon_start2): return -1
    else: return 0

def compare_exon_end(x, y):
    exon_end1=x[1]
    exon_end2=y[1]
    if int(exon_end1)>int(exon_end2): return 1
    elif int(exon_end1)<int(exon_end2): return -1
    else: return 0

def translate_sequence(sequence):
    AAs=[]
    for pos in range(0,len(sequence)-2,3):
        aa=DNA_tools.gencode.get(sequence[pos:pos+3].upper())
        if aa:
            AAs.append(aa)
        else:
            AAs.append("X")
    if len(sequence)%3==0:
        rest=''
    else:
        rest=sequence[-(len(sequence)%3):]
    return (''.join(AAs),rest)

def read_gff_file_old(gff_file):
    """Read gff file expects one feature per line in the following format gff3 format."""
    #II      Coding_transcript       CDS     3802    3984    .       +       1       ID=CDS:2L52.1;Note=Zinc finger%2C C2H2 type;Parent=Transcript:2L52.1;status=Partially_confirmed;wormpep=WP:CE32090
    nb_exon=0
    nb_transcript=0
    nb_cds=0
    nb_gene=0
    
    all_features_per_transcript={}
    all_gene_per_chr={}
    exonic_type=['CDS','exon','five_prime_UTR','three_prime_UTR','intron']
    transcript_type=['mRNA','ncRNA']
    gene_type=['gene']
    allowed_type=[]
    allowed_type.extend(exonic_type)
    allowed_type.extend(transcript_type)
    allowed_type.extend(gene_type)
    allowed_source=[]
    unused_type={}
    
    open_exon=utils_logging.open_input_file(gff_file, pipe=False)
    #allowed_source=['Coding_transcript']
    for line in open_exon:
        if line.startswith("#"):
            continue
        sp_line=line.strip().split('\t')
        if len(sp_line)<2:
            continue
        chr=sp_line[0]
        start=int(sp_line[3])
        end=int(sp_line[4])
        if sp_line[6]=='.' or DNA_tools.strand_is_positive(sp_line[6]):
            strand=1
        else:
            strand=-1
        source=sp_line[1]
        #skip thing that are not in the allowed source or not in the allowed type 
        if len(allowed_source)>0 and source not in allowed_source:
            continue
        type=sp_line[2]
        if type not in allowed_type:
            if unused_type.has_key(type):
                unused_type[type]+=1
            else:
                unused_type[type]=1
            continue
        
        attributes=sp_line[8].split(';')
        feature_id=None
        parents=[]
        name=None
        #Go through the attribute to find out what type of line it is
        for attribute in attributes:
            attribute_elements=attribute.split('=')
            type_attr=attribute_elements[0]
            attr_str="=".join(attribute_elements[1:])
            if type_attr=='ID':
                all_attr=attr_str.split(':')
                if len(all_attr)>1:
                    id_type=all_attr[0]
                    feature_id =all_attr[1]
                else:
                    feature_id =all_attr[0]
            elif type_attr=='Parent':
                for attr_elmt in attr_str.split(','):
                    all_attr = attr_elmt.split(':')
                    if len(all_attr)==2:
                        parents.append(all_attr)
                    elif len(all_attr)==1:
                        parents.append((None,all_attr[0]))
                    else:
                        logging.error("Unrecognized attribute %s"%(attr_elmt))
            elif type_attr=='Name':
                name=attr_str
            else:
                pass
        if type in exonic_type:
            for parent_type, parent_id in parents:
                if parent_type is None:
                    logging.debug("parent type can't be checked for %s."%parent_id)
                    logging.debug("line %s."%line.strip())
                elif parent_type!='Transcript' and parent_type!='mRNA':
                    logging.error('type %s has parent %s '%(type,parent_type))
                    logging.error("line %s."%line.strip())
                    continue
                
                parent_type='mRNA'
                info=all_features_per_transcript.get(parent_id)
                if info is None:
                    features_list=[]
                    all_features_per_transcript[parent_id]=(features_list,None,None,None,None,None,None)
                else:
                    features_list,dummy,dummy,dummy,dummy,dummy,dummy=info
                features_list.append((chr,start,end,strand,type,parent_id))
        elif type in transcript_type:
            if len(parents)>0:
                for parent_type, parent_id in parents:
                    if parent_type and parent_type.lower()!='gene':
                        logging.error('type %s has parent %s '%(type,parent_type))
                        logging.error("line %s."%line.strip())
                    else:
                        parent_type='gene'
                        info=all_features_per_transcript.get(feature_id)
                        if info is None:
                            all_features_per_transcript[feature_id]=([],chr,start,end,strand,feature_id,parent_id)
                        else:
                            features_list,old_chr,old_start,old_end,old_strand,transcript_id,gene_id=info
                            if (old_chr and old_chr!=chr) or (old_start and old_start!=start) or \
                            (old_end and old_end!=end) or (old_strand and old_strand!=strand) or \
                            (transcript_id and transcript_id!=feature_id) or(gene_id and gene_id!=parent_id):
                                logging.error('%s\t%s\t%s\t%s\t%s\t%s'%(old_chr,old_start,old_end,old_strand,transcript_id,gene_id))
                                logging.error('%s\t%s\t%s\t%s\t%s\t%s'%(chr,start,end,strand,feature_id,parent_id))
                            all_features_per_transcript[feature_id]=(features_list,chr,start,end,strand,feature_id,parent_id)
            else:
                parent_id=feature_id
                if not name is None:
                    parent_id = name
                info=all_features_per_transcript.get(feature_id)
                if info is None:
                    all_features_per_transcript[feature_id]=([],chr,start,end,strand,feature_id,parent_id)
                else:
                    features_list,old_chr,old_start,old_end,old_strand,transcript_id,gene_id=info
                    if (old_chr and old_chr!=chr) or (old_start and old_start!=start) or \
                    (old_end and old_end!=end) or (old_strand and old_strand!=strand) or \
                    (transcript_id and transcript_id!=feature_id) or(gene_id and gene_id!=parent_id):
                        logging.error('%s\t%s\t%s\t%s\t%s\t%s'%(old_chr,old_start,old_end,old_strand,transcript_id,gene_id))
                        logging.error('%s\t%s\t%s\t%s\t%s\t%s'%(chr,start,end,strand,feature_id,parent_id))
                    all_features_per_transcript[feature_id]=(features_list,chr,start,end,strand,feature_id,parent_id)
        elif type=='gene':
            pass
        else:
            logging.error('Unknown type %s'%(type))
    open_exon.close()
    for type in unused_type.keys():
        logging.warning('%d lines of type %s were unused'%(unused_type.get(type),type))
    all_genes_per_id={}
    #we've gone througt the whole file and gather exon for all transcript 
    #now create the structure that we'll pass on
    for transcript_id in all_features_per_transcript.keys():
        nb_transcript+=1
        all_exons=[]
        features_list,chr,transcript_start,transcript_end,strand,transcript_id,gene_id=all_features_per_transcript.get(transcript_id)
        if transcript_start is None or transcript_end is None or chr is None or strand is None:
            logging.error('transcript %s does not have valid information'%transcript_id)
            print features_list
            print chr,transcript_start,transcript_end,strand,transcript_id,gene_id
            continue
        all_genes_per_id[gene_id]=1
        if gene_id is None:
            gene_id=transcript_id
        features_list.sort(key=lambda info : info[1] )
        exon_start=0
        # take the first feature. Hopefully it's not an intron
        if len(features_list)==0:
            logging.error('transcript %s does not have features'%transcript_id)
            print features_list
            print chr,transcript_start,transcript_end,strand,transcript_id,gene_id
            continue
        feature_chr,curr_exon_start, curr_exon_end,strand,curr_type,dummy=features_list[0]
        if feature_chr!=chr:
            print features_list
            print chr,transcript_start,transcript_end,strand,transcript_id,gene_id
        if curr_type=='CDS':
            cds_start=curr_exon_start
            cds_end=curr_exon_end
        else:
            cds_start=None
            cds_end=None
        # merge with all the remaining pieces if they're not intron
        for feature in features_list[1:]:
            feature_chr,start,end,strand,type,parent_id=feature
            if feature_chr!=chr:
                print features_list
                print chr,transcript_start,transcript_end,strand,transcript_id,gene_id
            if type!='intron':
                if start<=curr_exon_end :
                    if end > curr_exon_end:
                        #logging.warning('Overlapping feature %s (%s-%s) and %s (%s-%s) for transcript %s'%(curr_type, curr_exon_start, curr_exon_end,
                        #                                                                                   type, start, end, transcript_id))
                        curr_exon_end=end
                else:
                    all_exons.append((curr_exon_start, curr_exon_end, transcript_id, gene_id, transcript_start,
                                     transcript_end, 0, 0, chr, strand))
                    curr_exon_start=start
                    curr_exon_end=end
                if type=='CDS':
                    if cds_start is None:
                        cds_start=start
                    if cds_end is None or cds_end<end:
                        cds_end=end
        #Add the last one to the array
        all_exons.append((curr_exon_start, curr_exon_end, transcript_id, gene_id, transcript_start,
                                     transcript_end, 0, 0, chr, strand)) 
        nb_exon+=len(all_exons)
        if cds_start is not None and cds_end is not None:
            nb_cds+=1
        #Propagate the cds info to all the exon and store them in a per chr storage
        for i in range(len(all_exons)):
            exon_start, exon_end, transcript_id, gene_id, transcript_start,transcript_end, dummy, dummy, chr, strand=all_exons[i]
            all_exons[i]=(exon_start, exon_end, transcript_id, gene_id, transcript_start,
                          transcript_end, cds_start, cds_end, chr, strand)
            exon = get_exon(reference=chr, gene_id=gene_id, transcript_id=transcript_id, strand=strand, exon_start=exon_start, exon_end=exon_end,
                            transcript_start=transcript_start, transcript_end=transcript_end,cds_start=cds_start, cds_end=cds_end)
            gene = exon.transcript.gene
            list_gene = all_gene_per_chr.get(chr)
            if list_gene is None:
                list_gene=set()
                all_gene_per_chr[chr]=list_gene
            
            list_gene.add(gene)
    nb_gene=len(all_genes_per_id)
    logging.info( 'retrieve %s genes %s transcripts with %s cds and %s exons'%(nb_gene,nb_transcript,nb_cds,nb_exon))
    return all_gene_per_chr

def read_gff_file(gff_file):
    """Read gff file expects one feature per line in the following format gff3 format."""
    #II      Coding_transcript       CDS     3802    3984    .       +       1       ID=CDS:2L52.1;Note=Zinc finger%2C C2H2 type;Parent=Transcript:2L52.1;status=Partially_confirmed;wormpep=WP:CE32090
    nb_exon=0
    nb_transcript=0
    nb_cds=0
    nb_gene=0
    increment=0
    all_gene_per_chr={}
    exonic_type=['CDS','exon','five_prime_UTR','three_prime_UTR','intron']
    transcript_type=['mRNA','ncRNA','transcript']
    gene_type=['gene']
    allowed_type=[]
    allowed_type.extend(exonic_type)
    allowed_type.extend(transcript_type)
    allowed_type.extend(gene_type)
    allowed_source=[]
    unused_type={}
    all_features_per_type={}
    all_features_per_type_count={}
    all_features={}
    open_exon=utils_logging.open_input_file(gff_file, pipe=False)
    for line in open_exon:
        if line.startswith("#"):
            continue
        sp_line=line.strip().split('\t')
        if len(sp_line)<2:
            continue
        chr=sp_line[0]
        start=int(sp_line[3])
        end=int(sp_line[4])
        if sp_line[6]=='.' or DNA_tools.strand_is_positive(sp_line[6]):
            strand=1
        else:
            strand=-1
        source=sp_line[1]
        #skip thing that are not in the allowed source or not in the allowed type 
        if len(allowed_source)>0 and source not in allowed_source:
            continue
        type=sp_line[2]
        if type not in allowed_type:
            if unused_type.has_key(type):
                unused_type[type]+=1
            else:
                unused_type[type]=1
            continue
        
        attributes=sp_line[8].split(';')
        
        #Go through the attribute to find out what type of line it is
        all_attributes={}
        for attribute in attributes:
            attribute_elements=attribute.split('=')
            type_attr=attribute_elements[0].strip()
            attr_str="=".join(attribute_elements[1:])
            all_attributes[type_attr]=attr_str
            
        children=[]
        if all_features.has_key(all_attributes.get('ID')):
            dummy_chr,dummy_type,dummy_start,dummy_end,dummy_strand,dummy_source,dummy_all_attributes,children=all_features.get(all_attributes.get('ID'))
            if not dummy_chr is None:
                logging.error("This feature has been seen more than once, create new id")
                logging.error(line.strip())
                all_attributes['ID']=all_attributes.get('ID')+"_%s"%(increment)
                increment+=1
        
        feature = (chr,type,start,end,strand,source,all_attributes,children)
        
        if not all_features_per_type_count.has_key(type):
            all_features_per_type_count[type]=1
        else:
            all_features_per_type_count[type]+=1
        
        if all_attributes.get('ID') is None and all_attributes.get('Parent') is None:
            logging.error("no id nor parent for this attribute")
            print line
            pprint(all_attributes)
        if all_attributes.has_key('ID'):
            all_features[all_attributes.get('ID')]=feature
            
        if all_attributes.get('Parent') is None:
            if all_features_per_type.has_key(type):
                all_features_per_type[type][all_attributes.get('ID')]=feature
            else:
                all_features_per_type[type]={all_attributes.get('ID'):feature}
        else:
            parent_feature = all_features.get(all_attributes.get('Parent'))
            if parent_feature is not None:
                (par_chr,par_type,par_start,par_end,par_strand,par_source,par_all_attributes,par_children) = parent_feature
                par_children.append(feature)
            else:
                all_features[all_attributes.get('Parent')]=(None,None,None,None,None,None,None,[feature])
    logging.info('load %s features'%(len(all_features)))
    for type in all_features_per_type_count.keys():
        logging.info('load %s %s features'%(all_features_per_type_count.get(type), type))
    
    def process_exons(gene_id, transcript_id, all_exonic_features):
        all_exonic_features.sort(key=lambda info : info[2] )
        exon_start=0
        # take the first feature. Hopefully it's not an intron
        if len(all_exonic_features)==0:
            logging.error('transcript %s does not have features'%transcript_id)
            return
        reference,curr_type,curr_start,curr_end,strand,source,all_attributes,children=all_exonic_features[0]
        if curr_type=='CDS':
            cds_start=curr_start
            cds_end=curr_end
        else:
            cds_start=None
            cds_end=None
        # merge with all the remaining pieces if they're not intron
        for feature in all_exonic_features[1:]:
            reference,type,start,end,strand,source,all_attributes,children=feature
            if type!='intron':
                if start<=curr_end+1 :
                    if end > curr_end:
                        curr_end=end
                else:
                    get_exon(reference=reference, gene_id=gene_id, transcript_id=transcript_id, strand=strand,
                             exon_start=curr_start, exon_end=curr_end, exon_cds_start=cds_start,
                             exon_cds_end=cds_end)
                    curr_start=start
                    curr_end=end
                    cds_start=None
                    cds_end=None
                if type=='CDS':
                    if cds_start is None:
                        cds_start=start
                    if cds_end is None or cds_end<end:
                        cds_end=end
        #Add the last one to the array
        exon=get_exon(reference=reference, gene_id=gene_id, transcript_id=transcript_id, strand=strand,
                 exon_start=curr_start, exon_end=curr_end, exon_cds_start=cds_start,
                 exon_cds_end=cds_end)
        exon.transcript._update_from_exons()
        return exon.transcript.gene
    
    def process_transcripts(gene_id, all_transcript_features):
        gene = None
        for transcript_feature in all_transcript_features:
            chr,curr_type,curr_start,curr_end,strand,source,all_attributes,children = transcript_feature
            if all_attributes.has_key('ID'):
                transcript_id = all_attributes.get('ID')
            else:
                transcript_id = gene_id
            if len(children)==0:
                exon=get_exon(reference=chr, gene_id=gene_id, transcript_id=transcript_id, strand=strand,
                             exon_start=curr_start, exon_end=curr_end, transcript_start=curr_start,
                             transcript_end=curr_end)
                gene = exon.transcript.gene
            else:
                gene = process_exons(gene_id, transcript_id, all_exonic_features=children)
        return gene
    
    def process_genes(feature):
        gene=None
        (chr,type,gene_start,gene_end,strand,source,all_attributes,children)=feature
        if len(children)==0:
            #create the transcript/exon from the gene
            exon=get_exon(reference=chr, gene_id=gene_id, transcript_id=gene_id, strand=strand,
                         exon_start=gene_start, exon_end=gene_end, transcript_start=gene_start,
                         transcript_end=gene_end)
            gene=exon.transcript.gene
        else:
            (dummy,child_type,dummy,dummy,dummy,dummy,dummy,dummy)=children[0]
            if child_type in transcript_type:
                gene = process_transcripts(gene_id, children)
            elif child_type in exonic_type:
                gene = process_exons(gene_id, gene_id, children)
        return gene
            
    #we've gone througt the whole file and gather all features
    #now create the structure that we'll pass on
    for gene_type_name in gene_type:
        all_genes_features = all_features_per_type.get(gene_type_name)
        
        if not all_genes_features is None:
            logging.debug("process %s features %s at the root"%(len(all_genes_features),gene_type_name))
            for gene_id in all_genes_features.keys():
                gene = process_genes(all_genes_features.get(gene_id))
                list_gene = all_gene_per_chr.get(gene.reference)
                if list_gene is None:
                    list_gene=[]
                    all_gene_per_chr[gene.reference]=list_gene
                list_gene.append(gene)
        else:
            logging.debug("No feature for %s at the root"%(gene_type_name))
    
    for transcript_type_name in transcript_type:
        all_transcript_features = all_features_per_type.get(transcript_type_name)
        if all_transcript_features:
            len_transcript=len(all_transcript_features)
        else:
            len_transcript=0
        logging.debug("process %s features %s at the root"%(len_transcript,transcript_type_name))
        if not all_transcript_features is None:
            for transcript_id in all_transcript_features.keys():
                gene = process_transcripts(transcript_id,[all_transcript_features.get(transcript_id)])
                list_gene = all_gene_per_chr.get(gene.reference)
                if list_gene is None:
                    list_gene=[]
                    all_gene_per_chr[gene.reference]=list_gene
                
                list_gene.append(gene)
    
     
    for chr in all_gene_per_chr.keys():
        list_gene = all_gene_per_chr.get(chr)
        for gene in list_gene:
            nb_gene+=1
            for transcript in gene.transcript_list:
                nb_transcript+=1
                if transcript.has_cds():
                    nb_cds+=1
                for exon in transcript.exon_list:
                    nb_exon+=1
    if len(unused_type)>0:
        logging.info("%s feature type unused in %s"%(len(unused_type), gff_file))
        logging.info(' - '.join(unused_type.keys()))
    logging.info( 'retrieve %s genes %s transcripts with %s cds and %s exons'%(nb_gene,nb_transcript,nb_cds,nb_exon))
    return all_gene_per_chr

def read_gtf_file(gtf_file):
    """Read gtf file expects one feature per line in the following format gff3 format."""
    #scaffold00002   Cufflinks       transcript      17153   18427   1000    +       .       gene_id "CUFF.7"; transcript_id "CUFF.7.1";
    nb_exon=0
    nb_transcript=0
    nb_cds=0
    nb_gene=0
    
    all_features_per_transcript={}
    all_gene_per_chr={}
    coding_type=['CDS','stop_codon','start_codon']
    exonic_type=['CDS','stop_codon','start_codon','exon','five_prime_UTR','three_prime_UTR','intron']
    transcript_type=['mRNA','ncRNA', 'transcript']
    gene_type=['gene']
    allowed_type=[]
    allowed_type.extend(exonic_type)
    allowed_type.extend(transcript_type)
    allowed_type.extend(gene_type)
    allowed_source=[]
    unused_type={}
    
    open_exon=utils_logging.open_input_file(gtf_file, pipe=False)
    #allowed_source=['Coding_transcript']
    for line in open_exon:
        if line.startswith("#"):
            continue
        sp_line=line.strip().split('\t')
        if len(sp_line)<2:
            continue
        reference=sp_line[0]
        start=int(sp_line[3])
        end=int(sp_line[4])
        if sp_line[6]=='.' or DNA_tools.strand_is_positive(sp_line[6]):
            strand=1
        else:
            strand=-1
        source=sp_line[1]
        #skip thing that are not in the allowed source or not in the allowed type 
        if len(allowed_source)>0 and source not in allowed_source:
            continue
        type=sp_line[2]
        if type not in allowed_type:
            if unused_type.has_key(type):
                unused_type[type]+=1
            else:
                unused_type[type]=1
            continue
        
        attributes=sp_line[8].split(';')
        gene_id=None
        transcript_id=None
        name=None
        #Go through the attribute to get gene_id and transcript_id
        feature={'reference':reference,'ft_start':start,'ft_end':end,'ft_strand':strand,
                'ft_type':type}
        for attribute in attributes:
            attribute = attribute.strip()
            attribute_elements=attribute.split(" ")
            type_attr=attribute_elements[0]
            attr_str=" ".join(attribute_elements[1:])
            attr_str=attr_str.replace('"','')
            feature[type_attr]=attr_str.strip()
            
        if not feature.has_key('gene_id') or not feature.has_key('transcript_id') :
            logging.error("line %s is missing gene_id or transcript_id"%line.strip())
            continue
        
        if type in exonic_type:
            info=all_features_per_transcript.get(feature['transcript_id'])
            if info is None:
                features_list=[]
                all_features_per_transcript[feature['transcript_id']]=(features_list,None,None,None,None,feature['transcript_id'],feature['gene_id'])
            else:
                features_list,dummy,dummy,dummy,dummy,dummy,dummy=info
            features_list.append(feature)
            
        elif type in transcript_type:
            info=all_features_per_transcript.get(feature['transcript_id'])
            if info is None:
                all_features_per_transcript[feature['transcript_id']]=([],reference,start,end,strand,feature['transcript_id'],feature['gene_id'])
            else:
                features_list,old_chr,old_start,old_end,old_strand,old_transcript_id,old_gene_id=info
                all_features_per_transcript[transcript_id]=(features_list,reference,start,end,strand,feature['transcript_id'],feature['gene_id'])
        elif type=='gene':
            #not bothering with gene features yet
            pass
        else:
            logging.error('Unknown type %s'%(type))
    open_exon.close()
    for type in unused_type.keys():
        logging.warning('%d lines of type %s were unused'%(unused_type.get(type),type))
    all_genes_per_id={}
    #we've gone through the whole file and gather features for all transcripts
    #now create the structure that we'll pass on
    for transcript_id in all_features_per_transcript.keys():
        nb_transcript+=1
        exon_start=None
        exon_end=None
        exon_cds_start = None
        exon_cds_end = None
        
        features_list, reference, transcript_start, transcript_end, strand, transcript_id, gene_id = all_features_per_transcript.get(transcript_id)
        if len(features_list)>0:
            features_list.sort(key=lambda info : info['ft_start'] )
            # take the first feature. Hopefully it's not an intron
            feature = features_list[0]
            if not reference:
                reference=feature['reference']
            if not gene_id:
                gene_id=feature['gene_id']
            if not strand:
                strand=feature['ft_strand']
            
            exon_start=feature['ft_start']
            exon_end=feature['ft_end']
            if feature['ft_type'] in coding_type:
                exon_cds_start = feature['ft_start']
                exon_cds_end = feature['ft_end']
            else:
                exon_cds_start = None
                exon_cds_end = None
        transcript=get_transcript(reference=reference, gene_id=gene_id, transcript_id=transcript_id,
                                  strand=strand, transcript_start=transcript_start, transcript_end=transcript_end)
        
        if len(features_list)>1:
            # merge with all the remaining pieces if they're not intron
            for feature in features_list[1:]:
                if reference and feature['reference']!=reference:
                    logging.error("feature has a different reference than parent %s -- %s "%(feature['reference'],reference))
                    logging.error(str(feature))
                    logging.error(str((reference,transcript_start,transcript_end,strand,transcript_id,gene_id)))
                if type!='intron':
                    if feature['ft_start'] <= exon_end :
                        if feature['ft_end'] > exon_end:
                            #logging.warning('Overlapping feature %s (%s-%s) and %s (%s-%s) for transcript %s'%(curr_type, curr_exon_start, curr_exon_end,
                            #                                                                                   type, start, end, transcript_id))
                            exon_end = feature['ft_end']
                    else:
                        exon=get_exon(reference=reference, gene_id=gene_id, transcript_id=transcript_id,strand=strand, 
                                      exon_start=exon_start, exon_end=exon_end, exon_cds_start=exon_cds_start,
                                      exon_cds_end=exon_cds_end, transcript_start=transcript_start, transcript_end=transcript_end)
                        exon_cds_start=None
                        exon_cds_end=None
                        exon_start=feature['ft_start']
                        exon_end=feature['ft_end']
                        nb_exon+=1
                        
                    if feature['ft_type'] in coding_type:
                        if exon_cds_start is None:
                            exon_cds_start=feature['ft_start']
                        if exon_cds_end is None or exon_cds_end<feature['ft_end']:
                            exon_cds_end=feature['ft_end']
        #Add the last one to the array
        exon=get_exon(reference=reference, gene_id=gene_id, transcript_id=transcript_id, strand=strand,
                      exon_start=exon_start, exon_end=exon_end, exon_cds_start=exon_cds_start,
                      exon_cds_end=exon_cds_end, transcript_start=transcript_start, transcript_end=transcript_end)
        nb_exon+=1
        transcript=exon.transcript
        transcript._update_from_exons()
        if transcript.cds_start:
            nb_cds+=1
        gene = transcript.gene
        all_genes_per_id[gene]=1
        list_gene = all_gene_per_chr.get(reference)
        if list_gene is None:
            list_gene=set()
            all_gene_per_chr[reference]=list_gene
        list_gene.add(gene)
    nb_gene=len(all_genes_per_id)
    logging.info( 'retrieve %s genes %s transcripts with %s cds and %s exons'%(nb_gene,nb_transcript,nb_cds,nb_exon))
    return all_gene_per_chr

def read_bed_file(bed_file):
    """Read bed file expects transcript per line in the following format:
chromosome, start, end start end score strand frame group 
Only the chromosome, start, end, name, strand and group are used, the other are ignored."""
    all_exons_per_chr={}
    all_gene_per_chr={}
    all_gene_id_seen={}
    record_number=1
    nb_transcript=0
    nb_cds=0
    nb_exon=0
    #chr1    16145919        16146865        ASO3599 0       +       16145919        16146865        0       1       946     0
    #chr1    16660029        16667563        ASO1842 0       +       16660029        16667563        0       2       2340,905        0,6629
    open_exon=utils_logging.open_input_file(bed_file, pipe=False)
    line = open_exon.readline()
    while line.startswith("#") or line.startswith("browser") or line.startswith("track"):
        line = open_exon.readline()
    while line:
        sp_line=line.strip().split('\t')
        if len(sp_line)>3:
            gene_id=sp_line[3]
            time_seen=all_gene_id_seen.get(gene_id)
            if not time_seen:
                time_seen=0
            time_seen+=1
            all_gene_id_seen[gene_id]=time_seen
            transcript_id=sp_line[3]+"_%s"%time_seen
        else:
            gene_id='record_%s'%record_number
            transcript_id='record_%s'%record_number
            record_number+=1
        nb_transcript+=1
        chr=sp_line[0]
        #Assuming bed file from UCSC so 
        transcript_start=int(sp_line[1])+1
        transcript_end=int(sp_line[2])
        
        if len(sp_line)>5 and not DNA_tools.strand_is_positive(sp_line[5]):
            strand=-1
        else:
            strand=1
        if len(sp_line)>7:
            cds_start=int(sp_line[6])+1
            cds_end=int(sp_line[7])
            nb_cds+=1
        else:
            cds_start=None
            cds_end=None
        if len(sp_line)>11:
            if sp_line[10].endswith(','):
                block_sizes=sp_line[10][:-1].split(',')
            else:
                block_sizes=sp_line[10].split(',')
            if sp_line[11].endswith(','):
                block_starts=sp_line[11][:-1].split(',')
            else:
                block_starts=sp_line[11].split(',')
        else:
            block_sizes=[transcript_end-transcript_start+1]
            block_starts=[0]
        for i,start in enumerate(block_starts):
            nb_exon+=1
            exon_start=transcript_start+int(start)
            exon_end=transcript_start+int(start)+int(block_sizes[i])-1
            exon=get_exon(reference=chr, gene_id=gene_id, transcript_id=transcript_id, strand=strand, 
                          exon_start=exon_start, exon_end=exon_end, transcript_start=transcript_start,
                          transcript_end=transcript_end, cds_start=cds_start, cds_end=cds_end)
            gene = exon.transcript.gene
            list_gene = all_gene_per_chr.get(chr)
            if list_gene is None:
                list_gene=set()
                all_gene_per_chr[chr]=list_gene
            list_gene.add(gene)
        line = open_exon.readline()
    nb_gene=0
    logging.info( 'retrieve %s genes %s transcripts with %s cds and %s exons'%(nb_gene,nb_transcript,nb_cds,nb_exon))
    return all_exons_per_chr

def read_ucsc_file_not_knowngene(ucsc_file):
    return read_ucsc_file(ucsc_file,knowngene=False)

def read_ucsc_file(ucsc_file,knowngene=True):
    """Read ucsc file expects transcript per line in the following format:
tanscript id chromosome, strand, transcript_start, transcript_end cds_start cds_end number_exon, 
coma_separated_exons_start, coma_separated_exons_end proteinID, and alignID.
The last two column are not used.
The coordinate are assumed to be 0-based"""
    all_genes_per_chr={}
    all_gene_id_seen={}
    record_number=1
    #uc007aet.1      chr1    -       3195984 3205713 3195984 3195984 2       3195984,3203519,        3197398,3205713,                uc007aet.1
    #uc007aeu.1      chr1    -       3204562 3661579 3206102 3661429 3       3204562,3411782,3660632,        3207049,3411982,3661579,        Q5GH67  uc007aeu.1
    #uc007aev.1      chr1    -       3638391 3648985 3638391 3638391 2       3638391,3648927,        3640590,3648985,                uc007aev.1
    nb_transcript=0
    nb_cds=0
    nb_exon=0
    all_genes=set()
    open_exon=utils_logging.open_input_file(ucsc_file, pipe=False)
    line = open_exon.readline()
    all_transcript_id={}
    while line.startswith("#") or line.startswith("browser") or line.startswith("track"):
        line = open_exon.readline()
    while line:
        nb_transcript+=1
        sp_line=line.strip().split('\t')
        if not knowngene:
            sp_line=sp_line[1:]
        if len(sp_line)>=12:
            gene_id=sp_line[11]
            all_genes.add(gene_id)
        else:
            gene_id=sp_line[0]
            all_genes.add(gene_id)
        
        #The transcript of UCSC databases are not always unique
        #Add a dot plus a digit as it's done with refseq
        transcript_id=sp_line[0]
        i=0
        while all_transcript_id.get(transcript_id) is not None:
            i+=1
            transcript_id=sp_line[0]+'.%s'%i
        all_transcript_id[transcript_id]=1
        
        chr=sp_line[1]
        if DNA_tools.strand_is_positive(sp_line[2]):
            strand=1
        else:
            strand=-1
        transcript_start=int(sp_line[3])+1
        transcript_end=int(sp_line[4])
        #1354    ENST00000444633    chr10    +    100807395    100808076    100808076    100808076    1    100807395,    100808076,    0    ENSG00000234109    none    none    -1,
        no_cds=False
        if len(sp_line)>=14:
            if sp_line[13]=='none' or sp_line[14]=='none' or sp_line[13]=='unk' or sp_line[14]=='unk':
                no_cds=True
        if sp_line[5].isdigit() and not no_cds:
            nb_cds+=1
            cds_start=int(sp_line[5])+1
        else:
            cds_start=None
        if sp_line[6].isdigit() and not no_cds:
            cds_end=int(sp_line[6])
        else:
            cds_end=None
        if sp_line[8].endswith(','):
            exons_starts=sp_line[8][:-1].split(',')
        else:
            exons_starts=sp_line[8].split(',')
        if sp_line[8].endswith(','):
            exons_ends=sp_line[9][:-1].split(',')
        else:
            exons_ends=sp_line[9].split(',')
        transcript = get_transcript(reference=chr, gene_id=gene_id, transcript_id=transcript_id, strand=strand,
                                    transcript_start=transcript_start, transcript_end=transcript_end, cds_start=cds_start, cds_end=cds_end)
        for i,exon_start in enumerate(exons_starts):
            nb_exon+=1
            exon_start=int(exon_start)+1
            exon_end=int(exons_ends[i])
            list_gene=all_genes_per_chr.get(chr)
            if not list_gene:
                list_gene=set()
                all_genes_per_chr[chr]=list_gene
            get_exon(reference=chr, gene_id=gene_id, transcript_id=transcript_id, strand=strand, 
                     exon_start=exon_start, exon_end=exon_end)
        transcript._update_from_exons()
        gene = transcript.gene
        list_gene.add(gene)
        line = open_exon.readline()
    
    logging.info( 'retrieve %s genes %s transcripts with %s cds and %s exons'%(len(all_genes), nb_transcript, nb_cds,nb_exon))
    return all_genes_per_chr

def recognize_format(input_file):
    """Recognize the format of a file in between gff, bed, and exon."""
    open_file=utils_logging.open_input_file(input_file, pipe=False)
    line = open_file.readline()
    while line and (line.startswith("#") or line.startswith("browser") or line.startswith("track")):
        line = open_file.readline()
    #chr1    16145919        16146865        ASO3599 0       +       16145919        16146865        0       1       946     0
    # bed has 3 mandatory field
    open_file.close()
    sp_line = line.strip().split('\t')
    if len(sp_line)>=3 and sp_line[1].isdigit() and sp_line[2].isdigit() and int(sp_line[1])<int(sp_line[2]):
        #probably bed
        is_bed=True
        if len(sp_line)>=6 and sp_line[5] not in ['+', '-', '.']:
            # not bed
            is_bed=False
        if len(sp_line)>=12 :
            if sp_line[10].endswith(','):
                block_sizes=sp_line[10][:-1].split(',')
            else:
                block_sizes=sp_line[10].split(',')
            if sp_line[11].endswith(','):
                block_starts=sp_line[11][:-1].split(',')
            else:
                block_starts=sp_line[11].split(',')
            if len(block_sizes)==len(block_starts):
                for block_size in block_sizes:
                    if not block_size.isdigit():
                        is_bed=False
                for block_start in block_starts:
                    if not block_start.isdigit():
                        is_bed=False
            else:
                is_bed=False
        if is_bed:
            logging.info("Recognize bed file for %s"%input_file)
            return read_bed_file
    
    #chr22  TeleGene    promoter  10010000  10010100  900    +  .  TG1
    # gff has 9 mandatory field
    if len(sp_line)==9 and sp_line[3].isdigit() and sp_line[4].isdigit() and int(sp_line[3])<int(sp_line[4]) and sp_line[6] in ['+','-','.']:
        
        #probably gff
        #check if gtf 
        if sp_line[8].find('gene_id')>=0 and sp_line[8].find('transcript_id')>=0:
            logging.info("Recognize gtf file for %s"%input_file)
            return read_gtf_file
        #Now I'm sure it's gff
        logging.info("Recognize gff file for %s"%input_file)
        return read_gff_file
    
    
    #uc007aeu.1      chr1    -       3204562 3661579 3206102 3661429 3       3204562,3411782,3660632,        3207049,3411982,3661579,        Q5GH67  uc007aeu.1
    # ucsc has  at least 10 mandatory field
    if len(sp_line)>9 and sp_line[2] in ['-','+'] and sp_line[3].isdigit() and sp_line[4].isdigit() and int(sp_line[3])<int(sp_line[4]):
        #probably ucsc
        logging.info("Recognize ucsc file for %s"%input_file)
        return read_ucsc_file
    elif len(sp_line)>10 and sp_line[3] in ['-','+'] and sp_line[4].isdigit() and sp_line[5].isdigit() and int(sp_line[4])<int(sp_line[5]):
        #probably ucsc
        logging.info("Recognize ucsc (not known gene) file for %s"%input_file)
        return read_ucsc_file_not_knowngene

        
    logging.error("Unrecognize format for file %s"%input_file)
    return None

class Annotation_Retriver():
    """This class can read from a from a file.
    It's used in Base_annotator to have a common way getting the exon information.
    database structure is defined in gene_annotation and exon file format is defined in read_exon_file's doc 
    """
    def __init__(self, annotation_file=None):
        self.all_gene_per_chr=None
        self.source=None
        if not annotation_file:
            logging.error("No annotation resource given")
            return 
        elif annotation_file:
            self.source=annotation_file
            reader=recognize_format(annotation_file)
            if reader is not None:
                self.all_gene_per_chr=reader(annotation_file)
    
    def get_annotation_from_chr(self,chr):
        if self.all_gene_per_chr:
            return self.all_gene_per_chr.get(chr)
        
    def get_chr_names(self):
        if self.all_gene_per_chr:
            return self.all_gene_per_chr.keys()
    
    def get_source(self):
        return self.source
    
class Exon_annotation_Retriver():
    """This class can read exon information from a database or from a file.
    It's used in Base_annotator to have a common way getting the exon information.
    database structure is defined in gene_annotation and exon file format is defined in read_exon_file's doc 
    """
    def __init__(self, annotation_file=None, annotation_database=None):
        self.anno_database_retriver=None
        if not annotation_file:
            logging.error("No annotation resource given")
            return 
        elif annotation_database:
            self.source=annotation_database
            logging.critical('The functionality of Exon_annotation_Retriver for databases is deprecated and will not work anymore')
            sys.exit(0)
        elif annotation_file:
            self.annotation_reader = Annotation_Retriver(annotation_file)
    
    def get_annotation_from_chr(self,chr):
        if self.annotation_reader:
            all_exons=[]
            all_genes= self.annotation_reader.get_annotation_from_chr(chr)
            if all_genes:
                for gene in all_genes:
                    for transcript in gene.transcript_list:
                        for exon in transcript.exon_list:
                            all_exons.append((exon.exon_start, exon.exon_end, transcript.transcript_id,
                                              gene.gene_id, transcript.transcript_start, transcript.transcript_end,
                                              transcript.cds_start, transcript.cds_end, gene.reference, transcript.strand))
            return all_exons
        
    def get_chr_names(self):
        if self.annotation_reader:
            return self.annotation_reader.get_chr_names()
    
    def get_source(self):
        return self.annotation_reader.get_source()

if __name__=="1__main__":
    utils_logging.init_logging()
    annotation_file='/home/tcezard/genomes/homo_sapiens/hg19/annotations/refGene_2011_05_27.txt.gz'
    #annotation_file='/home/tcezard/projects/2010053_Tom_Little_RNAseq_2/reference/daphmagna_201104m8.pasaupdate.gff.gz'
    #annotation_file='/home/tcezard/projects/2010053_Tom_Little_RNAseq_2/for tim/dmag_ep24augmap2an2.gff'
    annotation_reader = Annotation_Retriver(annotation_file=annotation_file)
    gene_list = annotation_reader.get_annotation_from_chr('chr1')
    #annotations = annotation_reader.get_annotation_from_chr('scaffold00002')
    for gene in gene_list:
        for transcript in gene.transcript_list:
            print transcript

if __name__=="1__main__":
    from utils.GenomeLoader import GenomeLoader
    utils_logging.init_logging()
    annotation_file='/home/tcezard/test_annotation.txt'
    genome_file='/home/tcezard/genomes/homo_sapiens/hg19/hg19.fa'
    genome_loader = GenomeLoader(genome_file)
    
    annotation_reader = Annotation_Retriver(annotation_file=annotation_file)
    
    gene_list = annotation_reader.get_annotation_from_chr('chr1')
    header, chr1_sequence = genome_loader.get_chr('chr1')
    #annotations = annotation_reader.get_annotation_from_chr('scaffold00002')
    for gene in gene_list:
        for transcript in gene.transcript_list:
            print transcript.transcript_id
            cdna_sequence = transcript.get_cdna_sequence(chr1_sequence)
            pos_in_cdna = transcript.get_position_in_cdna(894618)
            if  pos_in_cdna:
                print pos_in_cdna, cdna_sequence[pos_in_cdna-1]
            pos_in_cdna = transcript.get_position_in_cdna(894617)
            if  pos_in_cdna:
                print pos_in_cdna, cdna_sequence[pos_in_cdna-1]
            pos_in_cdna = transcript.get_position_in_cdna(894616)
            if  pos_in_cdna:
                print pos_in_cdna,cdna_sequence[pos_in_cdna-1]

if __name__=="1__main__":
    utils_logging.init_logging()
    gene_id='NR_024540'
    reference='chr1'
    transcript_id='NR_024540'
    transcript_start='14361'
    transcript_end='29370'
    cds_start='29370'
    cds_end='29370'
    strand='-'
    exon_number='1'
    exon_start='29320'
    exon_end='29370'
    
    print get_exon(gene_id, reference, transcript_id, transcript_start, transcript_end, cds_start, cds_end, strand, exon_start, exon_end)

if __name__=="__main__":
    utils_logging.init_logging()
    sys.argv[1]
    read_gff_file(sys.argv[1])
    read_gff_file_old(sys.argv[1])
