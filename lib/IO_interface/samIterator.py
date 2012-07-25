"""
@author:  Yaron S. Butterfield
@author: Timothee Cezard
 
 Create a SAM alignment format iterator
 
"""
import logging, re

from utils import utils_logging


class Sam_record(object):  
    def __init__(self, sam_line):
        self._sp_line=sam_line.strip().split()
    
    def get_query_name(self): return self._sp_line[0]
    def get_flag(self): return int(self._sp_line[1])
    def get_reference_name(self): return self._sp_line[2]
    def get_position(self): return int(self._sp_line[3])
    def get_mapping_quality(self): return int(self._sp_line[4])
    def get_cigar_string(self): return self._sp_line[5]
    def get_mate_reference_name(self): return self._sp_line[6]
    def get_mate_pos(self): return int(self._sp_line[7])
    def get_insert_size(self): return int(self._sp_line[8])
    def get_query_sequence(self): return self._sp_line[9]
    def get_query_quality(self): return self._sp_line[10]
    def get_tag(self,tag):
        res=None
        if len(self._sp_line)>10:
            for flag in self._sp_line[11:]:
                if flag.startswith(tag):
                    sp_flag=flag.split(':',2)
                    if sp_flag[1]=="i":
                        res=int(sp_flag[2])
                    elif sp_flag[1]=="A" or sp_flag[1]=="Z":
                        res=sp_flag[2]
                    elif sp_flag[1]=="f":
                        res=float(sp_flag[2])
                    else:
                        logging.error("Unsupported format for %s"%flag)
        return res
    
    def get_all_tags(self):
        res={}
        if len(self._sp_line)>10:
            for flag in self._sp_line[11:]:
                sp_flag=flag.split(':',2)
                if sp_flag[1]=="i":
                    res[sp_flag[0]]=int(sp_flag[2])
                elif sp_flag[1]=="A" or sp_flag[1]=="Z":
                    res[sp_flag[0]]=sp_flag[2]
                elif sp_flag[1]=="f":
                    res[sp_flag[0]]=float(sp_flag[2])
                else:
                    logging.error("Unsupported format for %s"%flag)
        return res
    
    def get_read_length(self): return len(self._sp_line[9])
    def get_alignment_length(self): 
        if not self.is_unmapped() and self.get_cigar_string() != '*':
            length=0
            all_cigars=re.findall('(\d+)([MXDINPSH])',self.get_cigar_string())
            for size,type in all_cigars:
                if type in ['M','X','D','N','P']:
                    length+=int(size)
            return length
        else:
            return 0
    def get_alignment_start(self): 
        """Equivalent to the position of the read except if the alignment starts with an insertion"""
        position = self.get_position()
        if not self.is_unmapped() and self.get_cigar_string() != '*':
            all_cigars=re.findall('(\d+)([MXDINPSH])',self.get_cigar_string())
            size,type = all_cigars[0]
            if type in ['I']:
                return position-int(size)
        return position
    
    
    
    def set_query_name(self,query_name): self._sp_line[0]=query_name
    def set_flag(self,flag): self._sp_line[1]=str(flag)
    def set_reference_name(self,reference_name): self._sp_line[2]=str(reference_name)
    def set_position(self,posistion): self._sp_line[3]=str(posistion)
    def set_mapping_quality(self, mapping_quality): self._sp_line[4]=str(mapping_quality)
    def set_cigar_string(self,cigar_string): self._sp_line[5]=cigar_string
    def set_mate_reference_name(self,mate_reference_name): self._sp_line[6]=str(mate_reference_name)
    def set_mate_pos(self, mate_pos):  self._sp_line[7]= str(mate_pos)
    def set_insert_size(self, insert_size): self._sp_line[8]=str(insert_size)
    def set_query_sequence(self, query_sequence): self._sp_line[9]=query_sequence
    def set_query_quality(self,query_quality): self._sp_line[10]=query_quality
    def is_paired(self):
        mask=1
        return self._check_bit(mask)
    
    def is_proper_pair(self):
        mask=2
        return self._check_bit(mask)
    
    def is_unmapped(self):
        mask=4
        return self._check_bit(mask)
    
    def is_mate_unmapped(self):
        mask=8
        return self._check_bit(mask)
    
    def is_reverse_strand(self):
        mask=16
        return self._check_bit(mask)
    
    def is_mate_reverse_strand(self):
        mask=32
        return self._check_bit(mask)
    
    def is_first_read(self):
        mask=64
        return self._check_bit(mask)
    
    def is_second_read(self):
        mask=128
        return self._check_bit(mask)
     
    def _check_bit(self, mask):
        if int(self.get_flag()) & mask == mask:
            return True
        else:
            return False
        
    def __str__(self):
        return '\t'.join(self._sp_line)+'\n'
    
    def set_tag(self,tag,type,new_value):
        added_new_value = False
        if len(self._sp_line)>11:
            for i in range(11,len(self._sp_line)):
                if self._sp_line[i].startswith(tag):
                    self._sp_line[i] = '%s:%s:%s' % (tag,type,new_value)
                    added_new_value = True
                    break
        if not added_new_value: 
            self._sp_line.append('%s:%s:%s' % (tag,type,new_value))
    
    def clear_all_tags(self):
        if len(self._sp_line)>11:
            self._sp_line = self._sp_line[:11]
    
    def clear_tag(self,requested_tag): #drop this tag from the record
        if len(self._sp_line)>10:
            for i in range(11,len(self._sp_line)):
                if (self._sp_line[i]).split(':')[0] == requested_tag:
                    self._sp_line.pop(i)
                    break
    #setting the strand, the mate's strand, and the read number are only required if we are reconstructing a read that went missing in the remapping
    def set_paired_flag(self, trueOrFalse = True):
        flag=1
        self._set_bit(flag,trueOrFalse)
    
    def set_proper_pair_flag(self, trueOrFalse = True):
        flag = 2
        self._set_bit(flag,trueOrFalse)
        
    def set_unmapped_flag(self, trueOrFalse = True):
        flag=4
        self._set_bit(flag,trueOrFalse)
        
    def set_mate_unmapped_flag(self, trueOrFalse = True):
        flag=8
        self._set_bit(flag,trueOrFalse)  
        
    def set_reverse_strand_flag(self, trueOrFalse = True):
        flag=16
        self._set_bit(flag,trueOrFalse)  

    def set_mate_reverse_strand_flag(self,trueOrFalse = True):
        flag=32
        self._set_bit(flag,trueOrFalse)  
    
    def set_first_read_flag(self, trueOrFalse = True):
        flag=64
        self._set_bit(flag,trueOrFalse)
    
    def set_second_read_flag(self, trueOrFalse = True):
        flag=128
        self._set_bit(flag,trueOrFalse)
        
    def set_duplicate_flag(self,trueOrFalse = True):
        flag=1024
        self._set_bit(flag,trueOrFalse) 
    
    def is_QC_failed(self):
        mask=512
        return self._check_bit(mask)
    
    def is_duplicate_read(self):
        mask=1024
        return self._check_bit(mask)
     
        
    def _set_bit(self,flag,trueOrFalse):
        #if trueOrFalse and (int(self.get_flag()) & flag <> flag): #request is to SET flag and the flag is currently unflagged
        #    self.set_flag(self.get_flag()+flag)
        #elif not trueOrFalse and (int(self.get_flag()) & flag == flag): #request to UNset the flag, which is currently true
        #    self.set_flag(self.get_flag()-flag)
        
        if trueOrFalse: #request is to set the flag
            self.set_flag(self.get_flag()|flag)
        elif not trueOrFalse and (int(self.get_flag()) & flag == flag): #request to UNset the flag, which is currently true
            self.set_flag(self.get_flag()-flag)
            
    def set_unaligned(self):
        self.set_unmapped_flag()
        self.set_proper_pair_flag(trueOrFalse=False)
        self.set_mapping_quality(0)
        self.set_cigar_string('*')
        self.set_insert_size('*')
        
        
            
    

def get_sam_from_mapview(mapview_line):
    """Create a sam line out of a mapview line. 
    The sam line may not be complete but contains all the information that can be gather from a maview line"""
    sp_line=mapview_line.strip().split()
    sp_line
    new_sam_sp_line=[]
    
    if sp_line[0].endswith('/1'):
        name=sp_line[0][:-2]
        mate_number=64
    elif sp_line[0].endswith('/2'):
        name=sp_line[0][:-2]
        mate_number=128
    else:
        name=sp_line[0]
        mate_number=0
    paired_flag=sp_line[5]
    flag=0;
    #construct the flag field
    if mate_number>0:
        flag+=1 #set the first byte to 1
        flag+=mate_number #set the seventh or eighth byte to 1
        if paired_flag==192:
            flag+=8 #set the forth byte to 1
        if paired_flag==64: # unmapped read
            flag+=4 #set the third byte to 1
        
    if sp_line[3] == "-":
        flag+=16 #set the fifth byte to 1

    new_sam_sp_line.append(name)
    new_sam_sp_line.append(str(flag))
    new_sam_sp_line.append(sp_line[1]) #chr
    new_sam_sp_line.append(sp_line[2]) #position
    new_sam_sp_line.append(sp_line[6]) #mapping quality
    new_sam_sp_line.append("%sM"%sp_line[13]) #Cigar
    new_sam_sp_line.append("*") #mate_chromsome (don't know)
    new_sam_sp_line.append(str( int(sp_line[2])+int(sp_line[4]))) #mate position
    new_sam_sp_line.append(sp_line[4]) #insert size
    new_sam_sp_line.append(sp_line[14]) #sequence
    new_sam_sp_line.append(sp_line[15]) #quality string
    new_sam_sp_line.append("MF:i:%s"%paired_flag)
    
    new_sam_sp_line.append("NM:i:%s"%sp_line[9])
    new_sam_sp_line.append("AM:i:%s"%sp_line[8])
    
    return '\t'.join(new_sam_sp_line)

def get_sam_record_from_mapview(mapview_line):
    """Create a Sam_record out of a mapview line. 
    The sam record may not be complete but contains all the information that can be gather from a maview line"""
    return Sam_record(get_sam_from_mapview(mapview_line))
    
class Sam_iterator(object):
    def __init__(self, input_stream):
        self.input_stream=input_stream
        self.iterator=None
        
    def read_file(self):
        for line in self.input_stream:
            yield Sam_record(line)
    
    def __iter__(self):
        if not self.iterator:
            self.iterator=self.read_file()
        return self.iterator

class Mapview_to_sam_iterator(object):
    def __init__(self, input_stream):
        self.input_stream=input_stream
        self.iterator=None
        
    def read_file(self):
        for line in self.input_stream:
            yield get_sam_record_from_mapview(line)
    
    def __iter__(self):
        if not self.iterator:
            self.iterator=self.read_file()
        return self.iterator



