'''
Created on 26 Mar 2010

@author: tcezard
'''
from overlap import test_overlap, merge_ranges, overlap_boundaries, test_in_range
from utils import utils_logging, benchmarck_timer


class Binnary_search(object):
    '''This class search for an overlap between a list of possibly overlapping range and a query.
    The requirements are:
     - The original array is an array of tuple where the first two elements are the start and the end of the range.
     - A start must always be lower or equal to its end.'''
    def __init__(self, array, sorted_query=False, close_query=False):
        ''''''
        self.array=array
        self.ranges=self.split_range(self.array)
        self.sorted_query=sorted_query
        self.low=0
        if close_query:
            self.sorted_query=close_query
        self.close_query=close_query
    
        
    def ovarlap_search(self,position):
        list_of_array_value=[]
        indice=self._search(position)
        if indice is not None:
            (s1,e1,set_indices)=self.ranges[indice]
            if set_indices is not None:
                for indice in set_indices:
                    list_of_array_value.append(self.array[indice])
        return list_of_array_value
    
    def is_ovarlap(self,position):  
        indice=self._search(position)
        if indice is not None:
            (s1,e1,set_indices)=self.ranges[indice]
            if set_indices is not None:
                return True
        return False
    
    def _search(self,position):
        '''Does the actual binary search in the range using an overlap function'''
        high = len(self.ranges) - 1
        original_low=self.low
        if not self.sorted_query:
            self.low = 0
        if self.close_query:
            mid = self.low
        else:
            mid = (self.low + high) / 2
        lastMid=-1
        while self.low <= high:
            res=self._overlap_range(mid, position)
            if res==0:
                self.low = mid
                return mid
            else :
                if res > 0 :
                    self.low = mid
                else :
                    high = mid
            if lastMid==mid:
                break;
            lastMid=mid
            mid = (self.low + high) / 2
        self.low=original_low
        return None;
    
    def split_range(self,array):  
        """Split an array of range (start-end) in a set of non-overlapping range linking to all the indices in the original array.
        It returns an array of tuple each containing 3 elements (start, end, set_of_annotation_id).""" 
        points = [] # list of (offset, plus/minus, id) tuples
        for pos in range(len(array)):
            start=array[pos][0]
            end=array[pos][1]
            points.append((start,'+',pos))
            points.append((end,'-',pos))
        points.sort()
        
        ranges = [] # output list of (start, end, symbol_set) tuples
        current_set = set()
        last_start = None
        for offset, pm, pos in points:
            if pm == '+':
                #TODO avoid outputting empty or trivial ranges
                if last_start<offset:
                    if len(current_set)==0:
                        ranges.append((last_start,offset-1,None))
                    else:
                        ranges.append((last_start,offset-1,current_set.copy()))
                current_set.add(pos)
                last_start = offset
            elif pm == '-':
                # Getting a minus without a last_start is impossible here, so not handled
                if last_start<offset:
                    ranges.append((last_start,offset,current_set.copy()))
                current_set.remove(pos)
                #if len(current_set)>0:
                last_start = offset+1
                #Now keep also the ranges with nothing in them
                #else:
                #    last_start=None
        # Finish off
        if last_start is not None:
            ranges.append((last_start,offset-1,current_set))
        return ranges
    
    def _overlap_range(self,indice, position):
        s1,e1,i= self.ranges[indice]
        return test_in_range(position,s1,e1)
                


class Binnary_search_overlap(object):
    '''This class search for an overlap between a list of possibly overlapping range and a query.
    The requirements are:
     - The original array is an array of tuple where the first two elements are the start and the end of the range.
     - A start must always be lower or equal to its end.'''
    def __init__(self, array, sorted_query=False, close_query=False):
        ''''''
        self.array=array
        self.ranges=self.split_range(self.array)
        self.sorted_query=sorted_query
        self.low=0
        if close_query:
            self.sorted_query=close_query
        self.close_query=close_query
    
    
        
    def ovarlap_search(self,s2,e2):
        array_of_indices=self._search(s2,e2)
        all_original_indices_and_length={}
        if array_of_indices is not None:
            for indice in array_of_indices:
                (s1,e1,set_indices)=self.ranges[indice]
                s, e = overlap_boundaries(s1, e1, s2, e2)
                if e<s:
                    print 'ERROR end %s is greater than start %s'%(e,s)
                    print s1, e1, s2, e2
                for original_i in set_indices:
                    if all_original_indices_and_length.has_key(original_i):
                        list_range=all_original_indices_and_length.get(original_i)
                        list_range.append((s,e))
                        list_range=merge_ranges(list_range)
                        all_original_indices_and_length[original_i]=list_range
                    else:
                        all_original_indices_and_length[original_i]=[(s,e)]
        list_of_array_value_and_length=[]
        for indice in all_original_indices_and_length.keys():
            list_of_array_value_and_length.append((self.array[indice],
                                                  all_original_indices_and_length.get(indice)))
        return list_of_array_value_and_length
            
    
    def _search(self,start,end):
        '''Does the actual binary search in the range using an overlap function'''
        array_of_indices=[]
        high = len(self.ranges) - 1
        if not self.sorted_query:
            self.low = 0
        if self.close_query:
            mid = self.low
        else:
            mid = (self.low + high) / 2
        lastMid=-1
        while self.low <= high:
            res=self._overlap_range(mid, start,end)
            if res==0:
                array_of_indices.append(mid)
                n = 1
                while mid-n > self.low and self._overlap_range(mid-n, start,end)==0:
                    array_of_indices.append(mid - n)
                    n+=1
                n = 1
                while  mid+n < high and self._overlap_range(mid+n, start,end)==0:
                    array_of_indices.append(mid + n)
                    n+=1
                return array_of_indices
            else :
                if res > 0 :
                    self.low = mid
                else :
                    high = mid
            if lastMid==mid:
                break;
            lastMid=mid
            mid = (self.low + high) / 2
        return None;
    
    def split_range(self,array):  
        """Split an array of range (start-end) in a set of non-overlapping range linking to all the indices in the original array.
        It returns an array of tuple each containing 3 elements (start, end, set_of_annotation_id).""" 
        points = [] # list of (offset, plus/minus, id) tuples
        for pos in range(len(array)):
            start=array[pos][0]
            end=array[pos][1]
            points.append((start,'+',pos))
            points.append((end,'-',pos))
        points.sort()
        
        ranges = [] # output list of (start, end, symbol_set) tuples
        current_set = set()
        last_start = None
        for offset, pm, pos in points:
            if pm == '+':
                if last_start is not None:
                    #TODO avoid outputting empty or trivial ranges
                    if last_start<offset:
                        ranges.append((last_start,offset-1,current_set.copy()))
                        
                current_set.add(pos)
                last_start = offset
            elif pm == '-':
                # Getting a minus without a last_start is impossible here, so not handled
                if last_start<offset:
                    ranges.append((last_start,offset,current_set.copy()))
                current_set.remove(pos)
                if len(current_set)>0:
                    last_start = offset+1
                else:
                    last_start=None
        # Finish off
        if last_start is not None:
            ranges.append((last_start,offset-1,current_set))
        return ranges
    
    def _overlap_range(self,indice, start,end):
        s1,e1,i= self.ranges[indice]
        return test_overlap(s1,e1,start,end)


def read_exon_capture_file(exon_capture_file, extension=0):
    open_exon_capture=utils_logging.open_input_file(exon_capture_file, pipe=False)
    all_segments_per_chr={}
    for line in open_exon_capture:
        if line.startswith('#') or line.startswith('track'):
            continue
        sp_line=line.strip().split()
        chr=sp_line[0]
        start=int(sp_line[1])-extension
        end=int(sp_line[2])+extension
        all_segments=all_segments_per_chr.get(chr)
        if all_segments is None:
            all_segments=[]
            all_segments_per_chr[chr]=all_segments
        all_segments.append((start,end))
    open_exon_capture.close()
    return all_segments_per_chr


if __name__=='1__main__':
    utils_logging.init_logging()
    exon_capture_file='/home/tcezard/projects/exons_capture_evaluation/agilent/source/SureSelect_All_Exon_G3362_with_names_hg19.bed'
    all_target_per_chr=read_exon_capture_file(exon_capture_file, extension=0)
    chr='chr1'
    #chr1    888551  888671
    positions=range(88550,888673)
    #positions=range(888551,888671)
    all_target=all_target_per_chr.get(chr)
    all_target_search=Binnary_search(all_target)
    #all_target_search=Binnary_search(all_target,close_query=True)
    #all_target_search=Binnary_search(all_target,sorted_query=True)
    for position in positions:
        benchmarck_timer.start('search')
        print position , all_target_search.is_ovarlap(position)
        all_target_search.is_ovarlap(position)
        benchmarck_timer.stop('search')
    benchmarck_timer.print_all()
    
if __name__=='__main__':
    from pprint import pprint
    all_target=[(10,20),
                (50,85),
                (104,221),
                (1100,1101),
                (1110,1120),
                (1109,1112)]
    positions_out=[]
    positions_out.extend(range(0,10))
    positions_out.extend(range(21,50))
    positions_out.extend(range(86,104))
    
    positions_in=[12,125,221,1109,1113]
    all_target_search=Binnary_search(all_target)
    for position_out in positions_out:
        try:
            assert not all_target_search.is_ovarlap(position_out)
        except AssertionError:
            print 'ERROR: position %s should not overlap the targets'%(position_out)
            pprint(all_target)
    for position_in in positions_in:
        try:
            assert all_target_search.is_ovarlap(position_in)
        except AssertionError:
            print 'ERROR position %s should overlap the targets'%(position_in)
            pprint(all_target)
    
        
    