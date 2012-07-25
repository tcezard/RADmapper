def test_overlap(s1,e1,s2,e2):
    """This function test if an overlap exists.
    if it does it return 0, if the first segment is before the second it returns 1 otherwise it returns -1 """
    if s1>e2: return -1
    if s2>e1: return 1
    else: return 0

def test_in_range(position, s,e):
    if position>e:return 1
    if position<s: return -1
    return 0
    
def overlap_boundaries(s1,e1,s2,e2):
    """calculate overlap of two segment of given coordinate.
    It assumes that s1<=e2 and s2<=e2 which means that the overlap exist.
    It also assumes that the coordinates are 1-based and end-inclusive. 
    So length of segment 1 = e1-s1+1 and length of segment 2 = e2-s2+1 """
    if s1>e2 or s2>e1: return 0
    else:
        return max(s1,s2),min(e1,e2)

def overlap_length(s1,e1,s2,e2):
    """calculate overlap of two segment of given coordinate.
    It assumes that s1<=e2 and s2<=e2.
    It also assumes that the coordinates are 1-based and end-inclusive. 
    So length of segment 1 = e1-s1+1 and length of segment 2 = e2-s2+1 """
    s,e=overlap_boundaries(s1,e1,s2,e2)
    return e-s+1

def get_overlap_if_exist(s1,e1,s2,e2):
    """This function returns the boundaries of an overlap."""
    if test_overlap(s1, e1, s2, e2)==0:
        return overlap_boundaries(s1, e1, s2, e2)
    else:
        return None

def merge_ranges(list_range, end_inclusive=True):
    """"Merge any overlapping range in the array and return a sorted version of the non overlapping range."""
    new_list_range=[]
    if len(list_range)>=2:
        list_range.sort()
        prev_range=list_range[0]
        for range in list_range:
            s1=prev_range[0]
            e1=prev_range[1]
            s2=range[0]
            e2=range[1]
            if (end_inclusive and (s1>e2+1 or s2>e1+1)) or \
            (not end_inclusive and (s1>e2 or s2>e1)): 
                new_list_range.append(prev_range)
                prev_range=range
            else:
                prev_range=min(s1,s2),max(e1,e2)
        new_list_range.append(prev_range)
    else:
        new_list_range=list_range
    return new_list_range

def get_length_of_list_range(array,end_inclusive=True):
    
    a=merge_ranges(array, end_inclusive)
    length=0
    for element in a:
        s=element[0]
        e=element[1]
        if end_inclusive:
            length+=e-s+1
        else:
            length+=e-s
    return length