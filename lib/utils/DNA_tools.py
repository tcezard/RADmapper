
STOP_CODON='*'

base_complement ={'N':'N','A':'T','C':'G','G':'C','T':'A','n':'n','a':'t','c':'g','g':'c','t':'a','*':'*'}
def rev_complements(sequence):
    """Method written by Johnson Pang to reverse complement the given dna sequence."""
    reverse_seq=sequence[::-1]
    return complements(reverse_seq)

def complements(sequence):
    """Method written by Johnson Pang to complement the given dna sequence."""
    #rev_comp_seq=[complement.get(base) for base in complement_seq ]
    comp_seq=[]
    for base in sequence:
        comp_base=base_complement.get(base)
        if not comp_base:
            print "ERROR %s is not valid a valid base"%base
        else :
            comp_seq.append(comp_base)
    
    return ''.join(comp_seq)
gencode = {
    'ATA':'I',    # Isoleucine
    'ATC':'I',    # Isoleucine
    'ATT':'I',    # Isoleucine
    'ATG':'M',    # Methionine
    'ACA':'T',    # Threonine
    'ACC':'T',    # Threonine
    'ACG':'T',    # Threonine
    'ACT':'T',    # Threonine
    'AAC':'N',    # Asparagine
    'AAT':'N',    # Asparagine
    'AAA':'K',    # Lysine
    'AAG':'K',    # Lysine
    'AGC':'S',    # Serine
    'AGT':'S',    # Serine
    'AGA':'R',    # Arginine
    'AGG':'R',    # Arginine
    'CTA':'L',    # Leucine
    'CTC':'L',    # Leucine
    'CTG':'L',    # Leucine
    'CTT':'L',    # Leucine
    'CCA':'P',    # Proline
    'CCC':'P',    # Proline
    'CCG':'P',    # Proline
    'CCT':'P',    # Proline
    'CAC':'H',    # Histidine
    'CAT':'H',    # Histidine
    'CAA':'Q',    # Glutamine
    'CAG':'Q',    # Glutamine
    'CGA':'R',    # Arginine
    'CGC':'R',    # Arginine
    'CGG':'R',    # Arginine
    'CGT':'R',    # Arginine
    'GTA':'V',    # Valine
    'GTC':'V',    # Valine
    'GTG':'V',    # Valine
    'GTT':'V',    # Valine
    'GCA':'A',    # Alanine
    'GCC':'A',    # Alanine
    'GCG':'A',    # Alanine
    'GCT':'A',    # Alanine
    'GAC':'D',    # Aspartic Acid
    'GAT':'D',    # Aspartic Acid
    'GAA':'E',    # Glutamic Acid
    'GAG':'E',    # Glutamic Acid
    'GGA':'G',    # Glycine
    'GGC':'G',    # Glycine
    'GGG':'G',    # Glycine
    'GGT':'G',    # Glycine
    'TCA':'S',    # Serine
    'TCC':'S',    # Serine
    'TCG':'S',    # Serine
    'TCT':'S',    # Serine
    'TTC':'F',    # Phenylalanine
    'TTT':'F',    # Phenylalanine
    'TTA':'L',    # Leucine
    'TTG':'L',    # Leucine
    'TAC':'Y',    # Tyrosine
    'TAT':'Y',    # Tyrosine
    'TAA':STOP_CODON,    # Stop
    'TAG':STOP_CODON,    # Stop
    'TGC':'C',    # Cysteine
    'TGT':'C',    # Cysteine
    'TGA':STOP_CODON,    # Stop
    'TGG':'W',    # Tryptophan
}

def strand_is_positive(strand):
    if str(strand)=="1" or str(strand)=="+":
        return True
    else:
        return False
    
def strand_is_same(strand1, strand2):
    return strand_is_positive(strand1) == strand_is_positive(strand2)


if __name__=="1__main__":
    seq='ATCGGGTC'
    print seq
    print rev_complements(seq)
    print rev_complements(seq)[::-1]
    
if __name__=="__main__":
    print "%s should be True"%(strand_is_same('+',1))
    print "%s should be False"%(strand_is_same('-',1))
    print "%s should be False"%(strand_is_same('-',"1"))
    print "%s should be True"%(strand_is_same('-',-1))
    print "%s should be False"%(strand_is_same(-1,"1"))
    

    
    