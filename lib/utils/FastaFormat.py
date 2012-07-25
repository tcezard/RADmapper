class FastaReader : 
    """
    This class reads a file in the fasta format.
    No format checking is done.
    """
    def __init__(self, handle):
        """Initialize a new FastaReader.
        """
        self.handle = handle
        self._lookahead = None
        #Skip any text before the first record (e.g. blank lines)
        while True:
            line = self.handle.readline()
            if not line or line.startswith(">"):
                break
        self._lookahead = line

    def __iter__(self):
        return iter(self.next, None)

    def next(self):
        """Return the next record in the file"""
        line = self._lookahead
        if not line:
            return None
        assert line[0]==">", line
        header = line.rstrip().strip('>')
        line = self.handle.readline()
        lines=[]
        while line:
            if line[0] == ">": break
            else :
                lines.append(line.rstrip())
            line = self.handle.readline()
        self._lookahead = line
        return (header,''.join(lines))

def fasta_reader(handle):
    header=None
    sequence=None
    for line in handle:
        if line.startswith('>'):
            if sequence and header:
                yield (header, sequence)
            header = line.rstrip().strip('>')
            sequence=''
        else:
            sequence+=line.strip()
    if sequence and header:
        yield (header, sequence)

class FastaWriter :
    def __init__(self, handle, line_length, header=None):
        """Initialize a new FastaWriter.
        """
        self.handle = handle
        self.line_length=line_length
        self.buffer=[]
        self.buffer_length=0
        self.max_buffer_length=2000
        self.in_sequence=False
        if header:
            self.start_sequence(header)
        
    def start_sequence(self, header):
        self._flush_buffer(True)
        self.handle.write('>%s\n'%header)
        self.in_sequence=True
    
    def write_sequence(self,sequence):
        self.buffer.append(sequence)
        self.buffer_length+=len(sequence)
        if self.buffer_length>self.max_buffer_length:
            self._flush_buffer()
    
    def _flush_buffer(self,flush_all=False):
        string=''.join(self.buffer)
        all_lines=[]
        for pos in range(0,len(string),self.line_length):
            all_lines.append(string[pos:pos+self.line_length])
        if len(all_lines)>0:
            if len(all_lines[-1])==self.line_length or flush_all:
                self.handle.write('\n'.join(all_lines)+'\n')
                self.buffer=[]
                self.buffer_length=0
            else :
                self.handle.write('\n'.join(all_lines[:-1])+'\n')
                self.buffer=[all_lines[-1]]
                self.buffer_length=len(all_lines[-1]) 
        
    
    def close(self):
        self._flush_buffer(True)
    
    def __del__(self):
        self.close()
        
        
if __name__=="__main__":
    import sys
    fasta_reader = FastaReader(open(sys.argv[1]))
    for fasta_record in fasta_reader:
        header, sequence = fasta_record
        print header
        print sequence
if __name__=="1__main__":
    handle=open('/home/tcezard/temp/test.fasta','w')
    fasta_writer=FastaWriter(handle, 20, 'test')
    print dir(fasta_writer)
    fasta_writer.write_sequence('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA')
    fasta_writer.write_sequence('TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT')
    fasta_writer.start_sequence('test2')
    fasta_writer.write_sequence('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA')
    fasta_writer.write_sequence('TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT')
    #fasta_writer.close()
    