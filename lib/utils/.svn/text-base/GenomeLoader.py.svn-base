import sys,logging
from utils.FastaFormat import FastaReader
from utils import utils_logging
from utils.utils_logging import open_input_file


class GenomeLoader:
    """Genome Loader take a fasta file and try to find a given chromosome in it.
    You can specify so it also keep the all the loaded unused chromosome into memory until they are required (keep_in_memory=True).
    You can also specify so it keep the all the loaded chromosome until the object is destroyed (keep_until_done=True).
    You can provide a prefix that will be added to the chromsome names if not already there.
    """
    def __init__(self,genome_file, keep_in_memory=True, keep_until_done=False, prefix=''):
        self.open_genome_file=open_input_file(genome_file)
        self.reader=FastaReader(self.open_genome_file)
        self.keep_in_memory=keep_in_memory
        self.keep_until_done=keep_until_done
        self.prefix=prefix
        self.all_chr={}
    
    def load_all(self):
        self.get_chr('***********************************************')
        return self.all_chr.keys()
    
    def get_chr(self, chr):
        if not chr.startswith(self.prefix):
            chr=self.prefix+chr
            
        if self.keep_until_done: #return if loaded already
            fasta_seq=self.all_chr.get(chr)
        else:               #remove if loaded already
            fasta_seq=self.all_chr.pop(chr,None)
        if fasta_seq:
            (header,seq)=fasta_seq
            logging.debug('return %s'%header)
            return fasta_seq
        curr_chr=''
        seq=''
        while not curr_chr==chr:
            fasta_seq=self.reader.next()
            if fasta_seq:
                (header,seq)=fasta_seq
                logging.debug('load %s'%header)
                curr_chr=header.split()[0]
                if not curr_chr.startswith(self.prefix):
                    curr_chr=self.prefix+curr_chr
                if (self.keep_in_memory and not curr_chr==chr) or self.keep_until_done:
                    logging.debug('keep %s'%header)
                    self.all_chr[curr_chr]=fasta_seq
            else: break
        return fasta_seq
    
    def next(self):
        fasta_seq=None
        if len(self.all_chr)>0:
            chr=self.all_chr.keys[0]
        if len(self.all_chr)>0:
            chr=self.all_chr.keys()[0]
            if self.keep_until_done: #return if loaded already
                fasta_seq=self.all_chr.get(chr)
            else:               #remove if loaded already
                fasta_seq=self.all_chr.pop(chr,None)
        if not fasta_seq:
            fasta_seq=self.reader.next()
            if fasta_seq:
                (header,seq)=fasta_seq
                curr_chr=header.split()[0]
                logging.debug('load %s'%header)
                if self.keep_until_done:
                    logging.debug('keep %s'%header)
                    self.all_chr[curr_chr]=fasta_seq
        if fasta_seq:
            logging.debug('return %s'%header)
            return fasta_seq
        return None
    
    def __iter__(self):
        return iter(self.next, None)
    
    def __del__(self):
        self.open_genome_file.close()
        self.all_chr=None
        
    def close(self):
        self.__del__()
        

if __name__=="1__main__":
    utils_logging.init_logging(logging.DEBUG)
    import benchmarck_timer
    genome_file='/home/tcezard/genomes/homo_sapiens/hg19/hg19.fa'
    genome_loader=GenomeLoader(genome_file)
    benchmarck_timer.start('load')
    sequence=genome_loader.get_chr('dfhgiduh')
    benchmarck_timer.stop('load')
    benchmarck_timer.start('retrieve')
    sequence=genome_loader.get_chr('1')
    benchmarck_timer.stop('retrieve')
    benchmarck_timer.print_all()
    
if __name__=="__main__":
    utils_logging.init_logging(logging.DEBUG)
    import benchmarck_timer
    genome_file='/home/tcezard/genomes/homo_sapiens/hg19/hg19.fa'
    genome_loader=GenomeLoader(genome_file)
    for fasta_rec in genome_loader:
        header, sequence=fasta_rec
        print header
    
    
    