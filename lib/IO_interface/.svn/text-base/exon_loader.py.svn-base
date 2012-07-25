'''
Created on 5 Mar 2010

@author: tcezard
'''
import logging
from utils import DNA_tools, utils_logging
from annotation_loader import Exon_annotation_Retriver


if __name__=="__main__":
    utils_logging.init_logging()
    gff_file='/home/tcezard/projects/2010053_Tom_Little_RNAseq_2/for tim/dmag_ep24augmap2an2.gff'
    exon_reader = Exon_annotation_Retriver(annotation_file=gff_file)
    annotations = exon_reader.get_annotation_from_chr('scaffold03376')
    for anno in annotations:
        start,end,transcript,gene,transcript_start,transcript_end,cds_start,cds_end,chr, dummy= anno
        #if transcript=='AUGep24bs03376g85t1':
        print anno
