from utils.parameters import WTSS_param
import logging
import utils
import os
from utils import utils_logging, utils_param

INPUT_MAP="map"
INPUT_MAPVIEW="mapview"
INPUT_BAM="bam"
INPUT_SAM="sam"

ALIGNEMENT_INPUT_TYPE=[INPUT_MAP,INPUT_MAPVIEW,INPUT_BAM,INPUT_SAM]


def open_alignement_file(input_file, input_file_type=None, pipe=True):
    """This method opens an alignment file in sam, bam, map, mapview format and returns an open file that will generate either mapview of sam.
    It returns a tuple (open_input, input_file_type) where open_input is an open file handle and input_file_type the type of line."""
    open_input=None
    if input_file_type is None:
        dummy,input_file_type=os.path.splitext(input_file)
        if input_file_type.startswith('.'):
            input_file_type=input_file_type[1:]
    if input_file_type==INPUT_MAPVIEW or input_file_type==INPUT_SAM:
        open_input=utils_logging.open_input_file(input_file,pipe)
    elif input_file_type==INPUT_MAP:
        if os.path.exists(input_file):
            wtss_param=utils_param.get_wtss_parameter()
            maq_bin=utils.findNewestMaqVersionForMapfile(wtss_param.get_maq_dir(), input_file)
            open_input=utils.get_mapview_stream(maq_bin, input_file)
            input_file_type=INPUT_MAPVIEW
        else:
            logging.error("file not found: %s"%input_file)
    elif input_file_type==INPUT_BAM:
        if os.path.exists(input_file):
            wtss_param=utils_param.get_wtss_parameter()
            samtools_bin=os.path.join(wtss_param.get_samtools_dir(),"samtools")
            open_input=utils.get_sam_stream(input_file,samtools_bin)
            input_file_type=INPUT_SAM
        else:
            logging.error("file not found: %s"%input_file)
    else:
        logging.error("Non recognize input type %s for %s"%(input_file_type,input_file))
    return open_input,input_file_type
        
    
    
    
    