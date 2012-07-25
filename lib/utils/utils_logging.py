'''
Created on 8 Apr 2010

@author: tcezard
'''
import sys
import gzip
import bz2
import logging



##########################################
#        logging and file opening        #
##########################################

all_handlers={}

def init_logging(output_level=logging.INFO, file_level=logging.INFO, log_file_name=None, overwrite=False, logger_name='', formatter=None):
    """initialise the logging with default parameter.
If use several time with same input file only the level of the logging will be change.
@param output_level: The logging level of output for the stdout. If set to None the stdout is not set or not changed. default: INFO.
@param file_level: The logging level of output for the log file. if set to None the logging in a file is not set or not changed. default: INFO.
@param log_file_name: The name of the file used for the logging. if set to None the logging in a file is not set or not changed. default: None.
@param overwrite: if set to true an existing log file will overwritten, otherwise the next log will be appended. default: None.
@param logger_name: Initialise non root logger that will be accessible through logging.getLogger(name). default: ''.
@param logger_name: set non default formatter that will be use to format the logging. default: ''.
Example:
init_logging() --> set the stdout only to INFO
init_logging(log_file_name="any_file") --> set the stdout and the log file to INFO
init_logging(output_level=None, log_file_name="any_file") --> set the log file only to INFO don't touch the stdout
"""
    if logger_name is not '':
        logging.root.setLevel(logging.NOTSET)
    
    if formatter is None:
        formatter=logging.Formatter('%(levelname)s %(message)s')
    
    if output_level is not None:
        #user wants to set the something to the output logger
        handler_info=all_handlers.get(logger_name+"std_out")
        
        if handler_info is None:
            #create a new handler
            console_handler=logging.StreamHandler(sys.stdout)
            current_logger=logging.getLogger(logger_name)
            current_logger.addHandler(console_handler)
            current_logger.setLevel(logging.NOTSET)
            current_logger.propagate=0
            all_handlers[logger_name+"std_out"]=(console_handler,output_level)
        else:
            (console_handler,current_level)=handler_info
        console_handler.setLevel(output_level)
        console_handler.setFormatter(formatter)

    if log_file_name and file_level is not None:
        handler_info=all_handlers.get(logger_name+log_file_name)
        if handler_info is None:
            #create a new handler
            if overwrite:
                file_handler = logging.FileHandler(log_file_name, 'w')
            else:
                file_handler = logging.FileHandler(log_file_name, 'a')
            current_logger=logging.getLogger(logger_name)
            current_logger.addHandler(file_handler)
            current_logger.setLevel(logging.NOTSET)
            current_logger.propagate=0
            all_handlers[logger_name+log_file_name]=(file_handler,file_level)
        else:
            (file_handler,current_level)=handler_info
        file_handler.setLevel(file_level)
        file_handler.setFormatter(formatter)
        
    return logging.getLogger(logger_name)

def remove_logging_std_out_handler(logger_name=''):
    """remove the handler that log on the standard out if it was set using init_logging"""
    console_output_level = all_handlers.pop(logger_name+"std_out",None)
    if console_output_level:
        (console, dummy)=console_output_level
        logging.getLogger(logger_name).removeHandler(console)
    return console_output_level

def remove_logging_file_handler(log_file_name, logger_name=''):
    """remove the handler that log on that specific log file if it was set using init_logging"""
    if log_file_name:
        logFile_file_level=all_handlers.pop(logger_name+log_file_name,None)
        if logFile_file_level:
            (logFile,dummy)=logFile_file_level
            logging.getLogger('').removeHandler(logFile)
    return logFile_file_level

def add_logging_std_err_handler(output_level=logging.INFO, logger_name=''):
    """add the handler that log on the standard err."""
    console = logging.StreamHandler(sys.stderr)
    console.setLevel(output_level)
    console.setFormatter(logging.Formatter('%(levelname)s %(message)s'))
    logging.getLogger(logger_name).addHandler(console)
    all_handlers[logger_name+"std_err"]=(console,output_level)
    
def change_log_stdout_to_log_stderr(logger_name=''):
    """Move the logging happening in stdout to stderr."""
    console_output_level=remove_logging_std_out_handler(logger_name)
    if console_output_level:
        (dummy,output_level)=console_output_level
        add_logging_std_err_handler(output_level,logger_name)
            
def open_input_file(input_file, pipe=True):
    """open the input depending on the value:
PIPE: take from standard input
file name ending with .gz use gzip module to open
otherwise open normally"""
    if pipe and input_file=="PIPE":
        return sys.stdin
    elif input_file.endswith('.gz'):
        return gzip.open(input_file)
    elif input_file.endswith('.bz2'):
        return bz2.BZ2File(input_file, 'r')
    else:
        return open(input_file)

def open_output_file(output_file, pipe=True):
    """open the input depending on the value:
PIPE: take from standard input
file name ending with .gz use gzip module to open
otherwise open normally"""
    if pipe and output_file=="PIPE":
        change_log_stdout_to_log_stderr()
        return sys.stdout
    elif output_file.endswith('.gz'):
        return gzip.open(output_file, 'w')
    elif output_file.endswith('.bz2'):
        return bz2.BZ2File(output_file, 'w')
    else:
        return open(output_file, 'w')
    
    
if __name__=='__main__':
    logger0=init_logging(logger_name='', formatter=logging.Formatter('%(levelname)s %(message)s -- root '))
    logger1=init_logging(logger_name='first', formatter=logging.Formatter('%(levelname)s %(message)s -- first '))
    logger2=init_logging(logger_name='second', formatter=logging.Formatter('%(levelname)s %(message)s -- second '))
    logger2=init_logging(logger_name='third', formatter=logging.Formatter('%(levelname)s %(message)s -- third '))
    
    logging.info('test1')
    logging.info('test2')
    logging.getLogger('first').info('test3')
    logging.getLogger('second').info('test4')
   