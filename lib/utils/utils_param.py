'''
Created on 8 Apr 2010

@author: tcezard
'''
import logging
import os
from utils.parameters import Config_file_error, Pipeline_places

#####################################
#        Exception utilities        #
#####################################
class Resource_File_Exception(Exception):
    """Raised when an error occur because a resource file doen't exists"""
    pass


#############################################
#        parameters and config files        #
#############################################

pipeline_parameters_obj={}
def get_pipeline_parameters(env="PIPELINE_CONFIG"):
    """This method retrieve the configuration file stored in the environment variable."""
    if not pipeline_parameters_obj.get(env):
        from utils.config import ConfigReader
        from utils.parameters import Pipeline_param
        ##################
        # get the config environment
        try :
            pipeline_config=ConfigReader(env)
        except IOError, e:
            logging.critical("%s environment variable is not set properly"%(env))
            raise Config_file_error("%s environment variable is not set properly"%(env))
        pipeline_param=Pipeline_param(pipeline_config)
        pipeline_parameters_obj[env]=pipeline_param
    return pipeline_parameters_obj.get(env)

def get_pipeline_places(env="PIPELINE_PLACES"):
    """This method retrieve the configuration file stored in the environment variable."""
    if not pipeline_parameters_obj.get(env):
        from utils.config import ConfigReader
        from utils.parameters import Pipeline_param
        ##################
        # get the config environment
        try :
            pipeline_config=ConfigReader(env)
        except IOError, e:
            try :
                pipeline_config=ConfigReader(file=get_pipeline_parameters().get_pipeline_place())
            except Config_file_error, e:
                message="%s environment variable is not set properly\n%s"%(env, str(e))
                logging.critical(message)
                raise Config_file_error(message)
        pipeline_places = Pipeline_places(pipeline_config)
        pipeline_parameters_obj[env]=pipeline_places
    return pipeline_parameters_obj.get(env)

wtss_param_obj={}
def get_wtss_parameter(env="SOLEXA_CONFIG"):
    """This method retrieve the configuration file stored in the environment variable."""
    if not wtss_param_obj.get(env):
        from utils.config import ConfigReader
        from utils.parameters import WTSS_param
        ##################
        # get the config environment
        try :
            wtss_config=ConfigReader(env)
            from utils.config import display_config
        except IOError, e:
            logging.critical("%s environment variable is not set properly"%(env))
            raise Config_file_error("%s environment variable is not set properly"%(env))
        wtss_param=WTSS_param(wtss_config)
        wtss_param_obj[env]=wtss_param
    return wtss_param_obj.get(env)


def check_input_file(input_file, pipe_allowed=False):
    arg_pass=True
    if not input_file:
        logging.error("You must specify an input file.")
        arg_pass=False
    elif input_file=="PIPE":
        if not pipe_allowed:
            logging.error("PIPE option is not allowed for the input file")
            arg_pass=False
    elif not os.path.exists(input_file):
        logging.error("You must specify a existing input file %s does not exist"%input_file)
        arg_pass=False
    return arg_pass

def check_output_file(output_file, pipe_allowed=False):
    arg_pass=True
    if not output_file:
        logging.error("You must specify an output file.")
        arg_pass=False
    elif output_file=="PIPE":
        if not pipe_allowed:
            logging.error("PIPE option is not allowed for the output file")
            arg_pass=False
    else:
        parent_dir=os.path.dirname(output_file)
        if parent_dir=='':
            parent_dir='.'
        if not os.path.exists(parent_dir):
            logging.error("directory %s does not exist cannot create %s"%(parent_dir,output_file))
            arg_pass=False
    return arg_pass


