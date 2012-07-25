"""
Reads the configuration file and makes it each
section available as an attribute. Each key/value pair within each
section then becomes a Dictionary.

By default each value is a string. Additional coercion may be done
to parse those strings into more useful data structures.
"""

import os
from exceptions import IOError
from ConfigParser import ConfigParser

class ConfigReader:
    """
    Provides a common, convenient way of reading configuration
    parameters for working within the SBS Solexa Bioinformatics
    environment.
    """
    
    def __init__( self, env=None, file=None):
        """
        Instantiating a SolexaConfig object causes the config
        file to be read and parsed.
        """
        cfg = self._read_config(env,file)
        if not cfg:
            raise IOError()
        for section in cfg.sections():
            section_dict = {}
            for k,v in cfg.items(section):
                section_dict[k] = v
            setattr(self, section, section_dict)
        self._coerce_values()
        
    def _coerce_values(self):
        """
        Do some additional munging of the config data
        Mash Until No Good or Modified Until Not Guessed Easily ???
        """
        # if we are going to 'munge' then we are going to write
        # our Python code as if it were Perl ...
        if hasattr(self, 'email') and self.email.get('notification_emails'):
            self.email['notification_emails'] = \
            [x.strip() for x in self.email['notification_emails'].split(',')]
    
    def _read_config(self, env, file):
        """
        Searches for the appropriate sbs-solexa-config.txt file
        and parses it.
        
        Currently the config is determined by the SOLEXA_CONFIG
        environment variable.
        """
        # XXX Put checks in place that help to ensure that the 
        # no one is accidently running dev code against the
        # production environment.
        if not file:
            config_path = os.getenv(env)
        else:
            config_path=file
        self.config_path=config_path
        if config_path and os.path.isfile(config_path):
            cfg = ConfigParser()
            cfg.read(config_path)
            return cfg
        else: return None
        
        

def display_config(config):
    """
    Return the contents of the config object as plain-text formatted string
    """
    # untested and fragile, only use for educational purposes!
    out = ''
    for k in dir(config):
        obj = getattr(config,k)
        if type(obj) == type({}):
            out += '[%s]\n' % (k)
            for k,v in obj.items():
                out += "%s = %s\n" % (k,v)
            out += '\n'
            
    return out
