import logging
from utils import utils_commands
from io import StringIO


TYPE_PRINT='print_out'
TYPE_RUN_LOCAL='run'
TYPE_RUN_SGE='sge'

RUNNER_TYPE=[TYPE_PRINT, TYPE_RUN_LOCAL]

class Command_runner():
    def __init__(self, run_type=TYPE_PRINT):
        if run_type in RUNNER_TYPE:
            self.type_r=run_type
        else:
            raise Exception()
    
    def run_command(self,command, **kwargs):
        """Run the command and return the code that this command returned"""
        if command.__class__==''.__class__:
            pass
        elif command.__class__==[].__class__:
            pass
        return_code=0
        if self.type_r==TYPE_PRINT:
            print command
        elif self.type_r==TYPE_RUN_LOCAL:
            logging.info(command)
            return_code=utils_commands.launchCommandLocally(command, nice=False, **kwargs)
        else:
            logging.error('The type %s is not handle by the Command_runner class'%self.type_r)
            return_code=1
        return return_code
            
    def run_command_no_wait(self,command):
        """Run the command and return the Process object that control that running command."""
        process=None
        if self.type_r==TYPE_PRINT:
            print command
        elif self.type_r==TYPE_RUN_LOCAL:
            logging.info(command)
            process = utils_commands._launchCommand_no_wait(command)
        else:
            logging.error('The type %s is not handle by the Command_runner class'%self.type_r)
            process = None
        return process
    
    def set_command_to_run_localy(self):
        self.type_r=TYPE_RUN_LOCAL
        
    def set_command_to_print(self):
        self.type_r=TYPE_PRINT
        
_default_instance=Command_runner()

def run_command(command):
    return _default_instance.run_command(command)

def run_command_no_wait(command):
    return _default_instance.run_command_no_wait(command)

def set_command_to_run_localy():
    _default_instance.set_command_to_run_localy()

def set_command_to_print():
    _default_instance.set_command_to_print()
    