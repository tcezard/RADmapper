'''
Created on 8 Apr 2010

@author: tcezard
'''
from threading import Thread
import logging
from subprocess import Popen, PIPE
import os

        
###################################
#        Command utilities        #
###################################
class Output_Logger(Thread):
    """Close the stdin and log the stdout and stderr using logging module.
    All the logging is done in a separate thread so it won't block because a buffer is full.
    You can specify the logger name that should be used but this logger should already be created in the logging module."""
    def __init__(self,stdin,stdout,stderr, logger_name=''):
        self.stdin=stdin
        self.stdout=stdout
        self.stderr=stderr
        if logger_name is not None:
            self.logger=logging.getLogger(logger_name)
        else:
            self.logger=None
        Thread.__init__(self)
        
    def run(self):
        import select
        read_set = []
        #the stdin need to be closed no information come through it or the select will freeze.
        if self.stdin:
            self.stdin.close()
        
        if self.stdout:
            read_set.append(self.stdout)
        if self.stderr:
            read_set.append(self.stderr)

        while read_set :
            rlist, wlist, xlist = select.select(read_set, [], [])

            if self.stdout in rlist:
                line = self.stdout.readline().strip()
                if not line:
                    self.stdout.close()
                    read_set.remove(self.stdout)
                elif self.logger:
                    self.logger.info(line)

            if self.stderr in rlist:
                line = self.stderr.readline().strip()
                if not line:
                    self.stderr.close()
                    read_set.remove(self.stderr)
                elif self.logger:
                    self.logger.error(line)


def _launchCommand(command, verbose=True, logger_name='', cwd=None):
    if not verbose:
        logging.getLogger(logger_name).setLevel(logging.CRITICAL)
    process=_launchCommand_no_wait(command, logger_name, cwd)
    return_code=process.wait()
    if not verbose:
        logging.getLogger('').setLevel(logging.NOTSET)
    logging.info("pid=%s returned %s"%(process.pid,return_code))
    return return_code

def _launchCommand_no_wait(command, logger_name='', cwd=None):
    process=Popen(command, stdin=PIPE, stdout=PIPE, stderr=PIPE, shell=True, executable='/bin/bash', cwd=cwd)
    logger=Output_Logger(process.stdin, process.stdout, process.stderr, logger_name)
    logger.setDaemon(True)
    logger.start()
    return process

def launchCommandOnServer(server,command, verbose=True):
    """
    Launch the command on the specified server. 
    """
    username=os.environ['USER']
    command='ssh %s@%s "nice %s ;"'%(username,server,command)
    return _launchCommand(command, verbose)


def mqsubCommandOnCluster(mqsub_file, cluster_name, indir=None):
    """
    Launch the command on the cluster and return the jobs id in an array. 
    """
    username=os.environ['USER']
    if indir and os.path.isdir(indir):
        command='cd %s; /opt/mqtools/bin/mqsub --file %s;'%(os.path.abspath(indir), mqsub_file)
    else:
        command='/opt/mqtools/bin/mqsub --file %s'%mqsub_file
    command='ssh %s@%s "%s"'%(username,cluster_name, command)
    process=Popen(command, stdout=PIPE, shell=True)
    all_id=[]
    for line in process.stdout:
        #Your job 3270533 ("tcezard_job.41.sh") has been submitted
        id=line.split()[2]
        if id.isdigit():
            all_id.append(int(id))
    return all_id


def analyse_qstat_output(list_job_id, qstat_file=True):
    username=os.environ['USER']
    if qstat_file:
        command="""ssh %s@apollo "less /var/log/qstat.snapshot | grep '%s'" """%(username,username)
    else:
        command="""ssh %s@apollo "qstat -u %s" """%(username,username)
    (dummy, handle_out) = os.popen4(command, 'r')
    id_to_info={}
    for line in handle_out.readlines():
        sp_line = line.split()
        if sp_line[0].isdigit():
            if int(sp_line[0]) in list_job_id :
                id_to_info[int(sp_line[0])]=(sp_line[3],sp_line[6])
            elif sp_line[0] in list_job_id:
                id_to_info[sp_line[0]]=(sp_line[3],sp_line[6])
    return id_to_info.keys()
    

def launchCommandLocally(command, nice=True, verbose=True):
    """
    Launch the command locally. 
    """
    if nice and not command.startswith('nice '):
        command='nice %s '%(command)
    return _launchCommand(command,verbose)


def launchCommandLocallyWithMemonitor(command, outputDir='.', verbose=True):
    """
    Launch the command locally using memonitor. 
    """
    return  None
    #TODO: Find a way to change The hard-coded path to a dynamic path
    command2="""python /projects/wtsspipeline/programs/code/trunk/src/utils/memonitor.py -c 'nice %s' -o %s"""%(command,outputDir)
    return _launchCommand(command2, verbose)


def createJarCommand(java_bin,memAlloc,jarFile, arguments):
    return "%s -Xmx%dG -jar %s %s"%(java_bin,memAlloc,jarFile, arguments)


def create_qsub_command(name,wtss_env, command):
    qsub_cmd = """#! /bin/sh
        #$ -S /bin/sh
        #$ -e "%s.e"
        #$ -o "%s.o"
        source %s; %s
        """ % (name, name, wtss_env, command)
    return qsub_cmd



def get_output_stream_from_command(command, logger_name=''):
    logging.debug(command)
    process=Popen(command, stdin=PIPE, stdout=PIPE, stderr=PIPE, shell=True)
    logger=Output_Logger(process.stdin, None, process.stderr, logger_name=logger_name)
    logger.setDaemon(True)
    logger.start()
    return process.stdout, process

def get_input_stream_from_command(command, logger_name=''):
    logging.debug(command)
    process=Popen(command, stdin=PIPE, stdout=PIPE, stderr=PIPE, shell=True)
    logger=Output_Logger(None, process.stdout, process.stderr, logger_name=logger_name)
    logger.setDaemon(True)
    logger.start()
    return process.stdin,process
    

def get_temporary_local_dir():
    import tempfile
    tmp_dir=tempfile.gettempdir()
    return tmp_dir


class Command_runner():
    def __init__(self, print_out=True):
        self.print_out=print_out
    
    def run_command(self,command):
        if self.print_out:
            print command
        else:
            logging.info(command)
            launchCommandLocally(command)
            
            
def shellize_file_name(string):
    s=string.replace(' ','_')
    s=s.replace('/','_')
    s=s.replace('(','')
    s=s.replace(')','')
    s=s.replace("'",'')
    return s