#!/usr/bin/env python
import os, sys, time
from optparse import OptionParser


def getVm_Size(pid):
    if os.path.exists("/proc/"+str(pid)+"/status"):
        lines=[line
                for line in open("/proc/"+str(pid)+"/status").readlines()
                    if line.find("VmSize")!=-1]
        return sum([int(line.split()[1]) for line in lines])
    else:
        return 0

def getVm_RSS(pid):
    
    if os.path.exists("/proc/"+str(pid)+"/status"):
        #for line in open("/proc/"+str(pid)+"/status").readlines():
        #    print line
        lines=[line
                for line in open("/proc/"+str(pid)+"/status").readlines()
                    if line.find("VmRSS")!=-1]
        return sum([int(line.split()[1]) for line in lines])
    else:
        return 0

def main(command, memOutputDir='',interval=5.0):
    if len(command)>10:
        abrev=command[0:10].replace(' ', '_')
    else :
        abrev=command.replace(' ', '_')
    interval=float(interval)
    pid= os.fork()
    if pid==0:
        #print "child"
        #print command
        for line in os.popen(command).readlines():
            print line.strip()
        
    elif pid:
        startAt=time.time()
        logFile='%smemonitor_%s.csv'%(memOutputDir,pid)
        openfile=open(logFile,'w');
        openfile.write('#%s\n'%command)
        openfile.write('#Elapsed_time\tNb_proc\tVm_Size\tVm_RSS\n')
        print "monitor command %s in %s"%(command,logFile);
        while os.waitpid(pid,os.WNOHANG)==(0,0):
            (allPID,Vm_Size,VM_RSS)=getChildrenMemUsage(pid)
            #print allPID
            if (Vm_Size==0 and VM_RSS==0):
                print "status file not found for any of the children"
            Elapsed_time=(time.time()-startAt)
            openfile.write('%s\t%s\t%s\t%s\n'%(int(Elapsed_time),len(allPID),Vm_Size,VM_RSS))
            openfile.flush()
            time.sleep(interval)
        openfile.close()
    else:
        print "Error Something broke: can't fork()!"
        return -1
    
def getChildrenMemUsage(pid):
    """ get the memory usage of the pid's children.
    @return an array of string containing all the children's PID, 
    the total Vm Size for all those process and the total Vm RSS for those process. 
    """
    sum_Vm_Size=0
    sum_VM_RSS=0
    allPID=[]
    #allPID.append(pid)
    for p in getChildrenPid(pid): allPID.append(p)
    for p in allPID:
        sum_Vm_Size+= getVm_Size(p)
        sum_VM_RSS += getVm_RSS(p)
    
    return (allPID,sum_Vm_Size,sum_VM_RSS)

def getChildrenPid(pid):
    """
    Recursive method that will get all the children, grand children ... of the specified pid 
    @param pid: The parent pid of the children your are looking for
    @return array of string containing the children pid 
    """
    command="ps x -o pid=,ppid= | sed -r 's/ +/\\t/g '| awk '{if ($2==%s) print $1}'"%pid
    childrenPID=[]
    for line in os.popen(command).readlines():
        #print line
        childrenPID.append(line.strip())
    if len(childrenPID)>0:
        tmpPID=[]
        for pid in childrenPID:
            tmp=getChildrenPid(pid);
            for pid in tmp:
                tmpPID.append(pid)
        for pid in tmpPID:
            childrenPID.append(pid)
    return childrenPID

def _prepare_optparser():
    """Prepare optparser object. New options will be added in this
    function first.
    """
    usage = """usage: %prog <-c "double_quoted_command"> <-o output_dirercotry>"""
    description = """Memonitor is a benchmarking script that follow the memory usage of a command.
    The command is launch in a shell and memonitor follow the process and its children."""

    optparser = OptionParser(version="%prog 0.1.2",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="show this help message and exit.")
    optparser.add_option("-c","--command",dest="command",type="string",
                         help="The command that will be launched and followed. It must be surrounded by quote or double quote.")
    optparser.add_option("-o","--output",dest="outputDir",type="string",
                         help="Path to a the directory where the memonitor file will be output.",
                         default='')
    optparser.add_option("-i","--interval",dest="interval",type="float",default=5.0,
                         help="The interval in between each memory call.")
    
    return optparser

if __name__=="__main__":
    optparser=_prepare_optparser()
    (options,args) = optparser.parse_args()
    
    main(options.command,options.outputDir,options.interval)
    
