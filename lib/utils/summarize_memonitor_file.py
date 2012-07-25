#!/usr/bin/env python
import os, sys, time


def summmarizeMemonitorFile(file):
    openFile=open(file)
    maxVmSize=0
    maxVmRSS=0
    maxNbProc=0
    outputLines=[]
    line=openFile.readline()
    #first line contain  the command 
    outputLines.append(line.strip()[1:])
   
    #remove the header lines
    line=openFile.readline()
    line=openFile.readline()
    
    while line:
        splitted=line.split()
        elapsed=int(splitted[0])
        nbProc=int(splitted[1])
        vmSize=int(splitted[2])
        vmRSS=int(splitted[3])
        if vmSize>maxVmSize :
            maxVmSize=vmSize
        if vmRSS>maxVmRSS :
            maxVmRSS=vmRSS
        if nbProc>maxNbProc:
            maxNbProc=nbProc
        line=openFile.readline()
            
    #print os.path.basename(file)
    outputLines.append('Maximum number of process=%s'%maxNbProc)
    outputLines.append('Maximum VM Size=%sMbi'%(maxVmSize/1000))
    outputLines.append('Maximum VM RSS=%sMbi'%(maxVmRSS/1000))
    outputLines.append('Total Elapsed time=%s:%s'%(elapsed/60,elapsed%60))
    outputLines.append('')
    return '\n'.join(outputLines)
        
def main(files):
    fileToOutput={}
    for file in files:
        if os.path.exists(file):
            if os.path.isfile(file):
                fileToOutput[os.stat(file)[8]]=summmarizeMemonitorFile(file)
            elif os.path.isdir(file):
                listFile=os.listdir(file)
                for file in listFile:
                    if file.find('memonitor'):
                        fileToOutput[os.stat(file)[8]]=summmarizeMemonitorFile(file)
    keys=fileToOutput.keys()
    keys.sort() ## sort by the last modification time
    for key in keys:
        print fileToOutput.get(key)
            

if __name__=="__main__":
    main(sys.argv[1:])
