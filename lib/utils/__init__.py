"""Package that contains utility method and object
@author: T. Cezard
@since: Oct 2 2008
"""
import os, sys, logging
import glob, re
from utils import utils_param, utils_logging, utils_commands
from utils.utils_param import Resource_File_Exception
import command_runner
from utils.parameters import Config_file_error


def pass_funct(*arg, **kwargs):
    """function that does nothing"""
    return

#######################################################
#        Directories and permissions utilities        #
#######################################################


def createDirectories(baseDir,directories, server=None):
    """
    This function create the directories listed in the directories array.
    @param baseDir: the parent directory.
    @param directories: the list of directory to create.
    @param server: the server on which the directory will be created.
    """
    for (directory) in directories:
        dir=os.path.join(baseDir,directory)
        if server is not None:
            command = 'ssh %s "ls -d1 %s"'%(server, dir)
            stream, process = utils_commands.get_output_stream_from_command(command, logger_name=None)
            line=None
            for line in stream:
                line=line.strip()
                if line == dir:
                    break
            if line != dir:
                logging.info('%s does not exists on %s: create it'%(dir, server))
                command = 'ssh %s "mkdir %s"'%(server, dir)
                utils_commands.launchCommandLocally(command)
        elif not os.path.exists(dir):
            logging.info('%s does not exists: create it'%dir)
            os.mkdir(dir,0775)


def chmodAndGroup(listPath, file_mode=0664, dir_mode=0775, group_id=548, quiet=False):
    """Change the permission and group of the given files and directories.
    The given integer value must be octal integer.
    To change an integer into an octal int apply int(value,8)."""
    from exceptions import OSError 
    umask=os.umask(0)
    for path in listPath:
        try:
            if os.path.isfile(path):
                os.chmod(path, file_mode)
            elif os.path.isdir(path):
                os.chmod(path, dir_mode)
        except OSError, e:
            if not quiet:
                logging.error('OSError chmod [%d]: %s on %s'% (e.args[0], e.args[1], path))
        try:
            os.chown(path, -1,group_id)
        except OSError, e:
            if not quiet:
                logging.error('OSError chown [%d]: %s on %s'% (e.args[0], e.args[1], path))
    os.umask(umask)


def change_mode_and_group_for_dir_recur(directory, mode_for_file=0664, mode_for_dir=0775,
                                        group_id=548, quiet=False):
    list_file=os.listdir(directory)
    list_path=[directory]
    for file in list_file:
        file=os.path.join(directory,file)
        list_path.append(file)
        if os.path.isdir(file):
            change_mode_and_group_for_dir_recur(file, mode_for_file, mode_for_dir, group_id, quiet)
    chmodAndGroup(list_path, mode_for_file, mode_for_dir, group_id,quiet)


###############################
#        Maq utilities        #
###############################


def getMaqFirstLine(maqBin,mapfile):
    """Get the first line of a map file using mapview"""
    from subprocess import Popen
    from subprocess import PIPE 
    process=Popen('%s mapview %s'%(maqBin,mapfile), stdin=PIPE, stdout=PIPE, stderr=PIPE,
                         shell=True)
    child_stderr=process.stderr
    child_stdout=process.stdout
    line=child_stdout.readline()
    os.system('kill -9 %s'%process.pid)
    errLine=child_stderr.readline()
    for line in child_stdout:
        pass
    for line in child_stderr:
        pass
    process.wait()
    if errLine.find('cannot execute binary file')>0:
        raise os.error('Cannot execute binary with %s on this system'%maqBin)
    return line

def findAllMaqBin(maqDir):
    maqBinList=glob.glob(os.path.join(maqDir,'*','maq'))
    maqBinList.extend(glob.glob(os.path.join(maqDir,'maq')))
    valid_maq_bin_list=[]
    for maq_bin in maqBinList:
        if os.path.exists(maq_bin):
            valid_maq_bin_list.append(maq_bin)
    return valid_maq_bin_list

def findNewestMaqVersion(maqDir):
    maqBinList=findAllMaqBin(maqDir)
    maqBinList.sort(reverse=True)
    return maqBinList[0]
    

def findNewestMaqVersionForMapfile(maqDir, mapfile):
    validMaqBin=None
    maqBinList=findAllMaqBin(maqDir)
    maqBinList.sort(reverse=True)
    for maqBin in maqBinList:
        line=getMaqFirstLine(maqBin, mapfile)
        if line and len(line.strip())>0:
            validMaqBin=maqBin
            break
    return validMaqBin

def get_maq2sam_bin(samtools_dir, maq_version):
    """Locate the maq2sam binnary in the samtools dir. returns None if enable to find it."""
    if maq_version.find('0.7')>=0:
        maq2sam_bin_name='maq2sam-long'
    else:
        maq2sam_bin_name='maq2sam-short'
    #there are two place where I saw this bin being placed
    maq2sam_binary=os.path.join(samtools_dir,maq2sam_bin_name)
    if not os.path.exists(maq2sam_binary):
        maq2sam_binary=os.path.join(samtools_dir,'misc',maq2sam_bin_name)
    if not os.path.exists(maq2sam_binary):
        return None
    else:
        return maq2sam_binary


def get_mapview_stream(maq_bin, map_file):
    """This method opens a .map file with Maq and returns an open file. The std error will be output in the console through another thread."""
    command='%s mapview %s'%(maq_bin, map_file)
    stdout, process = utils_commands.get_output_stream_from_command(command)
    return stdout

def get_sam_stream(bam_file, samtools_bin=None, options='', chomosome_and_position=''):
    """This method opens a .bam file with samtools and returns an open file. The std error will be output in the console through another thread."""
    if samtools_bin==None:
        try:
            pipeline_parm=utils_param.get_pipeline_parameters()
            samtools_bin=os.path.join(pipeline_parm.get_samtools_dir(),'samtools')
        except Config_file_error, e:
            logging.warning("Can't find the configuration file you'll need to have samtools in you path.")
            samtools_bin='samtools'
    
    command='%s view %s %s %s'%(samtools_bin, options, bam_file, chomosome_and_position)
    stdout, process = utils_commands.get_output_stream_from_command(command)
    return stdout, process


def get_pileup_from_bam(bam_file, genome_file=None, samtools_bin=None, options=''):
    if samtools_bin==None:
        try:
            pipeline_parm=utils_param.get_pipeline_parameters()
            samtools_bin=os.path.join(pipeline_parm.get_samtools_dir(),'samtools')
        except Config_file_error, e:
            logging.warning("Can't find the configuration file you'll need to have samtools in you path.")
            samtools_bin='samtools'
    if bam_file=='PIPE':
        bam_file='-'
    if genome_file:
        command = '%s pileup -f %s %s %s'%(samtools_bin, genome_file, bam_file, options)
    else:
        command = '%s pileup %s %s'%(samtools_bin, bam_file, options)
    logging.debug(command)
    stream, process = utils_commands.get_output_stream_from_command(command)
    return stream

def get_mpileup_from_bam(bam_file, genome_file=None, samtools_bin=None, options=''):
    if samtools_bin==None:
        try:
            pipeline_parm=utils_param.get_pipeline_parameters()
            samtools_bin=os.path.join(pipeline_parm.get_samtools_dir(),'samtools')
        except Config_file_error, e:
            logging.warning("Can't find the configuration file you'll need to have samtools in you path.")
            samtools_bin='samtools'
    if bam_file=='PIPE':
        bam_file='-'
    if genome_file:
        command = '%s mpileup -f %s %s %s'%(samtools_bin, genome_file, bam_file, options)
    else:
        command = '%s mpileup %s %s'%(samtools_bin, bam_file, options)
    logging.debug(command)
    stream, process = utils_commands.get_output_stream_from_command(command)
    return stream

def get_map_contigs(maqBin, mapfile):
    command= '%s mapview %s'%(maqBin,mapfile)
    from subprocess import Popen
    from subprocess import PIPE
    from subprocess import STDOUT
    import subprocess
    process=Popen(command, stdin=PIPE, stdout=PIPE, shell=True)
    line=process.stdout.readline()
    allcontig=[]
    while line:
        contig=line.split()[1]
        if contig not in allcontig:
            allcontig.append(contig)
        line=process.stdout.readline()
    return allcontig


##################################
#        Picard utilities        #
##################################

def create_sequence_dictionary(picard_dir, genome_file):
    """Create a sequence dictionary from the genome file provided"""
    name, dummy=os.path.splitext(genome_file)
    genome_dict=name+'.dict'
    if not os.path.exists(genome_dict) or os.path.getmtime(genome_file) > os.path.getmtime(genome_dict):
        CreateSequenceDictionary_jar=os.path.join(picard_dir,'CreateSequenceDictionary.jar')
        command='java -jar %s REFERENCE=%s O=%s'%(CreateSequenceDictionary_jar,genome_file,genome_dict)
        return_code=command_runner.run_command(command)
        if return_code!=0:
            genome_dict=None
    return genome_dict

def sort_bam_file_per_coordinate(picard_dir, input_bam, output_bam, overwrite=False,validation_stringency="LENIENT",**kwargs):
    return_code=1
    if picard_dir:
        options=[]
        for key in kwargs.keys():
            options.append("%s=%s"%(key,kwargs.get(key)))
        sort_jar=os.path.join(picard_dir,'SortSam.jar')
        command="java -Xmx4G -jar %s I=%s O=%s SO=coordinate VALIDATION_STRINGENCY=%s %s"%(sort_jar, input_bam, output_bam, 
                                                                                        validation_stringency, ' '.join(options))
        
        if (not os.path.exists(output_bam)) or overwrite:
            return_code = command_runner.run_command(command)
        else:
            logging.warning('The file %s exists, use overwrite option to overwrite if applicable.'%output_bam)
    return return_code



####################################
#        samtools utilities        #
####################################

def generate_fai_from_genome(genome_fasta, samtools_bin=None):
    if samtools_bin==None:
        try:
            pipeline_parm=utils_param.get_pipeline_parameters()
            samtools_bin=os.path.join(pipeline_parm.get_samtools_dir(),'samtools')
        except Config_file_error, e:
            logging.warning("Can't find the configuration file you'll need to have samtools in you path.")
            samtools_bin='samtools'
    genome_fai=None
    if genome_fasta:
        genome_fai=genome_fasta+'.fai'
        if not os.path.exists(genome_fai) or os.path.getmtime(genome_fasta) > os.path.getmtime(genome_fai):
            command='%s faidx %s'%(samtools_bin, genome_fasta)
            logging.info(command)
            return_code=command_runner.run_command(command)
            if return_code!=0:
                genome_fai=None
    return genome_fai


def index_bam_file(bam_file, samtools_bin=None):
    if samtools_bin==None:
        try:
            pipeline_parm=utils_param.get_pipeline_parameters()
            samtools_bin=os.path.join(pipeline_parm.get_samtools_dir(),'samtools')
        except Config_file_error, e:
            logging.warning("Can't find the configuration file you'll need to have samtools in you path.")
            samtools_bin='samtools'
    index_bam=bam_file+'.bai'
    if not os.path.exists(index_bam) or os.path.getmtime(bam_file) > os.path.getmtime(index_bam):
        command = '%s index %s'%(samtools_bin, bam_file)
        return_code = command_runner.run_command(command)
        if return_code!=0:
            index_bam=None
    return index_bam


def get_list_contig_file(sp_code,wtss_param=None):
    if not wtss_param:
        wtss_param=utils_param.get_wtss_parameter()
    list_contig_file=wtss_param.get_list_contig_file(sp_code)
    if not list_contig_file:
        samtools_bin=os.path.join(wtss_param.get_samtools_dir(),'samtools')
        genome_fasta=wtss_param.get_genome_fasta(sp_code)
        list_contig_file=generate_fai_from_genome(genome_fasta=genome_fasta, samtools_bin=samtools_bin)
        if not list_contig_file:
            logging.error("Can't create the .fai file from the genome fasta using %s on %s"%(samtools_bin,genome_fasta))
    return list_contig_file


######################################
#         Version utilities          #
######################################


def getMaqVersion(maqBin):
    dummy, child_stdout_stderr=os.popen4(maqBin)
    version=None
    for line in child_stdout_stderr:
        if line.startswith('Version:'):
            version=line.split()[1].strip()
    return version


def getSamtoolsVersion(samtools_bin):
    dummy, child_stdout_stderr=os.popen4(samtools_bin)
    version=None
    for line in child_stdout_stderr:
        if line.startswith('Version:'):
            version=line.split()[1].strip()
    child_stdout_stderr=None
    return version

def get_bwa_version(bwa_bin):
    dummy, child_stdout_stderr=os.popen4(bwa_bin)
    version=None
    for line in child_stdout_stderr:
        if line.startswith('Version:'):
            version=line.split()[1].strip()
    child_stdout_stderr=None
    return version

def getWtss_version():
    prog_version='UNVERSION'
    dir=os.path.dirname(sys.argv[0])
    version_path=os.path.abspath(os.path.join(dir,'version.txt'))
    if not os.path.isfile(version_path):
        version_path=os.path.abspath(dir+'/../version.txt')
    if not os.path.isfile(version_path):
        version_path=os.path.abspath(dir+'/../../version.txt')
    if not os.path.isfile(version_path):
        version_path=os.path.abspath(dir+'/../../../version.txt')
    if os.path.isfile(version_path):
        prog_version=open(version_path).readline().strip()
    return prog_version

def compare_version_number(v1, v2):
    def fixup(i):
        try:
            return int(i)
        except ValueError:
            return i
    a = map(fixup, re.findall("\d+|\w+", v1))
    b = map(fixup, re.findall("\d+|\w+", v2))
    return cmp(a, b) # -1 if a<b, 0 if a=b, 1 if a>b
    

###################################
#         File utilities          #
###################################


def checkFiles(files_path_list):
    """Check if the given set of file are files and if they are not empty
    @return True only if all the file exists and contain something. 
    """
    returnValue=True
    for filePath in files_path_list:
        returnValue=returnValue and checkFile(filePath)
    return returnValue
    
def checkFile(filePath, server=None):
    if check_file_or_dir(filePath, server) == 'file':
        return True
    else:
        return False

def check_file_or_dir(filePath, server=None):
    """ Check if the given file is a file and if its size is greater than 0."""
    if server:
        returnValue = False
        command = 'ssh %s "ls -ld %s"'%(server, filePath)
        stream,process = utils_commands.get_output_stream_from_command(command,logger_name=None)
        line=None
        for line in stream:
            line=line.strip()
            if line:
                if line.startswith('d'):
                    #It's a directory
                    returnValue = 'dir'
                else:
                    returnValue = 'file'
                break
            else:
                returnValue = False
    else:
        if os.path.isfile(filePath):
            returnValue = 'file'
        elif os.path.isdir(filePath):
            returnValue = 'dir'
        else:
            returnValue = False
        
    return returnValue

def check_one_file_true(files_path_list):
    """Check if the given set of file are files and if they are not empty
    @return True one of the files exist and contain something. 
    """
    returnValue=False
    for filePath in files_path_list:
        returnValue=returnValue or checkFile(filePath)
    return returnValue


def copy_file_across(source, destination, server_source=None, server_destination=None, overwrite=False):
    if check_file_or_dir(source,server_source) == 'file' and check_file_or_dir(destination,server_destination) == 'dir':
        destination_file=os.path.join(destination,os.path.basename(source))
    else:
        destination_file=destination
    if not checkFile(destination_file, server_destination) or overwrite:
        if server_source and server_destination:
            name = os.path.basename(source)
            tmp_file = '/tmp/%s'%(name)
            command='scp %s:%s %s'%(server_source, source, tmp_file)
            command_runner.run_command(command)
            command='scp %s %s:%s'%(tmp_file, server_destination, destination)
            command_runner.run_command(command)
            os.remove(tmp_file)
        else:
            if server_source:
                command='scp %s:%s %s'%(server_source, source, destination)
            elif server_destination:
                command='scp %s %s:%s'%(source, server_destination, destination)
            else:
                command='scp %s %s'%(source, destination)
            command_runner.run_command(command)
    else:
        if server_destination:
            logging.warning('%s exist on %s use force to overwrite'%(destination_file, server_destination))
        else:
            logging.warning('%s exist use force to overwrite'%destination_file)


##########################################
#        Resource files utilities        #
##########################################


def get_sp_code(specie, available_species=None):
    if specie:
        sp_lo=specie.lower()
        if not available_species:
            wtss_param=utils_param.get_wtss_parameter()
            available_species=wtss_param.get_available_species()
            if not available_species:
                #keep for backward compatibility
                available_species={'hs':['human','homo_sapiens','hs', '9606'],
                       'mm':['mouse', 'mus_musculus', 'mm', '10090'],
                       'pt':['poplar', 'populus_trichocarpa', 'pt', '3694'],
                       'ce':['worm', 'c_elegans', 'ce', '6239']}
        for sp in available_species.keys():
            if sp_lo in available_species.get(sp) or specie in available_species.get(sp):
                return sp
    return None

def check_sp_code_available(av_sp):
    species_keyword={}
    valid=True
    for sp in av_sp.keys():
        for keyword in av_sp.get(sp):
            other_sp=species_keyword.get(keyword)
            if other_sp and other_sp!=sp:
                logging.error("Keyword %s is in use for %s and %s"%(keyword,other_sp,sp))
                valid=False
            species_keyword[keyword]=sp
    return valid



def Find_analysis_base_tab_in_dir(ensembl_dir,genomic=False, intron=False):
    """
    """
    if intron:
        file='intronic_bases_tab.txt'
    else:
        file='exonic_bases_tab.txt'
    if genomic:
        file='genomic_'+file
    full_path=os.path.join(ensembl_dir,file)
    if os.path.exists(full_path):
        return full_path
    else: 
        logging.error( "Can't find required file in %s"%ensembl_dir)
        raise Resource_File_Exception("%s can't be found"%full_path)


def get_standard_chromosome(standard_chromosomes_file):
    """get the standard chromosome from the a file name. 
    """
    if not standard_chromosomes_file or not os.path.exists(standard_chromosomes_file):
        logging.error("invalid standard chromosome file: %s"%standard_chromosomes_file)
        return None
    CHR=[]
    infile = open(standard_chromosomes_file)
    for l in infile:
        c=l.strip()
        CHR.append(c)
    infile.close()
    return CHR



#################################
#        email utilities        #
#################################


def ask_for_help(mailhost='', reporter_email='', notification_emails=[], subject='', body=''):
    """
    Send out an email with subject and sender if specified.
    """
    import smtplib
    out = []
    if reporter_email: out.append('From: %s'%(reporter_email))
    if subject: out.append('Subject: %s'%(subject))
    out.append(body)
    smtpserver = smtplib.SMTP(mailhost, 25)
    smtpserver.sendmail(reporter_email, notification_emails, '\n'.join(out) )
    smtpserver.quit()
    return


## this function  sort dict by value
def sort_by_value(dict):
    items = [(v, k) for k, v in dict.items()]
    items.sort(reverse=True)
    items = [(k, v) for v, k in items]
    return items



def get_list_int_from_str(int_list_str):
    list_str=[]
    tmp=int_list_str.split()
    for i in range(len(tmp)):
        list_str.extend(tmp[i].strip().split(',')) 
    list_int=[]  
    for str in list_str:
        if str.isdigit():
            list_int.append(int(str))
        elif len(str.split('-'))==2:
            a,b= str.split('-')
            if a.isdigit() and b.isdigit() and int(a)<=int(b):
                list_int.extend(range(int(a),int(b)+1))
    return list_int

#####################################
# http://code.activestate.com/recipes/498181-add-thousands-separator-commas-to-formatted-number/
# Code from Michael Robellard's comment made 28 Feb 2010
# Modified for leading +, -, space on 1 Mar 2010 by Glenn Linderman
# 
# Tail recursion removed and  leading garbage handled on March 12 2010, Alessandro Forghieri
def split_thousands( s, tSep=',', dSep='.'):
    '''Splits a general float on thousands. GIGO on general input'''
    if s == None:
        return 0
    if not isinstance( s, str ):
        s = str( s )

    cnt=0
    numChars=dSep+'0123456789'
    ls=len(s)
    while cnt < ls and s[cnt] not in numChars: cnt += 1

    lhs = s[ 0:cnt ]
    s = s[ cnt: ]
    if dSep == '':
        cnt = -1
    else:
        cnt = s.rfind( dSep )
    if cnt > 0:
        rhs = dSep + s[ cnt+1: ]
        s = s[ :cnt ]
    else:
        rhs = ''

    splt=''
    while s != '':
        splt= s[ -3: ] + tSep + splt
        s = s[ :-3 ]

    return lhs + splt[ :-1 ] + rhs

def longest_common_substr_from_start(s1,s2):
    for pos in range(len(s1)):
        if len(s2)<=pos:
            break
        elif s1[pos]==s2[pos]:
            pass
        else:
            break
    return s1[:pos] 

   
   
_code2atgc={'A':['A'],'a':['A'],'T':['T'],'t':['T'],'G':['G'],'g':['G'],'C':['C'],'c':['C'],
           'R':['A','G'],'r':['A','G'],'Y':['C','T'],'y':['C','T'],'M':['C','A'],'m':['C','A'],
           'K':['T','G'],'k':['T','G'],'W':['T','A'],'w':['T','A'],'S':['C','G'],'s':['C','G'],
           'B':['C','T','G'],'b':['C','T','G'],'D':['A','T','G'],'d':['A','T','G'],'H':['A','T','C'],
           'h':['A','T','C'],'V':['A','C','G'],'v':['A','C','G'],'N':['A','C','T','G'],'n':['A','C','T','G']}


_atgc2iupac={'A':'A','T':'T','G':'G','C':'C','AG':'R','CT':'Y','AC':'M',
           'GT':'K','AT':'W','CG':'S','CGT':'B','AGT':'D','ACT':'H',
           'ACG':'V','ACGT':'N','R':'R','Y':'Y','M':'M','K':'K','W':'W',
           'S':'S','B':'B','D':'D','H':'H','V':'V','N':'N' }

def get_iupac_alphabet():
    return ''.join(_code2atgc.keys())

def get_nt_array_from_IUPAC(IUPAC_code):
    array=_code2atgc.get(IUPAC_code)
    return copy.copy(array)

def is_IUPAC_code(IUPAC_code):
    return _code2atgc.has_key(IUPAC_code)

def atgc2iupac(array_of_code):
    array_of_code.sort()
    s=''.join(array_of_code).upper()
    return _atgc2iupac.get(s)

def cplx_atgc2iupac(array_of_code):
    atgc={}
    for code in array_of_code:
        array=_code2atgc.get(code)
        for c in array: atgc[c]=1
    array=atgc.keys()
    array.sort()
    s=''.join(array).upper()
    return _atgc2iupac.get(s)


if __name__=="1__main__":
    maqDir='/home/pubseq/BioSw/Maq/'
    mapFile='/archive/solexa1_4/analysis/HS0615/30CG4AAXX_4/maq/30CG4AAXX_4.map'
    mapFile='/archive/solexa1_4/analysis/HS0616/30CHEAAXX_6/maq/30CHEAAXX_6.map'
    validMaqBin=findNewestMaqVersionForMapfile(maqDir, mapFile )
    print validMaqBin



if __name__=="1__main__":
    mapfile='/archive/solexa1_4/analysis/HS0727/309YUAAXX_7/maq/309YUAAXX_7.map'
    maqDir='/home/pubseq/BioSw/Maq/'
    maqBin=findNewestMaqVersionForMapfile(maqDir, mapfile)
    
    map_contig= get_map_contigs(maqBin, mapfile)
    for contig in map_contig:
        print contig
        

if __name__=="1__main__":
    mailhost='smtp.staffmail.ed.ac.uk'
    reporter_email='tcezard@staffmail.ed.ac.uk'
    notification_emails=['tcezard@staffmail.ed.ac.uk']
    subject='Test email'
    body='Test body'
    ask_for_help(mailhost, reporter_email, notification_emails,subject, body)
    
if __name__=="__main__":
    print compare_version_number("0.5.9", "0.5.9")
    print compare_version_number("0.5.8", "0.5.9")
    print compare_version_number("0.5.9", "0.5.8")
    print compare_version_number("0.5.8-r18", "0.5.9-r16")
