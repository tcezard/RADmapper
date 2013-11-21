'''
Created on Jul 14, 2009

@author: tcezard
'''
from utils.config import ConfigReader
import os
import logging
from utils import utils_logging


class Config_file_error(Exception):
    """Raised when an error because a config file is not found."""
    pass


class Pipeline_config():
    def __init__(self, config):
        """Create the Pipeline_param object. that allow access to the parameter in this config file"""
        self.pipeline_config=config
    
    def _get_file_location(self, attribute, key):
        """This function find in the config file the section corresponding to the attribute 
and the key corresponding to the parameter. It takes care of throwing exception if something goes wrong."""
        value = self._get_value(attribute, key)
        if value.startswith('~'):
            value = os.path.join(os.path.expanduser(value))
        if not os.path.exists(value):
            raise Config_file_error('File not found: %s'%value)
        return value
    
    def _get_value(self, attribute, key):
        if not hasattr(self.pipeline_config, attribute):
            raise Config_file_error("section %s does not exist in %s"%(attribute, self.pipeline_config.config_path))
        section = getattr(self.pipeline_config, attribute)
        if not section.has_key(key):
            raise Config_file_error("No parameter %s in section %s in %s"%(key, attribute, self.pipeline_config.config_path))
        return section.get(key)
    
    def _get_array_values(self, attribute, key):
        value = self._get_value(attribute, key)
        array_value = value.split(';')
        return array_value
    

class Pipeline_param(Pipeline_config):
    
    def get_bfast_dir(self):
        return self._get_file_location('bfast', 'bfast_dir')
    
    def get_bwa_dir(self):
        return self._get_file_location('bwa', 'bwa_dir')
    
    def get_gsnap_dir(self):
        return self._get_file_location('gsnap', 'gsnap_dir')
    
    def get_maq_dir(self):
        return self._get_file_location('maq', 'maq_dir')
    
    def get_samtools_dir(self):
        return self._get_file_location('samtools', 'samtools_dir')
    
    def get_picard_dir(self):
        return self._get_file_location('picard', 'picard_dir')
    
    def get_gatk_dir(self):
        return self._get_file_location('gatk', 'gatk_dir')
    
    def get_fastx_dir(self):
        return self._get_file_location('fastx', 'fastx_dir')
    
    def get_solexaqa_dir(self):
        return self._get_file_location('qc_tools', 'solexaqa_dir')
    
    def get_fastqc_dir(self):
        return self._get_file_location('qc_tools', 'fastqc_dir')
    
    def get_matrix2png_dir(self):
        return self._get_file_location('qc_tools', 'matrix2png_dir')
    
    def get_fasta2fastq_bin(self):
        return self._get_file_location('qc_tools', 'fasta2fastq')
    
    def get_blast_dir(self):
        return self._get_file_location('blast', 'blast_dir')
    
    def get_pipeline_place(self):
        return self._get_file_location('pipeline_param', 'pipeline_places')
    
    
class Pipeline_places(Pipeline_config):
    
    def get_illumina_adapter_sequences(self):
        return self._get_file_location('raw_data', 'illumina_adapter_sequences')
    
    def get_illumina_run_dir(self):
        return self._get_file_location('raw_data', 'illumina_machine_run_dir')
    
    def get_seq_archive(self):
        return self._get_array_values('raw_data', 'seq_archive')
    
    def get_illumina_project_dir(self):
        return self._get_file_location('analysed_data','illumina_project_dir')
    
    def get_webserver_run_info(self):
        return self._get_value('analysed_data','webserver_run_info')

    def get_webserver_dir(self):
        return self._get_value('delivered_data','webserver_dir')

    
class WTSS_param():
    
    def __init__(self, wtss_config, wtss_resources=None):
        """Create the WTSS_param object.
If the wtss_resources is not given it's extracted from the wtss_config"""
        self.wtss_config=wtss_config
        if wtss_resources:
            self.wtss_resources=wtss_resources
        else:
            try : 
                self.wtss_resources=ConfigReader(file=self.wtss_config.data_paths.get('resource_conf'))
            except IOError, e:
                logging.critical("Can't find the resources configuration file at %s: exit"%self.wtssConfig.data_paths.get('resource_conf'))

    def get_maq_dir(self):
        return self.wtss_resources.pileup_based_rp.get('maq_dir')
    def get_script_dir(self):
        return self.wtss_resources.pileup_based_rp.get('script_dir')
    def get_samtools_dir(self):
        return self.wtss_resources.pileup_based_rp.get('samtools_dir')
    def get_snvmix2_dir(self):
        return self.wtss_resources.pileup_based_rp.get('snvmix2_dir')
    def get_snvmix2_filter_type(self):
        return self.wtss_resources.pileup_based_rp.get('snvmix2_filter_type')
    def get_snvmix2_base_qual(self):
        return self.wtss_resources.pileup_based_rp.get('snvmix2_base_qual')
    def get_snvmix2_map_qual(self):
        return self.wtss_resources.pileup_based_rp.get('snvmix2_map_qual')
    def get_snvmix2_prob(self):
        return self.wtss_resources.pileup_based_rp.get('snvmix2_prob')
    def get_bfa_dir(self,sp_code):
        return self.wtss_resources.pileup_based_rp.get(sp_code+'_bfa_dir')
    def get_bfa_independent_dir(self, sp_code):
        return self.wtss_resources.pileup_based_rp.get(sp_code+'_bfa_independent_dir')
    def get_index_file(self,sp_code):
        return self.wtss_resources.pileup_based_rp.get(sp_code+'_junction_index_file')
    def get_ensembl_dir(self, sp_code):
        return self.wtss_resources.pileup_based_rp.get(sp_code+'_ensembl_dir')
    def get_intron_ensembl_dir(self, sp_code):
        return self.wtss_resources.pileup_based_rp.get(sp_code+'_intron_ensembl_dir')
    def get_genome_fasta(self, sp_code):
        return self.wtss_resources.pileup_based_rp.get(sp_code+'_genome_fasta')
    def get_known_snps(self, sp_code):
        return self.wtss_resources.pileup_based_rp.get(sp_code+'_known_snps')
    def get_list_contig_file(self, sp_code):
        return self.wtss_resources.pileup_based_rp.get(sp_code+'_genome_fai')
    def get_database_name(self, sp_code):
        return self.wtss_resources.pileup_based_rp.get(sp_code+'_database_name')
    def get_available_species_str(self):
        return self.wtss_resources.pileup_based_rp.get('available_species')
    def get_available_species(self):
        available_species=self.get_available_species_str()
        species_names={}
        if available_species:
            all_species=available_species.split(',')
            for species in all_species:
                all_names=species.split(":")
                if len(all_names)>0:
                    species_names[all_names[0]]=all_names
        if len(species_names)>0:
            return species_names
        else:
            return None
    def get_coding_regions(self,sp_code):
        return self.wtss_resources.useful_files.get(sp_code+'_coding_regions')
    def get_standard_chromosomes(self, sp_code):
        return self.wtss_resources.useful_files.get(sp_code+'_standard_chromosomes')
    def get_translate_chr_names(self, sp_code):
        return self.wtss_resources.useful_files.get(sp_code+'_translate_chr_names')
        
    def get_fp_dir(self):
        return self.wtss_resources.findpeaks.get('fp_dir')
    def get_fp_qualityfilter(self):
        return self.wtss_resources.findpeaks.get('qualityfilter')
    def get_fp_pet_flag_se(self):
        return self.wtss_resources.findpeaks.get('pet_flag_se')
    def get_fp_pet_flag_pe(self):
        return self.wtss_resources.findpeaks.get('pet_flag_pe')

    def get_input_solid_dir(self):
        return self.wtss_resources.solid_param.get('solid_dir')

    def get_env_name(self):
        return self.wtss_config.required_config_env.get('name')
    def get_env_shortname(self):
        return self.wtss_config.required_config_env.get('shortname')
    def get_env_version(self):
        return self.wtss_config.required_config_env.get('version')
        
    def get_bio_apps_analysis(self):
        return self.wtss_config.data_paths.get('bio_apps_analysis')
    def get_wtss_analysis(self):
        return self.wtss_config.data_paths.get('wtss_analysis')
    def get_resource_conf(self):
        return self.wtss_config.data_paths.get('resource_conf')
    def get_mysql_host(self):
        return self.wtss_config.data_paths.get('mysql_host')
    def get_mysql_user(self):
        return self.wtss_config.data_paths.get('mysql_user')
    def get_mysql_password(self):
        return self.wtss_config.data_paths.get('mysql_password')
    
    
def test_is_dir(dir, name):
    test_pass=True
    if dir is not None and os.path.isdir(dir):
        pass
    else:
        logging.warning("%s dir %s is not a directory"%(name,dir))
        test_pass=False
    return test_pass

def test_is_file(file, name):
    test_pass=True
    if file is not None and os.path.isfile(file):
        pass
    else:
        logging.warning("%s: %s is not a file"%(name,file))
        test_pass=False
    return test_pass

def test_is_valid_database(db, name):
    import utils_db
    test_pass=True
    database_names=utils_db.get_databases_name()
    if db in database_names:
        pass
    else:
        logging.warning("%s: %s is not a valid database"%(name,db))
        test_pass=False
    return test_pass

def test_is_int(value, name):
    test_pass=True
    if value is not None and str(value).isdigit():
        pass
    else:
        logging.warning("%s: %s is not an integer"%(name,value))
        test_pass=False
    return test_pass

def test_is_float(value, name):
    test_pass=True
    try:
        float(value)
    except ValueError:
        logging.warning("%s: %s is not a float"%(name,value))
        test_pass=False
    return test_pass

def test_is_not_None(value, name):
    test_pass=True
    if value is not None:
        pass
    else:
        logging.warning("%s: %s is None"%(name,value))
        test_pass=False
    return test_pass

def test_is_comma_separated_list_int(value,name):
    test_pass=True
    try:
        values=value.split(',')
        for v in values:
            int(v)
    except Exception:
        logging.warning("%s: %s is not a comma separated list of int"%(name,value))
        test_pass=False

def test_all():
    import utils_param
    utils_logging.init_logging()
    test_pass=True
    wtss_param=utils_param.get_wtss_parameter()
    test_pass and test_is_dir(wtss_param.get_maq_dir(),"maq")
    test_pass and test_is_dir(wtss_param.get_script_dir(),"script")
    test_pass and test_is_dir(wtss_param.get_samtools_dir(),"samtools")
    test_pass and test_is_dir(wtss_param.get_snvmix2_dir(),"snvmix2")
    test_pass and test_is_not_None(wtss_param.get_snvmix2_filter_type(),"snvmix2 filter type")
    test_pass and test_is_int(wtss_param.get_snvmix2_base_qual(),"snvmix2 base quality")
    test_pass and test_is_int(wtss_param.get_snvmix2_map_qual(),"snvmix2 mapping quality")
    test_pass and test_is_float(wtss_param.get_snvmix2_prob(),"snvmix2 probability")
    test_pass and test_is_dir(wtss_param.get_fp_dir(),"FindPeak")
    test_pass and test_is_int(wtss_param.get_fp_qualityfilter(),"FindPeak quality filter")
    test_pass and test_is_comma_separated_list_int(wtss_param.get_fp_pet_flag_se(),"FindPeak pet flag for single end")
    test_pass and test_is_comma_separated_list_int(wtss_param.get_fp_pet_flag_pe(),"FindPeak pet flag for paired end")
    if test_is_not_None(wtss_param.get_available_species(), "available species"):
        for sp_code in wtss_param.get_available_species().keys():
            test_pass and test_is_dir(wtss_param.get_bfa_dir(sp_code),"bfa dir for %s"%sp_code)
            test_pass and test_is_dir(wtss_param.get_bfa_independent_dir(sp_code),"bfa independent dir for %s"%sp_code)
            test_pass and test_is_file(wtss_param.get_index_file(sp_code),"index file for %s"%sp_code)
            test_pass and test_is_dir(wtss_param.get_ensembl_dir(sp_code),"ensembl dir for %s"%sp_code)
            test_pass and test_is_dir(wtss_param.get_intron_ensembl_dir(sp_code),"intron ensembl dir for %s"%sp_code)
            test_pass and test_is_file(wtss_param.get_known_snps(sp_code),"known snps file for %s"%sp_code)
            test_pass and test_is_file(wtss_param.get_list_contig_file(sp_code),"list of contig file for %s"%sp_code)
            test_pass and test_is_valid_database(wtss_param.get_database_name(sp_code),"database for %s"%sp_code)
            test_pass and test_is_file(wtss_param.get_coding_regions(sp_code),"coding region file for %s"%sp_code)
            test_pass and test_is_file(wtss_param.get_standard_chromosomes(sp_code),"standard chromosomes file for %s"%sp_code)
            test_pass and test_is_file(wtss_param.get_translate_chr_names(sp_code),"translate chr names file for %s"%sp_code)
            
        
    
        wtss_param.get_env_name()
        wtss_param.get_env_shortname()
        wtss_param.get_env_version()
            
        wtss_param.get_bio_apps_analysis()
        wtss_param.get_wtss_analysis()
        wtss_param.get_resource_conf()
        wtss_param.get_mysql_host()
        wtss_param.get_mysql_user()
        wtss_param.get_mysql_password()
        #wtss_param.get_input_solid_dir()
    
if __name__=="__main__":
    test_all()