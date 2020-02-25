RADmapper
=========

Set of scripts to create denovo consensus and map the read back

|WARNING: These scripts were not meant to be use outside of their original devlopment environment. They might still contain hard coded path. There are also much better approach nowdays.|
|---|


Download the scripts
--------------------

`git clone https://github.com/tcezard/RADmapper.git`


Prerequesite
------------
You'll need to have python 2.7 installed and set the PYTHONPATH
`export PYTHONPATH=<pathtoRADmapper>/lib:<pathtoRADmapper>/bin`

Set PIPELINE_CONFIG to point to config file

`export PIPELINE_CONFIG=<path to config file>`

Config file shoulg contain the path to executable similar to https://github.com/tcezard/RADmapper/blob/master/pipeline_tool.config

Scripts
-------

**RAD_assign_reads_to_consensus.py**

Align fastq file from Read1 to consensus sequence defined using stacks. Read2 are not align but are assigned the same consensus id and position as their read1 in the bam files. Read2 are also marked as unmapped using the sam flag `0x4`

`python RAD_assign_reads_to_consensus.py -g <consensus_fasta> -1 <first fastq file> -2 <second fastq file> -n <sample_name> -r RG\tID:uid\tSM:sample\tPL:Illumina`


**assemble_read2.py**

Exctract unmapped read2 from bam files and group them per consensus. Output the read2 as fastq files independend subdirectories and assemble them each independently.

`python assemble_read2.py --bam_files <sample1.bam> <sample2.bam> --assembler_name idba_ud --output_dir <output> --read1_consensus_file <consensus_fasta>`


**RAD_coverage_read1.py**

Collect the number of forward read aligned to each consensus sequence

`python RAD_coverage_read1.py -b <sample1.bam> <sample2.bam> -o coverage_read1.tsv`
