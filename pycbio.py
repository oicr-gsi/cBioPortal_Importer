# -*- coding: utf-8 -*-
"""
Created on Fri Sep 24 12:26:55 2021

@author: rjovelin
"""


import configparser
import argparse
import subprocess
import os
import shutil
import gzip
import sys
import uuid





def extract_files_from_map(mapfile, data_type):
    '''
    (str, str) -> list
    
    Returns a list of MAF files from the mapping file
        
    Parameters
    ----------
    - mapfile (str): Mapping file (map.csv) that contains paths to maf, seg, gep and mavis files    
    - data_type (str): File type to link in their respective folders.
                       Accepted values: maf, gep, fus, and seg
    '''

    # create input directories for each file type from map file    
    infile = open(mapfile)
    content = infile.read().rstrip().split('\n')
    infile.close()
    for i in range(len(content)):
        content[i] = list(map(lambda x: x.strip(), content[i].split(',')))
    
    # make a list of samples and files for which files exist
    if data_type == 'maf':
        # get the maf files
        j = 2
    elif data_type == 'seg':
        # get the cna files
        j = 3
    elif data_type == 'gep':
        # get the rna files
        j = 4
    elif data_type == 'fus':
        # get the fusion files
        j = 5
    
    files = [i[j] for i in content if i[j].upper() != 'NA' and os.path.isfile(i[j])]
    return files


def create_input_directories(outdir, mapfile):
    '''
    (str, str) -> None
    
    Create sub-directories in outdir for each type of file listed in map file if at least 1 such file exists
        
    Parameters
    ----------
    - outdir (str): Path to the output directory where mafdir, sgedir, fusdir and gepdir folders are located 
    - mapfile (str): Mapping file (map.csv) that contains paths to maf, seg, gep and mavis files
    '''
    
    for i in ['maf', 'seg', 'gep', 'fus']:
        files = extract_files_from_map(mapfile, i)
        if files:
            filedir = os.path.join(outdir, '{0}dir'.format(i))
            os.makedirs(filedir, exist_ok=True)
    

def write_meta_study(outputfile, study, genome, cancerType):
    '''
    (str, str, str, str) -> None
    
    Write file meta_study.txt in the cbioportal_import_data folder
    
    Parameters
    ----------
    - outputfile (str): Path to the outputfile 
    - study (str): Name of the study as it appears in cBioPortal
    - genome (str): Reference genome (hg19 or hg38)
    - cancerType (str): Cancer type as defined in http://oncotree.mskcc.org
    '''

    newfile = open(outputfile, 'w')
    L = ['cancer_study_identifier: {0}'.format(study),
         'description: {0}'.format(study),
         'groups: ',
         'name: {0}'.format(study),
         'reference_genome: {0}'.format(genome),
         'short_name: {0}'.format(study),
         'add_global_case_list: true',
         'type_of_cancer: {0}'.format(cancerType)]
    newfile.write('\n'.join(L))
    newfile.close()     
    

def write_meta_clinical(cbio_import_dir, study, data_type):
    '''
    (str, str, str) -> None
    
    Write file meta_clinical_patients.txt or meta_clinical_samples.txt in the cbioportal_import_data folder
    
    Parameters
    ----------
    - cbio_import_dir (str): Path to the cbioportal_import_data directory 
    - study (str): Name of the study as it appears in cBioPortal
    - data_type(str): Sample or patient
    '''
    
    if data_type.lower() not in ['sample', 'patient']:
        raise ValueError('ERROR. Data type must be "sample" or "patient"')
    
    filename = 'meta_clinical_{0}s.txt'.format(data_type.lower())
    outputfile = os.path.join(cbio_import_dir, filename)
    newfile = open(outputfile, 'w')
    L = ['cancer_study_identifier: {0}'.format(study),
         'data_filename: {0}'.format(filename.replace('meta', 'data')),
         'datatype: {0}_ATTRIBUTES'.format(data_type.upper()),
         'genetic_alteration_type: CLINICAL']
    newfile.write('\n'.join(L) + '\n')
    newfile.close()
    

def get_genome_version(mapfile):
    '''
    (str) -> str
    
    Returns the genome version (hg19 or hg38) using the ref of the first MAF file listed in the map file
    
    Parameters
    ----------
    - mapfile (str): Mapping file (map.csv) that contains paths to maf, seg, gep and mavis files    
    '''
    
    # get the 1st maf listed in mapfile
    maf = extract_files_from_map(mapfile, 'maf')[0]
    
    # grab genome column in maf. skipping header and commented lines
    infile = gzip.open(maf, 'rt')
    ncbiGenome = []
    for line in infile:
        if not line.startswith('#') and 'Hugo_Symbol' not in line:
            line = line.rstrip().split('\t')
            ncbiGenome.append(line[3])              
    infile.close()
    # reduce list to genome value
    ncbiGenome = list(set(ncbiGenome))[0]
    # convert genome identifier
    if ncbiGenome == "GRCh38":
        genomev="hg38"
    elif ncbiGenome == "GRCh37":
        genomev="hg19"
    return genomev        



def write_cases(outputfile, study, mapfile, data_type):
    '''
    (str, str, str, str) -> None
    
    Write list of samples profiled for data type 
            
    Parameters
    ----------
    - outputfile (str): Path to the outputfile
    - study (str): Name of the study as it appears in cBioPortal
    - mapfile (str): Mapping file (map.csv) that contains paths to maf, seg, gep and mavis files    
    - data_type (str): The type of data considered.
                       Values accepted are: seq, rna, cna, cna_seq, cna_seq_rna
    
    '''
    
    #read mapfile
    infile = open(mapfile)
    content = infile.read().rstrip().split('\n')
    infile.close()
    
    # make a list of samples
    if data_type == 'seq':
        # make a list of samples for which maf files are available
        samples = [i.split(',')[1] for i in content if i.split(',')[2].upper() != 'NA']
        name = 'Samples profiled for mutations'
        description = 'This is this case list that contains all samples that are profiled for mutations.'
        stable_id = '{0}_sequenced'.format(study)
    elif data_type == 'fusion':
        # make a list of samples for which fusion files are available
        samples = [i.split(',')[1] for i in content if i.split(',')[5].upper() != 'NA']
        name = 'Samples profiled for mutations'
        description = 'This is this case list that contains all samples that are profiled for fusion.'
        stable_id = '{0}_fusion'.format(study)
    elif data_type == 'rna':
        # make a list of samples for which rsem files are available
        samples = [i.split(',')[1] for i in content if i.split(',')[4].upper() != 'NA']
        name = 'Samples profiled for rnaseq'
        description = 'This is this case list that contains all samples that are profiled for rnaseq.'
        stable_id = '{0}_rna_seq_mrna'.format(study)    
    elif data_type == 'cna':
        # make a list of samples for which cna files are available
        samples = [i.split(',')[1] for i in content if i.split(',')[3].upper() != 'NA']
        name = 'Samples profiled for cnas'
        description = 'This is this case list that contains all samples that are profiled for cnas.'
        stable_id = '{0}_cna'.format(study)
    elif data_type == 'cna_seq':
        # make a list of samples for which cna and maf files are available
        samples = [i.split(',')[1] for i in content if i.split(',')[3].upper() != 'NA' and i.split(',')[2].upper() != 'NA' ]
        name = 'Samples profiled for cnas and sequencing'
        description = 'This is this case list that contains all samples that are profiled for mutations and cnas.'
        stable_id = '{0}_cnaseq'.format(study)
    elif data_type == 'cna_seq_rna':
        # make a list of samples for which cna and maf and rna files are available
        samples = [i.split(',')[1] for i in content if i.split(',')[2].upper() != 'NA' and i.split(',')[3].upper() != 'NA' and i.split(',')[4].upper() != 'NA']
        name = 'Samples profiled for all of mutation, cnas, and rnaseq'
        description = 'This is this case list that contains all samples that are profiled for mutations, cnas, and rnaseq.'
        stable_id = '{0}_3way_complete'.format(study)
    
    # write outputfile if samples exist
    if samples:
        newfile = open(outputfile, 'w')
        L = ['cancer_study_identifier: {0}'.format(study),
             'stable_id: {0}'.format(stable_id),
             'case_list_name: {0}'.format(name),
             'case_list_description: {0}'.format(description),
             'case_list_ids: {0}'.format('\t'.join(samples))]
        newfile.write('\n'.join(L))
        newfile.close()
    
    
    
def get_clinical_data(clinical_info):
    '''
    (str) -> dict
    
    Returns a dictionary with clinical information for patient,sample pairs
    
    Parameters
    ----------
    - clinical info (str): Path to the file with sample clinical information
    '''
    
    # create a dict to store clinical info {'patient;sample': 'field': value}
    D = {}
    
    infile = open(clinical_info)    
    header = infile.readline().rstrip().split('\t')
    for line in infile:
        line = line.rstrip()
        if line:
            line = line.split('\t')
            patient, sample = line[0], line[1]
            ID = patient + ';' + sample
            D[ID] = {}
            for i in range(len(line)):
                if i >=2:
                    D[ID][header[i]] = line[i]
    infile.close()
    return D


def map_column_data_type(sample_info):
    '''
    (dict) -> dict
    
    Returns a dictionary with the type of the data corresponding to each clinical field in sample information
        
    Parameters
    ----------    
    - sample_info (dict): Dictionary with patient and sample information. Populates the sample clinical file
    '''
    
    data_type = {}
    for i in sample_info:
        for j in sample_info[i]:
            try:
                float(sample_info[i][j])
                data_type[j] = 'NUMBER'
            except:
                data_type[j] = 'STRING'
    return data_type


def list_column_names(sample_info):
    '''
    (dict) -> list
    
    Returns a list of fields with clinical sample information from the sample_info dictionary
   
    Parameters
    ----------    
    - sample_info (dict): Dictionary with patient and sample information. Populates the sample clinical file
    '''
    
    # make a list of column names
    c = [list(sample_info[i].keys()) for i in sample_info]
    column_names = []
    for i in c:
        column_names.extend(i)
        column_names = list(set(column_names))
    return column_names


def check_column_names(column_names):
    '''
    (list) -> None
    
    Raise an Error if the column names of the user clinical data file includes MUTATION_COUNT or FRACTION_GENOME_ALTERED 
    These two fields are auto-populated and cannot be in the data_clinical_samples.txt file
    
    Parameters
    ----------
    - column_names (list): List of column names in the clinical samples file
    '''
    
    # check if column names include banned columns
    if 'mutation_count' in list(map(lambda x: x.lower(), column_names)) or 'fraction_genome_altered' in list(map(lambda x: x.lower(), column_names)):
        raise ValueError('MUTATION_COUNT and FRACTION_GENOME_ALTERED are auto-populated clinical attributes and are banned from clinical data files')

def map_columns_to_header(column_names, header):
    '''
    (list, list) -> dict
    
    Returns a dictionary with the index of the clinical fields of the user clinical sample information file
    in the header of the data_clinical_samples.txt if present
        
    Parameters
    ----------
    - column_names (list): List of column names in the user clinical sample information file
    - header (list): Lists of column names, data_types and priority from the header of the data_clinical_samples.txt file
    '''
    
    positions = {}
    for i in column_names:
        for j in range(len(header)):
            for k in range(len(header[j])):
                if i.lower() == header[j][k].replace('#', '').lower():
                    positions[i] = k
    return positions



def update_clinical_sample_header(column_names, positions, header, data_type):
    '''
    (list, dict, list, dict)

    Returns an updated header including the clinical fields from column_names and an updated dictionary positions
    with the index of these fields in header    
        
    Parameters
    ----------
    - column_names (list): List of column names in the user clinical sample information file
    - positions (dict): Dictionary with index of the fields in column_names in the data_clinical_samples.txt header
    - header (list): Lists of column names, data_types and priority from the header of the data_clinical_samples.txt file
    - data_type (dict): Dictionary with the type of the data corresponding to each field in column names    
   '''
    
    for i in column_names:
        if i not in positions:
            # add column name to header
            for k in [0,1,4]:
                header[k].append(i.upper())
            # add data type
            header[2].append(data_type[i])
            # add priority
            header[3].append('1')
            # record postion of column
            positions[i] = len(header[0]) - 1
    return header, positions


def write_minimal_clinical_information(outputfile, mapfile, data_type, centre, sample_info = None):
    '''
    (str, str, str, str, dict | None) -> None
    
    Write clinical files with minimal clinical information
    
    Parameters
    ----------    
    - outputfile (str): Path to the outputfile
    - mapfile (str): Mapping file (map.csv) that contains paths to maf, seg, gep and mavis files    
    - data_type(str): Sample or patient
    - centre (str): Genomic centre (eg TGL, OICR)
    - sample_info (dict | None): Dictionary with patient and sample information. Populates the sample clinical file
    '''
    
    # make a list with sample names and libraries
    infile = open(mapfile)
    content = infile.read().rstrip().split('\n')
    infile.close()
    S = [i.split(',')[0:2] for i in content]
        
    if data_type == 'patient':
        # make a list of unique records
        U = []
        for i in S:
            record = [i[0], centre]
            if record not in U:
                U.append(record)
        T = ['#Patient Identifier\tCentre\tAGE DIAGNOSIS\tSEX\tETHNICITY',
             '#Patient Identifier\tCentre\tAGE DIAGNOSIS\tSEX\tETHNICITY',
             '#STRING\tSTRING\tNUMBER\tSTRING\tSTRING']
        T.append('#1' + ('\t1' * (len(T[0].split('\t')) -1)))     
        T.append('PATIENT_ID\tCENTRE\tAGE\tSEX\tETHNICITY')
    elif data_type == 'sample':
        # build header
        T = ['#Patient Identifier\tSample Identifier\tCosmic Signature\tPRIMARY SITE\tCANCER TYPE\tCLOSEST TCGA\tSAMPLE ANATOMICAL SITE\tSAMPLE PRIMARY OR METASTASIS\tTREATMENT STATUS\tPATHOLOGICAL REVIEW\tPRIOR CLINCAL TEST RESULTS\tMEAN COVERAGE\tPCT V7 ABOVE 80X\tPCT CALLABILITY\tSEQUENZA PURITY FRACTION\tSEQUENZA PLOIDY\tTMB PER MB\tHRD SCORE\tMSI STATUS',
             '#Patient Identifier\tSample Identifier\tCosmic Signature\tPRIMARY SITE\tCANCER TYPE\tCLOSEST TCGA\tSAMPLE ANATOMICAL SITE\tSAMPLE PRIMARY OR METASTASIS\tTREATMENT STATUS\tPATHOLOGICAL REVIEW\tPRIOR CLINCAL TEST RESULTS\tMEAN COVERAGE\tPCT V7 ABOVE 80X\tPCT CALLABILITY\tSEQUENZA PURITY FRACTION\tSEQUENZA PLOIDY\tTMB PER MB\tHRD SCORE\tMSI STATUS',
             '#STRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tNUMBER\tNUMBER\tNUMBER\tNUMBER\tNUMBER\tNUMBER\tNUMBER\tSTRING']
        T.append('#1' + ('\t1' * (len(T[0].split('\t')) -1)))     
        T.append('PATIENT_ID\tSAMPLE_ID\tCOSMIC_SIGS\tCANCER_TYPE\tCANCER_TYPE_DETAILED\tCLOSEST_TCGA\tSAMPLE_ANATOMICAL_SITE\tSAMPLE_PRIMARY_OR_METASTASIS\tTREATMENT_STATUS\tPATHOLOGICAL_REVIEW\tPRIOR_CLINCAL_TEST_RESULTS\tMEAN_COVERAGE\tPCT_V7_ABOVE_80X\tFRAC_CALLABILITY\tSEQUENZA_PURITY_FRACTION\tSEQUENZA_PLOIDY\tTMB_PER_MB\tHRD_SCORE\tMSI_STATUS')
             
        # check if sample info was provided by user
        if sample_info:
            # convert header to list of lists
            for i in range(len(T)):
                T[i] = T[i].split('\t')
            # make a list of column names
            column_names = list_column_names(sample_info)
            # map column name with data type
            data_types = map_column_data_type(sample_info)
            # check if column names include banned columns
            check_column_names(column_names)
            # get the index of the column names in the clinical samples file if present
            positions = map_columns_to_header(column_names, T)
            # update header with the column names and record the column names positions in the new header
            T, positions = update_clinical_sample_header(column_names, positions, T, data_types)
            # reformat header
            for i in range(len(T)):
                T[i] = '\t'.join(T[i])
            # initialize all columns with empty values beside patient and sample Ids
            data = {}
            for i in S:
                ID = i[0] + ';' + i[1]
                data[ID] = [i[0], i[1]] + ['' for j in range((len(T[0]) - 2))]
            # add values to clinical fields
            for ID in sample_info:
                 for column in sample_info[ID]:
                     data[ID][positions[column]] = sample_info[ID][column]                
            # make a list of clinical data
            U = []
            for i in data:
                U.append('\t'.join(data[i]))
        else:
            # make a list of unique records
            U = []
            for i in S:
                record = [i[0], i[1]]
                if record not in U:
                    U.append(record)
    
    # add data to table        
    if sample_info:
        for i in U:
            T.append(i)
    else:
        for i in U:
            T.append('\t'.join(i + [''] * (len(T[0].split('\t')) - 2)))
    
    newfile = open(outputfile, 'w')
    newfile.write('\n'.join(T))
    newfile.close()
    

def write_clinical_oncokb(outputfile, mapfile, cancer_code):
    '''
    (str, str) -> None
    
    Write clinical file for oncokb-annotator
    
    Parameters
    ----------    
    - outputfile (str): Path to the outputfile
    - mapfile (str): Mapping file (map.csv) that contains paths to maf, seg, gep and mavis files    
    - cancer_code (str): Cancer code from OncoTree
    '''
       
    # get library and cancer type from the mapfile
    infile = open(mapfile)    
    content = infile.read().rstrip().split('\n')
    infile.close()
    # create a list of samples
    L = [i.split(',')[1].strip() for i in content]
    if L:
        newfile = open(outputfile, 'w')
        header = ['SAMPLE_ID', 'ONCOTREE_CODE']
        newfile.write('\t'.join(header))
        for i in L:
            newfile.write('\t'.join([i, cancer_code]) + '\n')
        newfile.close()



def link_files(outdir, mapfile, data_type):
    '''
    (str, str, str) -> None
    
    Link data_type files listed in mapfile in the corresponding sub-directory of outdir 
    
    Parameters
    ----------
    - outdir (str): Path to the output directory where mafdir, sgedir, fusdir and gepdir folders are located
    - mapfile (str): Mapping file (map.csv) that contains paths to maf, seg, gep and mavis files
    - data_type (str): File type to link in their respective folders.
                       Accepted values: maf, gep, fus, and seg
    Precondition: The folders for each file type already exist
    '''
    
    # create input directories for each file type from map file    
    infile = open(mapfile)
    content = infile.read().rstrip().split('\n')
    infile.close()
    for i in range(len(content)):
        content[i] = list(map(lambda x: x.strip(), content[i].split(',')))
    
    # make a list of samples and files for which files exist
    if data_type == 'maf':
        # get the maf files
        j = 2
        extension = '.maf.gz'
    elif data_type == 'seg':
        # get the cna files
        j = 3
        extension = '.seg'
    elif data_type == 'gep':
        # get the rna files
        j = 4
        extension = '.rsem'
    elif data_type == 'fus':
        # get the fusion files
        j = 5
        extension = '.fus'
        
    samples = [i[1] for i in content if i[j].upper() != 'NA' and os.path.isfile(i[j])]
    files = [i[j] for i in content if i[j].upper() != 'NA' and os.path.isfile(i[j])]
    assert len(files) == len(samples)
    
    if files:
        for i in range(len(files)):
            folder = os.path.join(outdir, '{0}dir'.format(data_type))
            os.makedirs(folder, exist_ok=True)
            target = os.path.join(folder,  samples[i] + extension)
            subprocess.call('ln -s {0} {1}'.format(files[i], target), shell=True)         
    else:
        print('Cannot link {0} files. No files exist in mapping file {1}'.format(data_type, mapfile))


def link_rsem_files_comparison(compdir, gepcomp):
    '''
    (str, str) -> None   
    
    Links the rsem files listed in gepcomp in the compdir directory
    
    Parameters
    ----------
    - compdir (str): 
    - gepcomp (str): Path to the gepcomp file with comma-separated list of samples, rsem paths    
    '''
 
    infile = open(gepcomp)
    samples, files = [], []
    for line in infile:
        line = line.rstrip()
        if line:
            line = line.split(',')
            # check that file is valid
            if os.path.isfile(line[1].strip()):
                # collect sample and file
                samples.append(line[0].strip())
                files.append(line[1].strip())
    infile.close()
    assert len(samples) == len(files)            
            
    if files:
        for i in range(len(files)):
            target = os.path.join(compdir, samples[i] + '.rsem')
            if os.path.exist(target) == False:
                subprocess.call('ln -s {0} {1}'.format(files[i], target), shell=True)    
            

def get_sample_from_filename(file):
    '''
    (str) -> str
    
    Returns the sample name from the file name
    Precondition: The file is named after sample and extension

    Parameters
    ----------
    - file (str): Path to the file
    '''

    # get the base name of the file
    filename = os.path.basename(file)
    # get file name without extension
    name = filename[:filename.index('.')]
    return name                
            
    
def concatenate_seg_files(segdir, outputfile):
    '''
    (str, str) -> None
    
    Concatenate seg files located in segdir into outputfile
    
    Parameters
    ----------
    - segdir (str): Directory with seg files
    - outputfile (str): Path to the concatenated seg file
    '''

    # make a list of seg files
    segfiles = [os.path.join(segdir, i) for i in os.listdir(segdir) if '.seg' in i]
    # get the header of the seg file
    infile = open(segfiles[0])
    header = infile.readline()
    infile.close()
    
    newfile = open(outputfile, 'w')
    newfile.write(header)
    
    for file in segfiles:
        # extract the sample name from file name
        sample = get_sample_from_filename(file)
        # get content of seg file
        infile = open(file)
        content = infile.read().rstrip().split('\n')
        infile.close()
        # remove header
        content.pop(0)
        # replace ID field with sample name
        for i in range(len(content)):
            content[i] = content[i].split('\t')
            content[i][0] = sample
            content[i] = '\t'.join(content[i])
        # write content of seg file to concatenated file
        newfile.write('\n'.join(content) + '\n')
    newfile.close()        


def write_metadata(outputfile, study, data_type, genome):
    '''
    (str, str, str, str) -> None

    Write CNA metadata for a given data type
        
    Parameters
    ----------
    - outputfile (str): Path to the output file
    - study (str): Name of the study as it appears in cBioPortal
    - data_type (str): Type of CNA metadata.
                       Accepted valued: discrete, log2-value, seg, expression, fusion, zscore, maf
    - genome (str): Reference genome (hg19 or hg38)
    '''
    
    if data_type == 'discrete':
        stable_id = 'gistic'
        description = 'profile_description: Putative copy-number calls:  Values: -2=homozygous deletion; -1=hemizygous deletion; 0=neutral/no change; 1=gain; 2=high level amplification'
        name = 'Putative copy-number alterations from GISTIC'
        filename = 'data_CNA.txt'
        alteration = 'COPY_NUMBER_ALTERATION'
        data = data_type.upper()
    elif data_type == 'log2-value':
        stable_id = 'log2CNA'
        description = 'profile_description: Log2 copy-number values'
        name = 'Log2 copy-number values' 
        filename = 'data_log2CNA.txt'
        alteration = 'COPY_NUMBER_ALTERATION'
        data = data_type.upper()
    elif data_type == 'seg':
        filename = 'data_segments.txt'
        description = 'description: Segment data'
        alteration = 'COPY_NUMBER_ALTERATION'
        data = data_type.upper()
    elif data_type == 'fusion':
        alteration = 'FUSION'
        stable_id = 'fusion' 
        name = 'Fusions'
        filename = 'data_fusions.txt'
        description = 'profile_description: Fusion data'
        data = data_type.upper()
    elif data_type == 'expression':
        alteration = 'MRNA_EXPRESSION'
        data = 'CONTINUOUS'
        stable_id = 'rna_seq_mrna'
        filename = 'data_expression.txt'
        name = 'mRNA expression RNA-Seq'
        description = 'profile_description: Expression levels RNA-Seq'
    elif data_type == 'zscore':
        alteration = 'MRNA_EXPRESSION'
        data = 'Z-SCORE'
        stable_id = 'rna_seq_mrna_median_Zscores'
        filename = 'data_expression_zscores.txt'
        description = 'profile_description: Expression levels z-scores'
        name = 'mRNA expression z-scores'
    elif data_type == 'maf':
        data = data_type.upper()
        alteration = 'MUTATION_EXTENDED'
        description = 'profile_description: Mutations'
        name = 'profile_name: Mutations'
        stable_id = 'mutations'
        filename = 'data_mutations_extended.txt'

    # collect file text commun to all metadata data types
    L = ['cancer_study_identifier: {0}'.format(study),
         'data_filename: {0}'.format(filename),
         'datatype: {0}'.format(data),
         description,
         'genetic_alteration_type: {0}'.format(alteration)]

    # add specific data type text
    if data_type in ['discrete', 'log2-value', 'fusion', 'expression', 'zscore', 'maf']:
        L.extend(['stable_id: {0}'.format(stable_id),
                  'show_profile_in_analysis_tab: true',
                  'profile_name: {0}'.format(name)])
    elif data_type == 'seg':
        L.append('reference_genome_id: {0}'.format(genome))

    newfile = open(outputfile, 'w')
    newfile.write('\n'.join(L))
    newfile.close()        
    


def concatenate_fusion_files(fusdir, outputfile):
    '''
    (str, str) -> None
    
    Concatenate fusion files located in fusdir into outputfile
    
    Parameters
    ----------
    - fusdir (str): Directory with fusion files
    - outputfile (str): Path to the concatenated seg file
    '''

    # make a list of fusion files
    fusfiles = [os.path.join(fusdir, i) for i in os.listdir(fusdir) if '.fus' in i]
    # get the header of the seg file
    infile = open(fusfiles[0])
    header = infile.readline().rstrip().split('\t')
    infile.close()
    
    # add Sample to header
    header.insert(0, 'Sample')
    header = '\t'.join(header)
    # write header to outputfile
    newfile = open(outputfile, 'w')
    newfile.write(header)
    
    for file in fusfiles:
        # extract sample name from file name
        sample = get_sample_from_filename(file)
        # get content of seg file
        infile = open(file)
        content = infile.read().rstrip().split('\n')
        infile.close()
        # remove header
        content.pop(0)
        # add sample to content
        for i in range(len(content)):
            content[i] = content[i].split('\t')
            content[i].insert(0, sample)
            content[i] = '\t'.join(content[i])
        # write content of seg file to concatenated file
        newfile.write('\n'.join(content) + '\n')
    newfile.close()        

    

def list_gep_samples(gepdir, outputfile):
    '''
    (str, str) -> None
    
    Write list of samples with rna data to outputfile
    
    Parameters
    ----------
    - gepdir (str): Directory with rsem files
    - outputfile (str): path to the outputfile
    '''

    # make a list of rsem files
    gepfiles = [os.path.join(gepdir, i) for i in os.listdir(gepdir) if '.rsem' in i]
    # make a list of sample names from the gep file names
    samples = [get_sample_from_filename(file) for file in gepfiles]
    # write samples to file
    newfile = open(outputfile, 'w')
    newfile.write('\n'.join(samples) + '\n')
    newfile.close()     
    


def extract_fpkm(gepfile):
    '''
    (str) -> dict
    
    Returns a dictionary of gene, fpkm key, value pairs
    
    Parameters
    ----------
    - gepfile (str): Path to the rsem file with fkpm counts
    '''
    
    # create a dict to store fpkm for each gene
    D = {}
    infile = open(gepfile)
    # skip header
    infile.readline()
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split('\t')
            gene = line[0]
            fpkm = line[6]
            assert gene not in D
            D[gene] = fpkm
    infile.close()
    return D
    

def collect_fpkm(gepdir):
    '''
    (str) -> dict
    
    Returns a dictionary with fpkm values for all genes and samples with rna data in directory gepdir
    
    Parameters
    ----------
    - gepdir (str): Directory with rsem files
    '''
    
    # make a list of rsem files
    gepfiles = [os.path.join(gepdir, i) for i in os.listdir(gepdir) if '.rsem' in i]
    # create a dict to store fpkm for each gene and sample
    D = {}
    for file in gepfiles:
        # get sample name from file name
        sample = get_sample_from_filename(file)
        # get the fkpm for each gene
        fpkm = extract_fpkm(file)
        assert sample not in D
        D[sample] = fpkm
    return D


def write_fpkm_to_file(D, outputfile):
    '''
    (dict, str) -> None
    
    Write the fpkm values for all samples and genes in dictionary D to outputfile
    
    Parameters
    ----------
    - D (dict): Dictionary with fpkm for all samples and gene {sample: {gene: fpkm}}
    - outputfile (str): Path to the outputfile
    '''
    
    # list all samples
    samples = list(D.keys())
    # list all genes
    genes = []
    for i in samples:
        genes.extend(list(D[i].keys()))
        genes = list(set(genes))
        
    # write fpkm to file
    newfile = open(outputfile, 'w')    
    header = ['gene_id'] + samples
    newfile.write('\t'.join(header) + '\n')    
    for gene in genes:
        line = [gene]
        for sample in samples:
            line.append(D[sample][gene])
        newfile.write('\t'.join(line) + '\n')    
    newfile.close()
    

def concatenate_fpkm_from_gep_files(gepdir, outputfile):
    '''
    (str, str) -> None
    
    Write fpkm for all genes and samples with rna data to outputfile
    
    Parameters
    ----------
    - gepdir (str): Directory with rsem files
    - outputfile (str): Path to the outputfile with fpkm for each sample and gene
    '''
    
    # collect all fpkm for each gene and sample 
    D = collect_fpkm(gepdir)
    # write fpkm to outputfile
    write_fpkm_to_file(D, outputfile)


def update_fpkm_file(fpkm_file, compdir):
    '''
    (str, str) -> None
    
    Update fpkm file with expression data of samples in compdir
    
    Parameters
    ----------
    - fpkm_file (str): Path to the file with fpkm for each gene and sample
    - compdir (str): Directory with rsem files located in gepcomp is option is used in config
    '''

    # collect fpkm from rsem files in compdir
    D = collect_fpkm(compdir)
    
    # extract fpkm from fpkm_file
    infile = open(fpkm_file)
    # make a list of samples as theu appear in header
    samples = infile.readline().rstrip().split('\t')[1:]    
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split('\t')
            gene = line[0]
            fpkm = line[1:]
            # updare D with sample if not already recorded
            for i in range(len(samples)):
                if samples[i] not in D:
                    D[samples[i]] = {}
                # keep fpkm from sample in current study if already recoded in gepcomp
                D[samples[i]][gene] = fpkm[i]

    # update fpkm_file with data from gepcomp 
    write_fpkm_to_file(D, fpkm_file)
    
    
def get_maf_header(maffile):
    '''
    (str) -> str
    
    Returns the header of the maf file
    
    Parameters
    ----------
    - maffile (str): Path to maf file
    '''
    
    header = ''
    infile = gzip.open(maffile, 'rt')
    for line in infile:
        # #version may or not be 1st line. loop until header is found 
        if 'Hugo_Symbol' in line:
            header = line.rstrip()
    infile.close()
    return header
    


def concatenate_maf_files(mafdir, outputfile):
    '''
    (str, str) -> None
    
    Concatenates all the gzipped maf files in mafdir into a text file outputfile
    with column Tumor_Sample_Barcode replaced by the sample name 
    
    Parameters
    ----------
    - mafdir (str): Directory containing the maf files
    - outputfile (str): Path to the concatenated maf file (unzipped)
    '''
     
    # make a list of maf files
    maffiles = [os.path.join(mafdir, i) for i in os.listdir(mafdir) if '.maf.gz' in i]
    # get the header of the maf file
    header = get_maf_header(maffiles[0])
    
    newfile = open(outputfile, 'w')
    newfile.write(header + '\n')
    for file in maffiles:
        # get the sample name from the file name
        sample = get_sample_from_filename(file)
        # read the content of the file
        infile = gzip.open(file, 'rt')
        # skip header
        for line in infile:
            if 'version' not in line and 'Hugo_Symbol' not in line:
                line = line.rstrip()
                if line != '':
                    line = line.split('\t')
                    # replace column 'Tumor_Sample_Barcode' with sample name
                    line[15] = sample
                    # write to temp file
                    newfile.write('\t'.join(line) + '\n')
    newfile.close()

                
def filter_mutations(maffile, outputfile, depth_filter, alt_freq_filter, gnomAD_AF_filter):
    '''
    (str, str, int, float, float) -> (int, int)
    
    Writes records from maffile to outputfile if mutations pass depth, alt_freq and gnomAd_AF filters
    Returns the total number of mutations before and after filtering
    Precondition: the maf file is unzipped.
    
    Parameters
    ----------
    - maffile (str): Path to the maf file (unzipped)
    - outputfile (str): Path to the output file
    - depth_filter (int): Minimum number of reads at a given position
    - alt_freq_filter (float): Minimum alternative allele frequency (t_alt_count / t_depth)
    - gnomAD_AF_filter (float): Maximum allele frequency is the Genome Aggregation Database
                                if Matched_Norm_Sample_Barcode is unmatched
    '''
    
    # count the number of mutations before and after filtering
    total, kept = 0, 0
    
    # open files
    newfile = open(outputfile, 'w')
    infile = open(maffile)
    
    # get file header
    header = infile.readline().rstrip().split('\t')
    
    # write header to outputfile
    newfile.write('\t'.join(header) + '\n')
    
    # make a list of accepted mutations
    valid_mutations = ['Frame_Shift_Del',  'Frame_Shift_Ins', 'In_Frame_Del',
                          'In_Frame_Ins', 'Missense_Mutation', 'Nonsense_Mutation',
                          'Nonstop_Mutation', 'Silent', 'Splice_Site', 'Translation_Start_Site']
    exclude = ['str_contraction', 't_lod_fstar']    

    # apply maf filters to all mutations
    for line in infile:
        # count total mutations
        total += 1
        line = line.rstrip()
        newline = ''
        if line != '':
            # check that mutation is valid and that excuded fields are not recorded
            mutations = [i in line for i in valid_mutations]     
            non_valid = [i in line for i in exclude] 
            if any(mutations) and not any(non_valid):
                # apply filters to mutations
                line = line.split('\t')
                # filter based on depth
                if int(line[header.index('t_depth')]) >= depth_filter:
                    # filter based on ratio t_alt_count / t_depth
                    if int(line[header.index('t_alt_count')]) / int(line[header.index('t_depth')]) >= alt_freq_filter:
                        # filter based on Matched_Norm_Sample_Barcode
                        if line[header.index('Matched_Norm_Sample_Barcode')] == "unmatched":
                            # check gnomAD_AF. field may be blank, check if value recorded
                            try:
                                float(line[header.index('gnomAD_AF')])
                            except:
                                # no value for gnomAD_AF, do not keep mutation
                                newline = ''
                            else:
                                # compare gnomAD_AF to folder
                                if float(line[header.index('gnomAD_AF')]) < gnomAD_AF_filter:
                                    newline = line
                                    kept += 1
                        else:
                            newline = line
                            kept += 1
                            
            if newline:
                newfile.write('\t'.join(newline) + '\n')
    newfile.close()                
    return total, kept                
  

def remove_indels(maffile, outputfile):
    '''
    (str, str) -> (int, int)
    
    Writes records from maffile to outputfile with indels removed.
    Returns the total number of mutations before and after filtering
    Precondition: the maf file is unzipped.
    
    Parameters
    ----------
    - maffile (str): Path to the maf file (unzipped)
    - outputfile (str): Path to the output file
    '''
    
    # count the number of mutations before and after filtering
    total, kept = 0, 0
    
    # open files
    newfile = open(outputfile, 'w')
    infile = open(maffile)
    
    # get file header
    header = infile.readline().rstrip().split('\t')
    
    # write header to outputfile
    newfile.write('\t'.join(header) + '\n')
    
    for line in infile:
        # count total mutations
        total += 1
        line = line.rstrip()
        if line != '':
            if line[header.index('Variant_Type')] not in ['DEL', 'INS']:
                # record mutations without indels and update counter
                newfile.write('\t'.join(line) + '\n')
                kept += 1
    newfile.close()                
    return total, kept                


def process_cna(segfile, genebed, oncolist, gain, amp, htz, hmz, ProcCNA, outdir, genelist):
    '''
    (str, str, str, int, int, int, int, str, str, str | None) -> None
    
    Process segmentation data to R script ProcCNA.r to generate data files data_segments.txt,
    data_log2CNA.txt, data_CNA.txt and data_CNA_short.txt

    Parameters
    ----------
    - segfile (str): Path to the concatenated segmentation file
    - genebed (str): Path to Tab-delimited 5 column bed file which defines the genomic positions of the canonical genes
    - oncolist (str): Path to list of cancer genes
    - gain (float): Threshold for CNA gain     
    - amp (float): Threshold for CNA amplification        
    - htz (float): Threshold for heterozygous deletion
    - hmz (float): Threshold for homozygoys deletion
    - ProcCNA (str): Path to the R script for CNA processing
    - outdir (str): - outdir (str): Path to the output directory where mafdir, sgedir, fusdir and gepdir folders are located
    - genelist (str | None): Path to list of Hugo Symbols. (Optional)
    '''    

    cmd = 'Rscript {0} {1} {2} {3} {4} {5} {6} {7} {8} {9}'.format(ProcCNA, segfile, genebed, genelist, oncolist, gain, amp, htz, hmz, outdir)
    print(cmd)
    
    if os.path.isfile(ProcCNA):
        exit_code = subprocess.call(cmd, shell=True)
        if exit_code:
            sys.exit('Could not process CNAs.')
    else:
        raise FileNotFoundError('Cannot find R script path {}'.format(ProcCNA))
       
        
def process_rna(gepfile, enscon, genelist, ProcRNA, outdir):
    '''
    (str, str, str | None, str, str) -> None
    
    Process RNAseq expression data through R script ProcRNA to generate data_expression.txt
    and data_expression_zscores.txt files
 
    Parameters
    ----------
    - gepfile (str): Path to concatenated file with expression data
    - enscon (str): path to tab-delimited 2 column file of ENSEMBLE gene ID and Hugo_Symbol 
    - genelist (str | None): Path to list of Hugo Symbols (optional)
    - ProcRNA (str) Path to R script ProcRNA.r
    - outdir (str): - outdir (str): Path to the output directory where mafdir, sgedir, fusdir and gepdir folders are located
    '''
    
    cmd = 'Rscript {0} {1} {2} {3} {4}'.format(ProcRNA, gepfile, enscon, genelist, outdir)
    print(cmd)
    
    if os.path.isfile(ProcRNA):
        exit_code = subprocess.call(cmd, shell=True)
        if exit_code:
            sys.exit('Could not process RNAseq expression.')
    else:
        raise FileNotFoundError('Cannot find R script path {}'.format(ProcRNA))
    

def process_fusion(fusfile, entcon, min_fusion_reads, ProcFusion, outdir):
    '''
    (str, str, int, str, str) -> None    

    Process fusion data through R script ProcFusion.r to generate fusion data data_fusions.txt
    
    Parameters
    ----------
    - fusfile (str): Path to concatenated fusion file
    - entcon (str): Path to tab-delimited 2 column file of ENTREZ gene ID and Hugo_Symbol 
    - min_fusion_reads (int): mininimum number of reads for fusion calls
    - ProcFusion (str) Path to the R script ProcFusion.r 
    - outdir (str): - outdir (str): Path to the output directory where mafdir, sgedir, fusdir and gepdir folders are located
    ''' 
    
    cmd = 'Rscript {0} {1} {2} {3} {4}'.format(ProcFusion, fusfile, entcon, min_fusion_reads, outdir)
    print(cmd)
    
    if os.path.isfile(ProcFusion):
        exit_code = subprocess.call(cmd, shell=True)
        if exit_code:
            sys.exit('Could not process fusions.')
    else:
        raise FileNotFoundError('Cannot find R script path {}'.format(ProcFusion))


def process_mutations(maffile, tglpipe, ProcMAF, outdir):
    '''
    (str, bool, bool, str, str) -> None
    
    Process variant calls through R script ProcMaf.r to generate mutaion data files
    data_mutations_extended.txt, unfiltered_data_mutations_extended.txt and weights.txt
    
    Parameters
    ----------
    - maffile (str): Path to concatenated maf file with variant calls
    - tglpipe (bool):  filter variants according to TGL specifications if True
    - ProcMaf (str): Path to R script ProcMaf.r
    - outdir (str): - outdir (str): Path to the output directory where mafdir, sgedir, fusdir and gepdir folders are located
    '''
    
    tglpipe = 'TRUE' if tglpipe else 'FALSE'
        
    cmd = 'Rscript {0} {1} {2} {3}'.format(ProcMAF, maffile, tglpipe, outdir)
    print(cmd) 
    
    if os.path.isfile(ProcMAF):
        exit_code = subprocess.call(cmd, shell=True)
        if exit_code:
            sys.exit('Could not process mutations.')
    else:
        raise FileNotFoundError('Cannot find R script path {}'.format(ProcMAF))
    

def get_sample_info(clinical_samples, outputfile):
    '''
    (str, str) -> None
    
    Write clinical data for samples compatible with IGV
    
    Parameters
    ----------
    - clinical_samples (str): Path to the file with clinical samples information
    - outputfile (str): Path to the outputfile with content re-formatted for IGV
    '''
    
    infile = open(clinical_samples)
    # find the header, skip commented lines
    samples = []
    for line in infile:
        if not line.startswith('#'):
            if line.startswith('PATIENT_ID'):
                header = line
            else:
                samples.append(line.rstrip())
    infile.close()
    while '' in samples:
        samples.remove('')            
    for i in range(len(samples)):
        samples[i] = samples[i].split('\t')
    # edit header
    header = header.rstrip().split('\t')
    header[0] = 'TRACK_ID'
    header[1] = 'PARTICIPANT_ID'
        
    newfile = open(outputfile, 'w')
    newfile.write('\t'.join(header) + '\n')
    for i in samples:
        newfile.write('\t'.join([i[1], i[0]]) + '\n')
    newfile.close()
    

def create_output_directories(outdir):
    '''
    (str) -> (str, str, str)
    
    Creates outdir directory and sub-folder structure and returns the path to the
    sub-folders within outdir. Removes outdir if it already exists
    
    Parameters
    ----------
    - outdir (str): Path to output directory 
    '''
    
    # remove old output directory if it exists
    if os.path.isdir(outdir):
        print('{0} exists already; removing'.format(outdir))
        shutil.rmtree(outdir)
    # create output directory
    print('creating output directory {0}'.format(outdir))
    os.makedirs(outdir, exist_ok=True)

    # create output folders
    cbiodir = os.path.join(outdir, 'cbioportal_import_data')
    casedir = os.path.join(cbiodir, 'case_lists')
    suppdir = os.path.join(outdir, 'supplementary_data')
    print('creating output sub-folders:', cbiodir, casedir, suppdir, sep = '\n')
    for i in [cbiodir, casedir, suppdir]:
        os.makedirs(i, exist_ok=True)
    return cbiodir, casedir, suppdir


def check_configuration(config):
    '''
    (configparser.ConfigParser) -> None
    
    Check the content of the config file and raise a ValueError if required information is missing or not as expected
    
    Parameters
    ----------
    - config (configparser.ConfigParser): Config file parsed with configparser
    '''
    
    # raise an error if a section is omitted
    missing_sections = [i for i in ['Resources', 'Options', 'Parameters', 'Filters'] if i not in config.sections()] 
    if missing_sections:
        raise ValueError('ERROR. Missing sections {0} from config'.format(', '.join(missing_sections)))
        
    # check paths from resources
    expected_resources = ['procmaf', 'proccna', 'procrna', 'procfusion', 'token', 'enscon_hg38', 'enscon_hg19', 'entcon', 'genebed_hg38', 'genebed_hg19', 'genelist', 'oncolist']
    missing_resources = [i for i in expected_resources if i not in list(config['Resources'].keys())]
    if missing_resources:
        raise ValueError('ERROR. Missing resources: {0}'.format(', '.join(missing_resources)))
    invalid_resource_files = [i for i in ['procmaf', 'proccna', 'procrna', 'procfusion', 'token', 'enscon_hg38', 'enscon_hg19', 'entcon', 'genebed_hg38', 'genebed_hg19', 'genelist', 'oncolist'] if config['Resources'][i] and os.path.isfile(config['Resources'][i]) == False]
    if invalid_resource_files:
        raise ValueError('ERROR. Provide valid path for {0}'.format(', '.join(invalid_resource_files)))
    
    # check options
    missing_options = [i for i in ['mapfile', 'outdir', 'study', 'center', 'cancer_code'] if config['Options'][i] is None]
    if missing_options:
        raise ValueError('ERROR. Provide values for {0} in the config'.format(', '.join(missing_options)))
    # check map file
    if os.path.isfile(config['Options']['mapfile']) == False:
        raise ValueError('ERROR: Provide valid path to mapfile in config')
    
    # check parameters
    expected_parameters = ['gain', 'amplification', 'heterozygous_deletion', 'homozygous_deletion', 'minfusionreads']
    missing_parameters = [i for i in expected_parameters if i not in list(config['Parameters'].keys())]
    if missing_parameters:
        raise ValueError('ERROR. Missing Parameters: {0}'.format(', '.join(missing_parameters)))
    # check value types
    for i in ['gain', 'amplification', 'heterozygous_deletion', 'homozygous_deletion', 'minfusionreads']:
        try:
            config['Parameters'].getfloat(i)
        except:
            raise ValueError('ERROR. {0} is not a number'.format(i))
            
    # check filters
    expected_filters = ['tglpipe', 'filter_variants', 'depth_filter', 'alt_freq_filter', 'gnomAD_AF_filter', 'filter_indels']
    missing_filters = [i for i in expected_filters if i not in list(config['Filters'].keys())]
    if missing_parameters:
        raise ValueError('ERROR. Missing Filters: {0}'.format(', '.join(missing_filters)))
    # check that filter parameters are provided if filtering of variants is expected
    try:
        filter_variants = config['Filters'].getboolean('filter_variants')
    except:
        raise ValueError('ERROR. {0} is not a boolean'.format('filter_variants'))
    finally:
        if filter_variants:
            missing_variant_filters = [i for i in ['depth_filter', 'alt_freq_filter', 'gnomAD_AF_filter'] if not config['Filters'][i]]
            if missing_variant_filters:
                raise ValueError('ERROR. Expecting variant filtering. Provide values for {0}'.format(', '.join(missing_variant_filters)))
    # check value types
    for i in ['tglpipe', 'filter_variants', 'filter_indels']:
        try:
            config['Filters'].getboolean(i)
        except:
            raise ValueError('ERROR. {0} is not a boolean'.format(i))
    for i in ['depth_filter', 'alt_freq_filter', 'gnomAD_AF_filter']:
        try:
            config['Filters'].getfloat(i)
        except:
            raise ValueError('ERROR. {0} is not a number'.format(i))
            


def extract_resources_from_config(config):
    '''
    (configparser.ConfigParser) -> (str, str, str, str, str, str, str, str)
    
    Returns the variables listed in the Resources section of the config
    
    Parameters
    ----------
    - config (configparser.ConfigParser): Config file parsed with configparser
    '''
    
    resources = ['procmaf', 'proccna', 'procrna', 'procfusion', 'token', 'enscon_hg38', 'enscon_hg19', 'entcon', 'genebed_hg38', 'genebed_hg19', 'genelist', 'oncolist']
    L = [config['Resources'][i] for i in resources]
    ProcMAF, ProcCNA, ProcRNA, ProcFusion, token, enscon_hg38, enscon_hg19, entcon, genebed_hg38, genebed_hg19, genelist, oncolist = L
    return ProcMAF, ProcCNA, ProcRNA, ProcFusion, token, enscon_hg38, enscon_hg19, entcon, genebed_hg38, genebed_hg19, genelist, oncolist
    

def extract_options_from_config(config):
    '''
    (configparser.ConfigParser) -> (str, str, str)
    
    Returns the variables listed in the Options section of the config
    
    Parameters
    ----------
    - config (configparser.ConfigParser): Config file parsed with configparser
    '''    
    
    options = ['mapfile', 'outdir', 'study', 'center', 'cancer_code']
    L = [config['Options'][i] for i in options]
    mapfile, outdir, study, center, cancer_code = L
    return mapfile, outdir, study, center, cancer_code


def extract_parameters_from_config(config):
    '''
    (configparser.ConfigParser) -> (int, int, int, int, int, bool, str)
    
    Returns the variables listed in the Parameters section of the config
    
    Parameters
    ----------
    - config (configparser.ConfigParser): Config file parsed with configparser
    '''
    
    # get the numeric variables
    nums = ['gain', 'amplification', 'heterozygous_deletion', 'homozygous_deletion', 'minfusionreads']
    L = [config['Parameters'].getint(i) if '.' not in config['Parameters'][i] else config['Parameters'].getfloat(i) for i in nums]
    gain, amplification, heterozygous_deletion, homozygous_deletion, minfusionreads = L
    return gain, amplification, heterozygous_deletion, homozygous_deletion, minfusionreads


def extract_filters_from_config(config):
    '''
    (configparser.ConfigParser) -> (float | None, float | None, float | None, bool, bool, bool)
    
    Returns the variables listed in the Filters section of the config
    
    Parameters
    ----------
    - config (configparser.ConfigParser): Config file parsed with configparser
    '''
    
    # get the numeric variables
    L = [config['Filters'].getfloat(i) for i in ['depth_filter', 'alt_freq_filter', 'gnomAD_AF_filter']]
    # add booleans
    L.extend([config['Filters'].getboolean(i) for i in ['tglpipe', 'filter_variants', 'filter_indels']])
    depth_filter, alt_freq_filter, gnomAD_AF_filter, tglpipe, filter_variants, filter_indels = L
    return depth_filter, alt_freq_filter, gnomAD_AF_filter, tglpipe, filter_variants, filter_indels


def check_cancer_type(cancer_code):
    '''
    (str, str) -> None
    
    Raises ValueError if cancer code is not in listed in oncoTree
    
    Parameters
    ----------
    - cancer_code (str): Cancer type as defined in http://oncotree.mskcc.org
    '''
    
    # list all the cancer types obtained with the OncoTree api
    cancers = 'MMB,GCB,SBLU,OHNCA,PAOS,TMDS,ARMS,SCST,ITLPDGI,MBC,AWDNET,AMLCBFBMYH11,\
    ROCY,LAM,SBL,AWM,AMLMRC,CABC,PCGP,SCCO,LGT,GBM,BL,MFS,AN,SELT,AMOL,PERL,AMLCEBPA,\
    UASC,LNM,PXA,PMBL,BIALCL,SRAP,LXSC,ADRENAL_GLAND,WM,RHM,WDLS,PEMESO,CNL,HNSC,BLLRGA,\
    DESM,OEC,SCLC,EMPSGC,ACBC,OPHSC,CERMS,GCEMU,ENCG,PTLD,SUBE,SCGBC,FIOS,HGSOC,THHC,LUNE,\
    MGUSIGG,UPDC,OGBL,CCPRC,PSC,RSCC,BCCA,SRCBC,PEL,CCS,UPA,MFH,FIBS,COM,BMGCT,PVMF,SCHW,\
    BTBEOV,CERVIX,SIC,NSCHL,NSCLC,OGCT,CEVG,MAC,VDYS,PHCH,BPLL,GMN,GCTSTM,BCAC,OM,BPSCC,BLL,\
    PCNSMT,MLNFGFR1,OSMBT,MDSID5Q,PHPTLD,ASPS,NETNOS,IOPN,AMPULLA_OF_VATER,MPALBCRABL1,MCBCL,\
    ALKLBCL,EATL,BLADDER,TCCA,MMBL,SWDNET,MUCC,USTUMP,USARC,OCS,CMPT,DSRCT,ALUCA,SPTCL,SPC,EBOV,\
    THPD,SGO,MBGN,MRC,CSCHW,EVN,USC,PERITONEUM,NECNOS,MIDDA,WDTC,AHCD,CSCC,BYST,OSGCT,BLLHYPER,\
    SBOV,EMBCA,ADPA,LGSOC,HGESS,EBVDLBCLNOS,SEBA,BLLIAMP21,HNNE,MCHS,NPC,APMF,LGESS,PLSMESO,\
    HCL-V,PPCT,OMGCT,OPE,TPLL,TAC,UA,CEAD,ALCLALKN,HEMA,CMCD,PLLS,DDLS,LUMEC,MDLC,ISMCL,OVT,\
    NSGCT,OYST,SPN,MDSMPNU,UMNC,IPN,MNM,SDCA,AGA,PTFL,BRSRCC,TYST,PAAD,MGCT,CECC,LMS,AMLBCRABL1,\
    MPN,SGAD,SCBC,DIG,DMBL,MRT,IMTB,MELC,BLLETV6RUNX1,LDCHL,APXA,CCHM,PINC,AMLRGA,MT,MCSL,BMGT,\
    SPB,CCM,SM,LECLC,BCLU,BLLHYPO,CUPNOS,LIMNET,IFS,RMS,OCNOS,THPA,CHRCC,PNS,ANM,PLBL,LAIS,UEC,\
    SFTCNS,ENKL,VMA,ASTR,CMML,VA,GNBL,ICPN,MDSMD,MIDDO,MOV,TNET,GNC,APE,CCHDM,BTBOV,GSARC,PTCL,\
    CTAAP,BGCT,RBL,SPCC,ANGL,AMLRARA,NST,MTNN,CCRCC,UUC,AMLMLLT3KMT2A,EGCT,BMT,UCCA,SCUP,TAM,\
    LATL,AFH,IAMPCA,MLADS,USTAD,CESE,LYMPH,LGGNOS,APLPMLRARA,MLYM,SRCC,BCCP,BRAME,AITL,PAAC,\
    PINT,MDSRS,HS,CHLPTLD,SRCCR,OIMT,VGCE,MPALKMT2A,UELMS,PRSC,HCL,CESC,SYNS,ERMS,EYE,ASM,UM,\
    CHOL,ESST,MDS/MPN,FL,SGTTL,PTH,BOWEL,WT,TET,IUP,TLL,SS,HGNEC,AODG,PORO,HTAT,LDD,HGSFT,\
    NMCHN,UMLMS,HEAD_NECK,UDMN,ACRM,MLNPDGFRA,PHC,ATLL,THRLBCL,PTAD,ESMM,MNG,LIPO,LBLIRF4,\
    DCIS,TT,PCSMTPLD,AMML,ESCA,MCC,URCC,MCD,URCA,OMT,LYG,SPDAC,NUTCL,SEF,AECA,ETMF,GINETES,\
    OAST,DES,PAASC,SNSC,ALCLALKP,MDSRSSLD,BRCA,IHCH,ODGC,EBVMCU,LCLC,STAS,COAD,SCT,EHAE,EMBT,\
    READ,CNC,UCS,BLLBCRABL1,HDCN,MTSCC,MP,MDS,VPE,FHPTLD,URMM,BLL11Q,SNUC,EGC,MDSEB,MYXO,DF,\
    CCLC,NBL,RDD,DDCHDM,IDCT,LUAD,SLCT,INTS,UAD,PROSTATE,FTCL,HCC,ACN,MDSEB2,PLMESO,NCCRCC,\
    LBGN,COADREAD,VGCT,MPE,CLNC,PGNG,MYELOID,EMALT,BONE,AMLRBM15MKL1,IPMN,MGUS,BILIARY_TRACT,\
    EMPD,CHOS,RLCLC,OSACA,PANET,GINET,SFT,EPM,HMBL,HGSOS,TLYM,PNET,TMESO,PAST,DFSP,THAP,SKLMM,\
    TGCT,PTCY,ASCT,BNNOS,MACR,IBC,GTD,TMN,MNET,MEITL,PRAD,MRTL,NKCLL,OSOS,AIS,AASTR,VIMT,AMLDEKNUP214,\
    MGUSIGA,KIDNEY,MDEP,SSM,CMML2,HNSCUP,IMS,SMN,PRSCC,DSTAD,SECOS,NMZL,PT,MSCHW,ACPG,SCCE,NPTLTFH,VMM,\
    MGST,BTMOV,CHDM,ACPP,PSTAD,MEL,CM,HGNEE,MAMPCA,PRCC,BLPT,LUPC,VMGCT,GCLC,OFMT,BPDCN,PACT,EHCH,LUSC,\
    FLC,TAML,CMC,APTAD,ATM,MYCF,MASCC,CHOM,LUNG,ALT,VPDC,AFX,MBOV,PMA,PSCC,CEAIS,ODYS,FHRCC,PD,FRCT,AUL,\
    TLGL,ARMM,DLBCLNOS,TISSUE,PCLPD,PMFOFS,PCNSM,THME,IDC,TSTAD,NSCLCPD,UNEC,CHL,OSMAD,BLCA,CSNOS,PPTLD,\
    USCC,CELI,EMCHS,HHV8DLBCL,ONBL,HGONEC,SPZM,IMTL,BLLIL3IGH,PBS,UMC,VMT,SMMCL,MPC,ASTB,ICEMU,DDCHS,\
    SACA,PADA,PLRMS,PANCREAS,HNMUCM,AML,RWDNET,CCOC,MBN,HGBCLMYCBCL2,LNET,CLLSLL,LPL,LUACC,MAAP,AMKL,\
    ABC,PMHE,PRNET,STOMACH,BA,SKIN,HGNET,CELNOS,ECAD,SOC,BLSC,MBL,CHBL,ES,UCP,GBC,JXG,DIA,PBL,OVARY,\
    PTCA,SNA,CML,UDDC,TSCST,HDCS,MASC,AOAST,APAD,LIVER,MCCE,EMYOCA,LUAS,CEMN,HNMASC,BLCLC,BIMT,MPTLD,\
    NHL,MNGLP,SMZL,CDRCC,LGCOS,CHM,RCSNOS,ISTAD,GHCD,SCEMU,MCL,MBT,PCLBCLLT,CMML0,AMLGATA2MECOM,UTUC,\
    GBAD,STSC,SCB,UMEC,SCCRCC,CMLBCRABL1,EPIS,SCOAH,MCHSCNS,AMLNOS,BTOV,ETANTR,BPT,SAAD,MDSSLD,FDCS,\
    RGNT,WPSCC,ACA,ECD,PAC,SEBVTLC,UAS,LGFMS,PB,SEM,HSTCL,DIFG,STAD,PCALCL,DTE,TESTIS,CAEXPA,POCA,\
    UESL,LCIS,ANSC,UTERUS,SCLG,OS,PMF,MPALTNOS,SARCL,ODG,EP,ESS,NFIB,MRLS,GCTB,GS,GIST,MZL,MIXED,\
    BRAIN,PCM,SPIR,UCA,THYROID,MPT,MLNPDGFRB,PSTT,BRCNOS,PPB,OSMCA,BFN,LCS,SCOS,SMAHN,ANKL,AA,GBASC,\
    LRCHL,ATRT,CLPDNK,PRNE,GB,LM,FT,GNG,DFL,AMLMD,MPALBNOS,SARCNOS,FA,MHCD,BLAD,EOV,STMYEC,AGNG,\
    RCYC,LYP,MSCC,SOFT_TISSUE,PTHC,CEAS,THYMUS,UCEC,BLLBCRABL1L,CSCLC,PAMPCA,MLNPCM1JAK2,JMML,\
    MLNER,PCGDTCL,ET,PHM,CPC,SKCM,ACC,ADNOS,GCCAP,TRCC,SCCNOS,SCSRMS,CACC,IHM,HGBCL,PLEURA,\
    VOEC,MDSU,PPTID,CCOV,CUP,GNOS,RAML,HVLL,OTHER,OCSC,SDRPL,ANGS,GEJ,CPP,IVBCL,NLPHL,THYC,\
    ILC,PCNSL,HL,LAMN,PBT,ULM,ALCL,AMPCA,RNET,ALAL,DNT,GRCT,PENIS,ACYC,VPSCC,MCCHL,CEEN,CHS,\
    BLLTCF3PBX1,BEC,HGGNOS,CPT,HPHSC,OAT,LIHB,MNGT,IMPTLD,CCBOV,DIPG,USMT,PMFPES,ULMS,GRC,\
    OUSARC,DLBCLCI,CENE,LUCA,TEOS,MBEN,THYM,ACLG,THFO,EPDCA,ADMA,PECOMA,MDSEB1,UUS,BLLNOS,\
    ACCC,SKCN,BLLKMT2A,MMBC,GN,SCGBM,PCAECTCL,PGNT,ESCC,SBMOV,IMMC,MS,RAS,SBWDNET,MF,LIAS,\
    AMLNPM1,UCCC,AMLRUNX1RUNX1T1,VYST,HPCCNS,OUTT,IDCS,LGNET,PLEMESO,HCCIHCH,BREAST,LIAD,\
    MDSMPNRST,ETT,IMT,PDC,AMLRUNX1,PV,SBC,CAIS,MXOV,MUP,JSCB,UCU,HGNES,ISM,MATPL,MPNST,PTPR,\
    PCATCL,PCFCL,DA,CMML1,ABL,MYCHS,CCSK,MCS,PLBMESO,MYEC,ITPN,MIDD,CEGCC,CEMU,MPNU,EPMT,GCT,\
    PEOS,ACML,SKAC,DCS,CCE,EMBC,MDSRSMD,BCC,ISFN,UPECOMA,MPRDS,SCRMS,MCN,AMBL,CHGL,RCC,MGUSIGM,\
    BRCANOS,VULVA,SSRCC,DASTR,MSTAD,ETPLL,PTES,PPM,VSC,OOVC,LCH,PSEC,AM'
    
    cancers = list(map(lambda x: x.strip(), cancers.split(',')))
    # check that the cancer type is a valid cancer type in OncoTree
    if cancer_code not in cancers:
        raise ValueError('ERROR. Cancer type is not a valid cancer code defined in OncoTree')
   
    
    
def copy_segmentation_data(cbiodir, suppdir):
    '''
    (str, str) -> None
    
    Copy file processed segmentation file data_segments.txt from the cbioportal
    import folder cbiodir to the supplementary folder suppdir and rename it to data_segments.seg
    
    Parameters
    ---------
    - cbiodir (str): Path to the cbioportal import folder
    - suupdir (str): Path the supplementary folder 
    '''
    
    if os.path.isfile(os.path.join(cbiodir, 'data_segments.txt')):
        shutil.copy(os.path.join(cbiodir, 'data_segments.txt'), os.path.join(suppdir, 'data_segments.seg'))


def get_token(token_file):
    '''
    (str) -> str
    
    Returns the oncoKb token for variant and CNA annotation
    
    Parameters
    ----------
    - token_file (str): File containing the oncoKb token
    '''
    
    # get oncokb token
    infile = open(token_file)
    oncokb_token = infile.read().rstrip()
    infile.close()
    return oncokb_token



def make_import_folder(args):
    '''
    (list) -> None
    
    Generate folder cbioportal_import_folder with metadata and processed data files
    to upload to cBioPortal
        
    Parameters
    ----------
    - config (str): Path to the config file
    '''
    
    # parse config file
    config = configparser.ConfigParser(allow_no_value=True)
    config.read(args.config)
    # check config
    check_configuration(config)
    print('read and checked config')
    
    # extract variables from config
    ProcMAF, ProcCNA, ProcRNA, ProcFusion, token, enscon_hg38, enscon_hg19, entcon, genebed_hg38, genebed_hg19, genelist, oncolist = extract_resources_from_config(config)
    mapfile, outdir, study, center, cancer_code = extract_options_from_config(config)
    gain, amplification, heterozygous_deletion, homozygous_deletion, minfusionreads = extract_parameters_from_config(config)
    depth_filter, alt_freq_filter, gnomAD_AF_filter, tglpipe, filter_variants, filter_indels = extract_filters_from_config(config)
    print('extracted variables from config')
    
    # check that at least 1 maf file exists in map file
    mafs = extract_files_from_map(mapfile, 'maf')
    if len(mafs) == 0:
        raise ValueError('ERROR. At least 1 MAF file is required')
    print('checked for MAF files') 

    # check that cancer type is correctly defined
    check_cancer_type(cancer_code)
    print('checked cancer code')

    # create output directory and output sub-folders. remove output directory if it exists
    cbiodir, casedir, suppdir = create_output_directories(outdir)
    print('created output directories')

    # create input directories for each file type from map file    
    create_input_directories(outdir, mapfile)
    print('created input directories')
    
    # check for genome version in the first mutation file encountered
    genome = get_genome_version(mapfile)
    print('got reference genome from maf file: {0}'.format(genome))
    
    # get genome specific variables
    if genome == 'hg38':
        enscon = enscon_hg38
        genebed = genebed_hg38
    elif genome == 'hg19':
        enscon = enscon_hg19
        genebed = genebed_hg19
    
    # write meta study and clinical files
    write_meta_study(os.path.join(cbiodir, 'meta_study.txt') , study, genome, cancer_code)
    write_meta_clinical(cbiodir, study, 'sample')
    write_meta_clinical(cbiodir, study, 'patient')
    print('wrote study and clinical metadata')    
    
    # write cases
    write_cases(os.path.join(casedir, 'cases_sequenced.txt'), study, mapfile, 'seq')
    write_cases(os.path.join(casedir, 'cases_rna_seq_mrna.txt'), study, mapfile, 'rna')
    write_cases(os.path.join(casedir, 'cases_cna.txt'), study, mapfile, 'cna')
    write_cases(os.path.join(casedir, 'cases_cnaseq.txt'), study, mapfile, 'cna_seq')
    write_cases(os.path.join(casedir, 'cases_3way_complete.txt'), study, mapfile, 'cna_seq_rna')
    write_cases(os.path.join(casedir, 'cases_fusion.txt'), study, mapfile, 'fusion')
    print('wrote cases')

    # write patient and sample clinical information
    write_minimal_clinical_information(os.path.join(cbiodir, 'data_clinical_patients.txt'), mapfile, 'patient', center)
    # get the user defined sample clinical information
    write_minimal_clinical_information(os.path.join(cbiodir, 'data_clinical_samples.txt'), mapfile, 'sample', center, sample_info = args.clinical)
    print('wrote clinical information')
    
    # write clinical input file for oncokb-annotator
    write_clinical_oncokb(os.path.join(suppdir, 'oncokb_clinical_info.txt'), mapfile, cancer_code)
    print('wrote oncoKb clinical information') 
        
    # link files
    for i in ['maf', 'seg', 'gep', 'fus']:
        link_files(outdir, mapfile, i)
    print('linked files to input directories')
    
    # concatenate input files
    # concatenate mafs    
    mafdir = os.path.join(outdir, 'mafdir')
    assert os.path.isdir(mafdir)
    # concatenate all maf files into a plain text file
    concatenate_maf_files(mafdir, os.path.join(mafdir, 'all_mutations.maf.txt'))
    print('concatenated maf files')
    
    # check if seg dir and seg files exist
    segdir = os.path.join(outdir, 'segdir')
    if os.path.isdir(segdir):
        segfile = os.path.join(segdir, 'input.seg.txt')
        # concatenate seg files if they exist
        concatenate_seg_files(segdir, segfile)
    else:
        segfile = ''

    # check if fusion dir and fusion files exist
    fusdir = os.path.join(outdir, 'fusdir')
    if os.path.isdir(fusdir):
        fusfile = os.path.join(fusdir, 'input.fus.txt')
        # concatenate fusion files, if they exist
        concatenate_fusion_files(fusdir, fusfile)
    else:
        fusfile = ''
    
    # extract and concatenate fpkm from gep files
    gepdir = os.path.join(outdir, 'gepdir')
    if os.path.isdir(gepdir):
        gepfile = os.path.join(gepdir, 'input.fpkm.txt')
        # write fpkm to file
        concatenate_fpkm_from_gep_files(gepdir, gepfile)
    else:
        gepfile = ''
    
    
    # filter maf files and write metadata
    # define maffile, output of MafAnnotator
    maffile = os.path.join(mafdir, 'input.maf.txt')
    # filter mutations and indels if option is activated
    if filter_variants:
        total, kept = filter_mutations(os.path.join(mafdir, 'all_mutations.maf.txt'), os.path.join(mafdir, 'filtered.mutations.txt'), depth_filter, alt_freq_filter, gnomAD_AF_filter)
        print("before mutations filtering: ", total)
        print("after mutations filtering: ", kept)
        print('filtered variants')
    # filter indels if option is activated
    if filter_indels:
        # check if variants are filtered
        if filter_variants:
            maf_to_filter = os.path.join(mafdir, 'filtered.mutations.txt')
            maf_filtered = os.path.join(mafdir, 'filtered.mutations.indels.txt')
        else:
            maf_to_filter = os.path.join(mafdir, 'all_mutations.maf.txt')
            maf_filtered = os.path.join(mafdir, 'filtered.indels.txt')
        # output file for MafAnnotator
        maf_input_annotation = maf_filtered
        total, kept = remove_indels(maf_to_filter, maf_filtered)
        print("before indel filtering: ", total)
        print("after indel filtering: ", kept)
        print('filtered indels')
    else:
        if filter_variants:
            maf_input_annotation = os.path.join(mafdir, 'filtered.mutations.txt')    
        else:
            print('kept all variants')
            maf_input_annotation = os.path.join(mafdir, 'all_mutations.maf.txt')
    
    # get oncokb token
    oncokb_token = get_token(token)
    
    # annotate mafs with oncokb-annotate        
    maf_annotation = subprocess.call('MafAnnotator -i {0} -o {1} -c {2} -b {3}'.format(maf_input_annotation, maffile, os.path.join(suppdir, 'oncokb_clinical_info.txt'), oncokb_token), shell=True)
    # check exit code
    if maf_annotation:
        sys.exit('Error when running MafAnnotator.')
    else:
        print('Annotated variants with MafAnnotator')
    
    # generate mutations data
    process_mutations(maffile, tglpipe, ProcMAF, outdir)
    # write metadata file
    write_metadata(os.path.join(cbiodir, 'meta_mutations_extended.txt'), study, 'maf', genome)
    print('wrote mutations metadata')
    
    # generate CNA data and metadata files if input segmentation file exists
    if segfile:
        # generate metadata files
        # write cna metadata
        write_metadata(os.path.join(cbiodir, 'meta_CNA.txt'), study, 'discrete', genome)    
        write_metadata(os.path.join(cbiodir, 'meta_log2CNA.txt'), study, 'log2-value', genome)    
        write_metadata(os.path.join(cbiodir, 'meta_segments.txt'), study, 'seg', genome) 
        print('wrote CNA metadata files')
        if genelist:
            print('Restricting CNAs to the list of genes provided in {0}'.format(genelist))
        # generate data files
        process_cna(segfile, genebed, oncolist, gain, amplification, heterozygous_deletion, homozygous_deletion, ProcCNA, outdir, genelist)
        print('wrote CNA data files')
    # generate expression data and metadata file if input file exists
    if gepfile:
        # write metadata files
        write_metadata(os.path.join(cbiodir, 'meta_expression.txt'), study, 'expression', genome)
        write_metadata(os.path.join(cbiodir, 'meta_expression_zscores.txt'), study, 'zscore', genome)
        print('wrote expression metadata files')
        # write all samples with rna data to file 
        list_gep_samples(gepdir, os.path.join(outdir, 'gep_study.list'))
        # generate expression data files
        process_rna(gepfile, enscon, genelist, ProcRNA, outdir)
        print('wrote expression data files')        
    # generate fusion data and metadata if inout file exists
    if fusfile:
        # write fusion metadata
        write_metadata(os.path.join(cbiodir, 'meta_fusions.txt'), study, 'fusion', genome)
        print('wrote fusion metadata')
        # generate fusion data files
        process_fusion(fusfile, entcon, minfusionreads, ProcFusion, outdir)
        print('wrote fusion data files')
    
    # annonate CNA files with oncoKb for supplementary interpretation data
    # check that CNA data file is generated
    if os.path.isfile(os.path.join(suppdir, 'data_CNA_short.txt')):
        cna_annotation = subprocess.call('CnaAnnotator -i {0} -o {1} -c {2} -b {3}'.format(os.path.join(suppdir, 'data_CNA_short.txt'), os.path.join(suppdir, 'data_CNA_oncoKB.txt'), os.path.join(suppdir, 'oncokb_clinical_info.txt'), oncokb_token), shell = True)  
        if cna_annotation:
            sys.exit('Error when running CnaAnnotator.')
        else:
            print('Annotated CNAs with CnaAnnotator')
    else:
        print('WARNING. File {0} does not exist. Skipping CNA annotation'.format(os.path.join(suppdir, 'data_CNA_short.txt')))

    print('Success! Data in the cbioportal import folder is ready for upload.')



def import_cbioportal_project(args):
    '''
    (list) -> None
    
    Import a study into the GSI cBioPortal instance
    
    Parameters
    ----------
    - folder (str): Path to the cbioportal import folder with data to upload
    - key (str): Path to the file with cbioportal authentification key
    - user (str): User name. Default is ubuntu
    '''

    
    if args.key:
        key = '-i ' + args.key
    else:
        key = args.key
    
    base_folder = os.path.basename(os.path.abspath(args.folder))
    
    # create directory on cbioportal server 
    copying_message = 'Copying import folder to GSI cBioPortal instance' 
    print(copying_message + '\n' + len(copying_message) * '=' + '\n\n')
        
    unique_id = (uuid.uuid1()).int
    new_dir = "~/gsi/{0}.{1}".format(base_folder, unique_id)
    exit_code =  subprocess.call("ssh {0} ubuntu@cbio ' mkdir {1}'".format(key, new_dir), shell=True)
    if exit_code:
        sys.exit('Could not create folder on GSI cBioPortal instance')
    else:
        print('Created folder {0} on GSI cBioPortal instance'.format(os.path.basename(new_dir)))
    
    # copy data over to directory on cbioportal
    import_files = [os.path.join(args.folder, i) for i in os.listdir(args.folder)]
    for i in import_files:
        exit_code  = subprocess.call('scp -r {0} {1} ubuntu@cbio:{2}'.format(key, i, new_dir), shell=True)
        if exit_code:
            sys.exit('Could not copy file {0} to folder {1} on GSI cBioPortal instance'.format(os.path.basename(i), new_dir))
        else:
            print('copied file {0} to {1} on GSI cBioPortal instance'.format(os.path.basename(i), new_dir))
        
    import_message = 'Importing study to GSI cBioPortal instance'
    print(import_message + '\n' + len(import_message) * '=' + '\n\n')
    
    # Import study with import_study.sh. Precondition. cbioportal is installed with docker image
    cmd = "ssh {0} -t ubuntu@cbio ' /home/ubuntu/import_study.sh {1}'".format(key, new_dir)
    output = subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True).decode('utf-8').rstrip()
    
    print(output + '\n\n\n')


if __name__ == '__main__':

    # create parser
    parser = argparse.ArgumentParser(prog = 'pycBio.py', description='A script to generate the cbiportal import folder')
    # create subparser
    subparsers = parser.add_subparsers(help='sub-command help', dest='subparser_name')
    
    # generate import folder
    g_parser = subparsers.add_parser('generate', help="Generate cbio import folder")
    g_parser.add_argument('-cf', '--Config', dest='config', help='Path to the config file', required = True)
    g_parser.add_argument('-cl', '--Clinical', dest='clinical', help='Path to the sample clinical file')
    g_parser.set_defaults(func=make_import_folder)
    
    # import folder to gsi cbioportal instance
    i_parser = subparsers.add_parser('import', help="Import cbioportal folder to GSI cBioPortal instance")
    i_parser.add_argument('-f', '--folder', dest='folder', help='Path to cbioportal import folder', required = True)
    i_parser.add_argument("-k", "--key", dest='key', default = '', help="Path to the cBioPortal Key")
    i_parser.add_argument("-u", "--user", dest= 'user', default='ubuntu', help="The linux distribution. Default is ubuntu")
    i_parser.set_defaults(func=import_cbioportal_project)
 
    # get arguments from the command line
    args = parser.parse_args()
    # pass the args to the default function
    args.func(args)
