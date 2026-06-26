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
import numpy as np
import json
import glob
import zipfile
import pandas as pd 
import gc 
import warnings
from rpy2 import robjects 
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
from rpy2.robjects.conversion import localconverter
from scipy.stats import zscore



cntools = importr('CNTools')



def extract_files_from_map(mapfile, data_type):
    '''
    (str, str) -> list
    
    Returns a list of input files from the mapping file
        
    Parameters
    ----------
    - mapfile (str): Mapping file (map.csv) that contains paths to maf, seg, gep and mavis files    
    - data_type (str): File type to link in their respective folders.
                       Accepted values: snv, expr, fus, and cna
    '''

    # create input directories for each file type from map file    
    infile = open(mapfile)
    content = infile.read().rstrip().split('\n')
    infile.close()
    for i in range(len(content)):
        content[i] = list(map(lambda x: x.strip(), content[i].split(',')))
    
    # make a list of samples and files for which files exist
    if data_type == 'snv':
        # get the maf files
        j = 2
    elif data_type == 'cna':
        # get the cna files
        j = 3
    elif data_type == 'expr':
        # get the rna files
        j = 4
    elif data_type == 'fus':
        # get the fusion files
        j = 5
    
    files = [i[j] for i in content if i[j].upper() != 'NA' and os.path.isfile(i[j])]
    return files


def extract_purple_ploidy_purity(purple_purity_file):
    '''
    (str) -> (float, float)
    
    Returns a tuple with purity and ploidy extracted from file .purple.purity.tsv
 
    Parameters
    ----------
    - purple_purity_file (str): Path to the Purple .purple.purity.tsv file
    '''
    
    infile = open(purple_purity_file)
    header = infile.readline().rstrip().split('\t')
    content = infile.readline().rstrip().split('\t')
    infile.close()
    
    purity = float(content[header.index('purity')])
    ploidy = float(content[header.index('ploidy')])

    return purity, ploidy
    

def convert_purple_to_seg(somatic_file, purity, ploidy, sample_name):
    '''
    (str, str) -> list
    
    Returns a list of fields extracted from the purple somatic file
    corresponding to the sequenza segmentation output file 
    
    Fields are mapped and renamed as:
    chromosome -> chrom
    start -> loc.start
    end -> loc.end 
    bafCount -> num.mark
    copyNumber -> seg.mean adjusted with purity and ploidy  
    
    Parameters
    ----------
    - somatic_file (str): Path to the purple "purple.cnv.somatic.tsv" file
    - sample_name (str): Name of the sample
    - purity (float): Purity estimated from purple
    - ploidy (float): Ploidy estimated from purple
    '''
    
    L = [['ID', 'chrom', 'loc.start', 'loc.end', 'num.mark', 'seg.mean']]
    
    infile = open(somatic_file)
    header = infile.readline().rstrip().split('\t')
    for line in infile:
        line = line.rstrip()
        if line:
            line = line.split('\t')
            # remove records when bafcount is 0
            if int(line[header.index('bafCount')]) != 0:
                # compute seg.mean       
                copyNumber = float(line[header.index('copyNumber')])
                # only includes seg.mean when it is defined
                with warnings.catch_warnings():
                    warnings.simplefilter("error", RuntimeWarning)
                    try:
                        adjusted_copy_number = np.log2(1 + (purity *(copyNumber - ploidy)/ploidy))
                    except:
                        adjusted_copy_number = 'NA'
                    finally:
                        if adjusted_copy_number != 'NA':
                            L.append([sample_name, line[header.index('chromosome')],
                                      line[header.index('start')],
                                      line[header.index('end')],
                                      line[header.index('bafCount')],
                                      str(adjusted_copy_number)])
    infile.close()
    
    return L


def generate_purple_segmentation(somatic_file, purple_purity_file, sample_name, outputfile):
    '''
    (str, str, str) -> None
    
    Convert the purple somatic file to a segmentation file and write to outputfile
    
    Parameters
    ----------
    - somatic_file (str): Path to the purple "purple.cnv.somatic.tsv" file
    - purple_purity_file (str): Path to the Purple .purple.purity.tsv file
    - sample_name (str): Name of the sample
    - outputfile (str): Path to the outputfile
    '''

    purity, ploidy = extract_purple_ploidy_purity(purple_purity_file)
    
    L = convert_purple_to_seg(somatic_file, purity, ploidy, sample_name)
    newfile = open(outputfile, 'w')
    for i in L:
        newfile.write('\t'.join(i) + '\n')
    newfile.close()


def split_column_take_max(df, columns):
    '''
    Splits columns by ';' and takes the maximum value from a list of columns
    Parameters
    ----------
    - df (pd.DataFrame) : The input data that contains the columns to be processed
    - columns (list) : A list of columns names (str) to be processed
    Returns
    -------
    - df (pd.DataFrame) : The modified dataframe 
    '''
    df[columns] = df[columns].fillna(0).replace('', 0)

    for column in columns:
        df[column] = df[column].apply(
             lambda x: (
                max([int(t) if t.strip().lstrip('+-').isdigit() else 0 for t in str(x).split(';')])
                if isinstance(x, str) and ';' in str(x)
                else (0 if x is None or (isinstance(x, str) and x.strip().lower() in {'none', 'na', 'nan', 'null', ''}) else int(x)))        
        )
    
    for column in columns:
        df[column] = pd.to_numeric(df[column], errors='coerce')  

    return df


def preProcFus(datafile, readfilt, entrfile):
    '''
    Preprocesses the fusion input data and formats it for CBioPortal output. 
    Parameters
    -----------
    - datafile (str) : Path to the input fusion data file
    - readfilt (int) : Minimum number of reads for fusion calls
    - entrfilw (str) : Path to the Entrez gene ID file
    Returns
    --------
    - df_cbio (pd.DataFrame) : Processed and formatted data
    '''
    data = pd.read_csv(datafile, sep="\t")
    entr = pd.read_csv(entrfile, sep="\t")

    # reformat the filtering columns to split and take the max value within cell
    columns = ['contig_remapped_reads', 'flanking_pairs', 'break1_split_reads', 'break2_split_reads', 'linking_split_reads']
    data = split_column_take_max(data, columns)

    # add a column which pulls the correct read support column    
    data['read_support'] = data.apply(
    lambda row: row['contig_remapped_reads'] if 'contig' in row['call_method'] 
    else row['flanking_pairs'] if 'flanking reads' in row['call_method'] 
    else (row['break1_split_reads'] + row['break2_split_reads'] + row['linking_split_reads']) if 'split reads' in row['call_method'] 
    else 0, 
    axis=1)

    # filter by minimum read support
    data = data[data['read_support'] > readfilt]

    # sort descending read support 
    data = data.sort_values(by='read_support', ascending=False)

    # get unique fusions for each sample
    data['gene1_aliases'] = data['gene1_aliases'].fillna('')
    data['gene2_aliases'] = data['gene2_aliases'].fillna('')

    data['fusion_tuples'] = data.apply(
    lambda row: '-'.join(sorted([str(x) if x else 'None' for x in [row['gene1_aliases'], row['gene2_aliases']]])), axis=1)

    # add index which is sample, tuple
    data['index'] = data['Sample'] + data['fusion_tuples']

    # deduplicate
    data_dedup = data.drop_duplicates(subset='index', keep='first')

    # gene1 should not equal gene2
    data_dedup = data_dedup[data_dedup['gene1_aliases'] != data_dedup['gene2_aliases']]

    # merge in entrez gene ids
    data_dedup = pd.merge(data_dedup, entr, how='left', left_on='gene1_aliases', right_on='Hugo_Symbol')
    data_dedup = pd.merge(data_dedup, entr, how='left', left_on='gene2_aliases', right_on='Hugo_Symbol', suffixes=('.x', '.y'))

    # add some missing columns
    data_dedup['DNA_support'] = 'no'
    data_dedup['RNA_support'] = 'yes'
    data_dedup['Center'] = 'TGL'
    data_dedup['Frame'] = 'frameshift'
    data_dedup['Fusion_Status'] = 'unknown'

    # write out the nice header 
    header = ['Hugo_Symbol', 'Entrez_Gene_Id', 'Center', 'Tumor_Sample_Barcode', 'Fusion', 'DNA_support', 'RNA_support', 'Method', 'Frame', 'Fusion_Status']

    # get left gene data
    col_left = ['gene1_aliases', 'Entrez_Gene_Id.x', 'Center', 'Sample', 'fusion_tuples', 'DNA_support', 'RNA_support', 'tools', 'Frame', 'Fusion_Status']
    data_left = data_dedup[col_left]
    data_left.columns = header

    # get right gene data
    col_right = ['gene2_aliases', 'Entrez_Gene_Id.y', 'Center', 'Sample', 'fusion_tuples', 'DNA_support', 'RNA_support', 'tools', 'Frame', 'Fusion_Status']
    data_right = data_dedup[col_right]
    data_right.columns = header

    data_left.loc[:, 'Hugo_Symbol'] = data_left['Hugo_Symbol'].replace("", float('nan')).fillna(data_right['Hugo_Symbol'])
    data_left.loc[:, 'Entrez_Gene_Id'] = data_left['Entrez_Gene_Id'].fillna(data_right['Entrez_Gene_Id'])
    data_right.loc[:, 'Hugo_Symbol'] = data_right['Hugo_Symbol'].replace("", float('nan')).fillna(data_left['Hugo_Symbol'])
    data_right.loc[:, 'Entrez_Gene_Id'] = data_right['Entrez_Gene_Id'].fillna(data_left['Entrez_Gene_Id'])

    # append it all together
    df_cbio = pd.concat([data_left, data_right])
    df_cbio['Entrez_Gene_Id'] = (    
        pd.to_numeric(df_cbio['Entrez_Gene_Id'], errors='coerce').fillna(0).astype(int)
    )

    # remove rows where gene is not known (this still keeps the side of the gene which is known)
    df_cbio = df_cbio.dropna()

    return df_cbio



def preProcCNA(segfile, genebed, gain, amp, htz, hmz, oncolist, genelist=None):
    '''
    Processes CNA data by applying thresholds and performing gene-level segmentation.
    Parameters
    ----------
    - segfile (str) : Path to the input concatenated segmentation file from Sequenza
    - genebed (str) : Path to a tab-delimited bed file defining the genomic positions of canonical genes
    - gain (float) : Threshold for gains in CNA data
    - amp (float) : Threshold for amplification in CNA data
    - htz (float) : Threshold for heterozygous deletion in CNA data
    - hmz (float) : Threshold for homozygous deletion in CNA data
    - oncolist (str) : Path to a tab-delimited file containing the list of cancer genes
    - genelist (str or None) : (Optional) Path to a file containing a list of Hugo Symbols to filter the genes, or None if no filtering is required
    Returns
    -------
    - segData (pd.DataFrame) : Dataframe containing the processed segmentation data
    - df_cna (pd.DataFrame) : Dataframe with gene-level log2 copy number alterations
    - df_cna_thresh (pd.DataFrame) : Dataframe with thresholded CNA values (5-state matrix)
    '''
    # check if genelist is None when preProcCNA.py is called by pycbio.py and genelist is omitted
    if genelist:
        print('genelist is used during CNA processing')
    else:
        print('genelist is not used during CNA processing')
        
    
    # read oncogenes
    oncogenes = pd.read_csv(oncolist, sep='\t')

    # set thresholds
    print('setting thresholds')
    gain, amp, htz, hmz = float(gain), float(amp), float(htz), float(hmz)

    # small fix segmentation data
    segData = pd.read_csv(segfile, sep='\t')
    segData['chrom'] = segData['chrom'].str.replace('chr', '')
    # remove NA
    segData.dropna(inplace=True)

    # get the gene info
    print('getting gene info')
    geneInfo = pd.read_csv(genebed, sep='\t')

    # convert pandas dataframes to R dataframes
    with localconverter(robjects.default_converter + pandas2ri.converter):
        segData_r = robjects.conversion.py2rpy(segData)
        print('converted segData to R df')
        geneInfo_r = robjects.conversion.py2rpy(geneInfo)
        print('converted geneInfo to R df')
        
        # make CN matrix gene level
        print('converting seg')
        cnseg = cntools.CNSeg(segData_r)
        print('get segmentation by gene')
        rdByGene = cntools.getRS(cnseg, by='gene', imput=False, XY=False, geneMap=geneInfo_r, what='median', mapChrom='chrom', mapStart='start', mapEnd='end')
        print('get reduced segmentation data')
        reducedseg_df = cntools.rs(rdByGene)

        # convert data from R back to Py
        reducedseg = robjects.conversion.rpy2py(reducedseg_df)    
        print('converted reducedSegData back to pandas df')

    # some reformatting and return log2cna data
    df_cna = reducedseg.iloc[:, [4] + list(range(5, reducedseg.shape[1]))]
    df_cna = df_cna.drop_duplicates(subset=[df_cna.columns[0]])
    df_cna.columns = ['Hugo_Symbol'] + list(df_cna.columns[1:])

    # set thresholds and return 5-state matrix
    print('thresholding cnas')
    df_cna_thresh = df_cna.copy()
    df_cna_thresh.iloc[:, 1:] = df_cna_thresh.iloc[:, 1:].apply(pd.to_numeric)

    # threshold data
    for col in df_cna_thresh.columns[1:]:
        df_cna_thresh[col] = df_cna_thresh[col].apply(lambda x: 2 if x > amp
                                                else (-2 if x < hmz
                                                      else (1 if gain < x <= amp
                                                            else (-1 if hmz <= x < htz
                                                                  else 0)
                                                        )
                                                    )
                                                )
    
    # fix rownames of log2cna data
    df_cna.set_index('Hugo_Symbol', inplace=True)
    df_cna = df_cna.round(4)
    df_cna_thresh.set_index(df_cna_thresh.columns[0], inplace=True)

    # subset of oncoKB genes
    df_cna_thresh_onco = df_cna_thresh.loc[df_cna_thresh.index.isin(oncogenes.index)]

    # subset if gene list is given
    if genelist is not None:
        with open(genelist, 'r') as file:
            keep_genes = [line.strip() for line in file]
            
        df_cna = df_cna.loc[df_cna.index.isin(keep_genes)]
        df_cna_thresh = df_cna_thresh.loc[df_cna_thresh.index.isin(keep_genes)]
      
    return segData, df_cna, df_cna_thresh



def addVAFtoMAF(maf_df, alt_col, dep_col, vaf_header):
    '''
    Adds a VAF column to a MAF DataFrame.
    Parameters
    ----------
    - maf_df (pd.DataFrame): The MAF DataFrame containing mutation data.
    - alt_col (str): The name of the column representing the alternate allele count.
    - dep_col (str): The name of the column representing the depth. 
    - vaf_header (str): The name for the new column that will hold the VAF values.
    Returns
    -------
    - maf_df (pd.DataFrame) : The modified dataframe with the VAF column. 
    '''
    # print a warning if any values are missing (shouldn't happen), but change them to 0
    if maf_df[alt_col].isnull().any() or maf_df[dep_col].isnull().any():
        print('Warning! Missing values found in one of the count columns')
        maf_df[alt_col] = maf_df[alt_col].fillna(0)
        maf_df[dep_col] = maf_df[dep_col].fillna(0)
    
    # ensure factors end up as numeric
    maf_df[alt_col] = pd.to_numeric(maf_df[alt_col], errors='coerce')
    maf_df[dep_col] = pd.to_numeric(maf_df[dep_col], errors='coerce')

    # ensure position comes after alternate count field 
    bspot = maf_df.columns.get_loc(alt_col)
    maf_df.insert(bspot + 1, vaf_header, maf_df[alt_col] / maf_df[dep_col])

    # check for any NAs
    if maf_df[vaf_header].isnull().any():
        print('Warning! There are missing values in the new vaf column')
        maf_df[vaf_header] = maf_df[vaf_header].fillna(0)
    
    return maf_df
    

def procVEP(datafile):
    '''
    Processes the input MAF file, adding various computed columns, and applying multiple filters to prepare the data for further analysis.
    Parameters
    ----------
    - datafile (str): The file path to the input in tab-separated format.
    Returns
    -------
    - df_anno (pd.DataFrame) : The modified dataframe. 
    '''
    print("--- reading data ---")
    data = pd.read_csv(datafile, sep="\t")

    print('--- doing some formatting ---')

    # add vaf columns 
    print('add tumor_vaf')
    data = addVAFtoMAF(data, 't_alt_count', 't_depth', 'tumor_vaf')
    print('add normal_vaf')
    data = addVAFtoMAF(data, 'n_alt_count', 'n_depth', 'normal_vaf')

    # clear memory (important when the mafs are huge with millions and millions of lines)
    df_anno = data
    del data
    gc.collect()

    # add oncogenic yes or no columns 
    print('add oncogenic status')
    df_anno['oncogenic_binary'] = df_anno['oncogenic'].apply(lambda x: 'YES' if x in ['Oncogenic', 'Likely Oncogenic'] else 'NO')

    # add common_variant yes or no columns
    df_anno['ExAC_common'] = df_anno['FILTER'].apply(lambda x: 'YES' if 'common_variant' in x else 'NO')

    # add POPMAX yes or no columns 
    print('add population level frequency')
    gnomad_cols = ['gnomAD_AFR_AF', 'gnomAD_AMR_AF', 'gnomAD_ASJ_AF', 'gnomAD_EAS_AF', 'gnomAD_FIN_AF', 'gnomAD_NFE_AF', 'gnomAD_OTH_AF', 'gnomAD_SAS_AF']
    df_anno[gnomad_cols] = df_anno[gnomad_cols].fillna(0)
    df_anno['gnomAD_AF_POPMAX'] = df_anno[gnomad_cols].max(axis=1)

    # caller artifact filters 
    print('apply filters')
    df_anno['FILTER'] = df_anno['FILTER'].replace('clustered_events', 'PASS')
    df_anno['FILTER'] = df_anno['FILTER'].replace('common_variant', 'PASS')

    # some specific filter flags should be rescued if oncogenic (i.e. EGFR had issues here)
    print("rescue filter flags if oncogenic")
    df_anno['FILTER'] = df_anno.apply(
        lambda row: 'PASS' if row['oncogenic_binary'] == 'YES' and row['FILTER'] in ['triallelic_site', 'clustered_events;triallelic_site', 'clustered_events;homologous_mapping_event'] else row['FILTER'],
        axis=1
    )

    # Artifact Filter
    print('artifact filter')
    df_anno['TGL_FILTER_ARTIFACT'] = df_anno['FILTER'].apply(lambda x: 'PASS' if x == 'PASS' else 'Artifact')

    # ExAC Filter
    print('exac filter')
    df_anno['TGL_FILTER_ExAC'] = df_anno.apply(
        lambda row: 'ExAC_common' if row['ExAC_common'] == 'YES' and row['Matched_Norm_Sample_Barcode'] == 'unmatched' else 'PASS',
        axis=1
    )

    # gnomAD_AF_POPMAX Filter
    print('population frequency filter')
    df_anno['TGL_FILTER_gnomAD'] = df_anno.apply(
        lambda row: 'gnomAD_common' if row['gnomAD_AF_POPMAX'] > 0.001 and row['Matched_Norm_Sample_Barcode'] == 'unmatched' else 'PASS',
        axis=1
    )

    # VAF Filter
    print('VAF Filter')
    df_anno['TGL_FILTER_VAF'] = df_anno.apply(
        lambda row: 'PASS' if row['tumor_vaf'] >= 0.1 or 
        (row['tumor_vaf'] < 0.1 and row['oncogenic_binary'] == 'YES' and 
         ((row['Variant_Classification'] in ['In_Frame_Del', 'In_Frame_Ins']) or 
          row['Variant_Type'] == 'SNP')) 
        else 'low_VAF', axis=1
    )

    # Mark filters 
    print('Mark filters')
    df_anno['TGL_FILTER_VERDICT'] = df_anno.apply(
        lambda row: 'PASS' if row['TGL_FILTER_ARTIFACT'] == 'PASS' and 
                          row['TGL_FILTER_ExAC'] == 'PASS' and 
                          row['TGL_FILTER_gnomAD'] == 'PASS' and 
                          row['TGL_FILTER_VAF'] == 'PASS' 
                  else ';'.join([row['TGL_FILTER_ARTIFACT'], row['TGL_FILTER_ExAC'], 
                                 row['TGL_FILTER_gnomAD'], row['TGL_FILTER_VAF']]), axis=1
    )

    return df_anno




def preProcRNA(fpkmfile, enscon, genelist=None):
    '''
    Preprocesses RNA expression data by merging with gene symbol annotations, optionally subsetting by a gene list.
    Parameters
    ----------
    - fpkmfile (str): Path to the input concatenated FPKM data from RSEM workflow.
    - enscon (str): Path to a tab-delimited file with ENSEMBLE gene ID and Hugo_Symbol.
    - genelist (str, optional): Path to a file with list of Hugo Symbols to report in the final results. Default is None.
    Returns
    -------
    - df (pd.DataFrame): Preprocessed DataFrame with gene expression data.
    '''
    # check if genelist is None when preProcRNA.py is called by pycbio.py and genelist is omitted
    if genelist:
        print('genelist is used during RNA processing')
    else:
        print('genelist is not used during RNA processing')
        
    # read in data
    gep_data = pd.read_csv(fpkmfile, sep="\t")
    ens_conv = pd.read_csv(enscon, sep="\t", header=None)

    # rename columns
    ens_conv.columns = ["gene_id", "Hugo_Symbol"]

    # merge in Hugo's, re-order columns, deduplicate
    df = pd.merge(gep_data, ens_conv, on="gene_id", how="left")
    df = df.iloc[:, [-1] + list(range(1, df.shape[1] - 1))]
    df = df[~df.duplicated(subset=[df.columns[0]])]

    df.set_index("Hugo_Symbol", inplace=True)

    # subset if gene list is given
    if genelist is not None:
        with open(genelist, 'r') as file:
            keep_genes = [line.strip() for line in file]

        df = df[df.index.isin(keep_genes)]
    
    # return the data frame
    return df


def compZ(df):
    '''
    Computes row-wise z-scores for a given DataFrame and fills NaN values with zero.
    Parameters
    ----------
    - df (pd.DataFrame): Input DataFrame with gene expression data.
    Returns
    -------
    - df_zscore (pd.DataFrame): DataFrame with z-scores for each gene.
    '''
    # scale row-wise
    df_zscore = df.apply(lambda x: (x - x.mean()) / x.std(), axis=1)
    df_zscore = df_zscore.fillna(0)

    return df_zscore


def create_input_directories(outdir):
    '''
    (str) -> list
    
    Create sub-directories in outdir for each type of data file (even when data file are not present)
    and returns the list of directories    
            
    Parameters
    ----------
    - outdir (str): Path to the output directory where mafdir, sgedir, fusdir and gepdir folders are located 
    '''

    L = []

    # create input directories regardless of input data
    for i in ['maf', 'seg', 'gep', 'fus']:
        filedir = os.path.join(outdir, '{0}dir'.format(i))
        os.makedirs(filedir, exist_ok=True)
        L.append(filedir)
        
    return L
    

def write_meta_study(outputfile, study, project_name, description, genome, cancerType):
    '''
    (str, str, str, str) -> None
    
    Write file meta_study.txt in the cbioportal_import_data folder
    
    Parameters
    ----------
    - outputfile (str): Path to the outputfile 
    - study (str): Long study name, following the format: ACRONYM: Top-level-OncoTree, Concept (PI, Centre)
    - project_name: Short project name, field project_name in the configuration file 
    - description: Short description of the study
    - genome (str): Reference genome (hg19 or hg38)
    - cancerType (str): Cancer type as defined in http://oncotree.mskcc.org
    '''

    newfile = open(outputfile, 'w')
    L = ['cancer_study_identifier: {0}'.format(project_name),
         'description: {0}'.format(description),
         'groups: ',
         'name: {0}'.format(study),
         'reference_genome: {0}'.format(genome),
         'short_name: {0}'.format(project_name),
         'add_global_case_list: true',
         'type_of_cancer: {0}'.format(cancerType)]
    newfile.write('\n'.join(L))
    newfile.close()     
    

def write_meta_clinical(cbio_import_dir, project_name, data_type):
    '''
    (str, str, str) -> None
    
    Write file meta_clinical_patients.txt or meta_clinical_samples.txt in the cbioportal_import_data folder
    
    Parameters
    ----------
    - cbio_import_dir (str): Path to the cbioportal_import_data directory 
    - project_name (str): project_name in the config file
    - data_type(str): Sample or patient
    '''
    
    if data_type.lower() not in ['sample', 'patient']:
        raise ValueError('ERROR. Data type must be "sample" or "patient"')
    
    filename = 'meta_clinical_{0}s.txt'.format(data_type.lower())
    outputfile = os.path.join(cbio_import_dir, filename)
    newfile = open(outputfile, 'w')
    L = ['cancer_study_identifier: {0}'.format(project_name),
         'data_filename: {0}'.format(filename.replace('meta', 'data')),
         'datatype: {0}_ATTRIBUTES'.format(data_type.upper()),
         'genetic_alteration_type: CLINICAL']
    newfile.write('\n'.join(L) + '\n')
    newfile.close()
    

def check_genome_version(mapfile, genome):
    '''
    (str, str) -> None
    
    Check if the genome version (hg19 or hg38) found in each MAF file listed in the map file,
    if they exist, are the same as the expected genome ref from the config file
        
    Parameters
    ----------
    - mapfile (str): Mapping file (map.csv) that contains paths to maf, seg, gep and mavis files    
    - genome (str): Genome reference from the config file
    '''
    
    # get mafs listed in mapfile
    mafs = extract_files_from_map(mapfile, 'snv')
    L = set()
    if mafs:
        for i in mafs:
            # grab genome column in maf. skipping header and commented lines
            infile = gzip.open(i, 'rt')
            for line in infile:
                if not line.startswith('#') and 'Hugo_Symbol' not in line:
                    line = line.rstrip().split('\t')
                    L.add(line[3])              
            infile.close()
    
    # convert genome identifier
    L = ';'.join(list(L))
    if L == "GRCh38":
        genomev="hg38"
    elif L == "GRCh37":
        genomev="hg19"
    else:
        genomev = L
        
    if genomev:
        if genome != genomev:
            raise ValueError('ERROR. Reference in MAF file does not match reference in config: {0} vs {1}'.format(genome, genomev))
        else:
            print('validated reference genome: {0}'.format(genome))



def write_cases(outputfile, project_name, mapfile, data_type):
    '''
    (str, str, str, str) -> None
    
    Write list of samples profiled for data type 
            
    Parameters
    ----------
    - outputfile (str): Path to the outputfile
    - project_name (str): Field project_name in configuration file
    - mapfile (str): Mapping file (map.csv) that contains paths to maf, seg, gep and mavis files    
    - data_type (str): The type of data considered.
                       Values accepted are: seq, rna, cna, cna_seq, cna_seq_rna, sv
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
        stable_id = '{0}_sequenced'.format(project_name)
    elif data_type == 'sv':
        # make a list of samples for which SV files are available
        samples = [i.split(',')[1] for i in content if i.split(',')[5].upper() != 'NA']
        name = 'Samples profiled for structural variants'
        description = 'This is this case list that contains all samples that are profiled for structural variants.'
        stable_id = '{0}_sv'.format(project_name)
    elif data_type == 'rna':
        # make a list of samples for which rsem files are available
        samples = [i.split(',')[1] for i in content if i.split(',')[4].upper() != 'NA']
        name = 'Samples profiled for rnaseq'
        description = 'This is this case list that contains all samples that are profiled for rnaseq.'
        stable_id = '{0}_rna_seq_mrna'.format(project_name)    
    elif data_type == 'cna':
        # make a list of samples for which cna files are available
        samples = [i.split(',')[1] for i in content if i.split(',')[3].upper() != 'NA']
        name = 'Samples profiled for cnas'
        description = 'This is this case list that contains all samples that are profiled for cnas.'
        stable_id = '{0}_cna'.format(project_name)
    elif data_type == 'cna_seq':
        # make a list of samples for which cna and maf files are available
        samples = [i.split(',')[1] for i in content if i.split(',')[3].upper() != 'NA' and i.split(',')[2].upper() != 'NA' ]
        name = 'Samples profiled for cnas and sequencing'
        description = 'This is this case list that contains all samples that are profiled for mutations and cnas.'
        stable_id = '{0}_cnaseq'.format(project_name)
    elif data_type == 'cna_seq_rna':
        # make a list of samples for which cna and maf and rna files are available
        samples = [i.split(',')[1] for i in content if i.split(',')[2].upper() != 'NA' and i.split(',')[3].upper() != 'NA' and i.split(',')[4].upper() != 'NA']
        name = 'Samples profiled for all of mutation, cnas, and rnaseq'
        description = 'This is this case list that contains all samples that are profiled for mutations, cnas, and rnaseq.'
        stable_id = '{0}_3way_complete'.format(project_name)
    
    # write outputfile if samples exist
    if samples:
        newfile = open(outputfile, 'w')
        L = ['cancer_study_identifier: {0}'.format(project_name),
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



def get_clinical_fields(header):
    '''
    (list) -> list
    
    Returns  a non-redendant list of fields from the clinical sample file header
    
    Parameters
    ----------
    - header (list): List representation of the clinical sample file header
    '''
    
    header_columns = []
    for i in range(len(header)):
        if i in [0, 1, 4]:
            for j in header[i]:
                header_columns.append(j.replace('#', ''))
    header_columns = list(set(header_columns))
        
    return header_columns


def update_clinical_sample_header(sample_info, header):
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
    
    # make a list of column names
    column_names = list_column_names(sample_info)
    # map column name with data type
    data_types = map_column_data_type(sample_info)
    # make a list of existing column names in header
    header_columns = list(lambda x: x.lower(), get_clinical_fields(header))
    
    for i in column_names:
        if i.lower() not in header_columns:
            # add column name to header
            for k in [0, 1, 4]:
                header[k].append(i.upper())
            # add data type
            header[2].append(data_types[i])
            # add priority
            header[3].append('1')
    return header



def write_patient_minimal_clinical_information(outputfile, mapfile, centre):
    '''
    (str, str, str) -> None
    
    Write clinical files with minimal clinical information
    
    Parameters
    ----------    
    - outputfile (str): Path to the outputfile
    - mapfile (str): Mapping file that contains paths to maf, segments, expression and fusion files    
    - centre (str): Genomic centre (eg TGL, OICR)
    '''
    
    # make a list with sample names and libraries
    infile = open(mapfile)
    content = infile.read().rstrip().split('\n')
    infile.close()
    S = [i.split(',')[0:2] for i in content]
        
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
        
    for i in U:
        T.append('\t'.join(i + [''] * (len(T[0].split('\t')) - len(i))))
        
    newfile = open(outputfile, 'w')
    for i in T:
        newfile.write(i + '\n')
    newfile.close()



def write_sample_minimal_clinical_information(outputfile, mapfile, centre, sample_info = None):
    '''
    (str, str, str, str, dict | None, dict | None, list | None) -> None
    
    Write clinical files with minimal clinical information
    
    Parameters
    ----------    
    - outputfile (str): Path to the outputfile
    - mapfile (str): Mapping file (map.csv) that contains paths to maf, segments, expression and fusion files    
    - centre (str): Genomic centre (eg TGL, OICR)
    - sample_info (dict | None): Dictionary with patient and sample information. Populates the sample clinical file
    '''
    
    # make a list with sample names and libraries
    infile = open(mapfile)
    content = infile.read().rstrip().split('\n')
    infile.close()
    S = [i.split(',')[0:2] for i in content]
    
    # build header
    T = ['#Patient Identifier\tSample Identifier\tCosmic Signature\tPRIMARY SITE\tCANCER TYPE\tCLOSEST TCGA\tSAMPLE ANATOMICAL SITE\tSAMPLE PRIMARY OR METASTASIS\tTREATMENT STATUS\tPATHOLOGICAL REVIEW\tPRIOR CLINCAL TEST RESULTS\tMEAN COVERAGE\tPCT V7 ABOVE 80X\tPCT CALLABILITY\tSEQUENZA PURITY FRACTION\tSEQUENZA PLOIDY\tTMB PER MB\tHRD SCORE\tMSI STATUS',
         '#Patient Identifier\tSample Identifier\tCosmic Signature\tPRIMARY SITE\tCANCER TYPE\tCLOSEST TCGA\tSAMPLE ANATOMICAL SITE\tSAMPLE PRIMARY OR METASTASIS\tTREATMENT STATUS\tPATHOLOGICAL REVIEW\tPRIOR CLINCAL TEST RESULTS\tMEAN COVERAGE\tPCT V7 ABOVE 80X\tPCT CALLABILITY\tSEQUENZA PURITY FRACTION\tSEQUENZA PLOIDY\tTMB PER MB\tHRD SCORE\tMSI STATUS',
         '#STRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tNUMBER\tNUMBER\tNUMBER\tNUMBER\tNUMBER\tNUMBER\tNUMBER\tSTRING']
    T.append('#1' + ('\t1' * (len(T[0].split('\t')) -1)))     
    T.append('PATIENT_ID\tSAMPLE_ID\tCOSMIC_SIGS\tCANCER_TYPE\tCANCER_TYPE_DETAILED\tCLOSEST_TCGA\tSAMPLE_ANATOMICAL_SITE\tSAMPLE_PRIMARY_OR_METASTASIS\tTREATMENT_STATUS\tPATHOLOGICAL_REVIEW\tPRIOR_CLINCAL_TEST_RESULTS\tMEAN_COVERAGE\tPCT_V7_ABOVE_80X\tFRAC_CALLABILITY\tSEQUENZA_PURITY_FRACTION\tSEQUENZA_PLOIDY\tTMB_PER_MB\tHRD_SCORE\tMSI_STATUS')
    
    # convert header to lists of lists
    for i in range(len(T)):
        T[i] = T[i].split('\t')

    # update header with clinical information                         
    if sample_info:
        T = update_clinical_sample_header(sample_info, T)
    
    # check if column names include banned columns
    header_columns = get_clinical_fields(T)
    check_column_names(header_columns)
    
    # map columns to positions
    positions = map_columns_to_header(header_columns, T)
                           
    # initialize all columns with empty values beside patient and sample Ids
    data = {}
    for i in S:
        ID = i[0] + ';' + i[1]
        data[ID] = [i[0], i[1]] + ['' for j in range((len(T[0]) - 2))]
    if sample_info:
        for ID in sample_info:
            # add values to clinical fields
            for field in sample_info[ID]:
                data[ID][positions[field]] = sample_info[ID][field]      

    for ID in data:
        assert len(data[ID]) == len(T[0])
    
    # write sample clinical file    
    newfile = open(outputfile, 'w')
    for i in T:
        newfile.write('\t'.join(i) + '\n')
    # sort Ids
    IDs = sorted(list(data.keys()))
    for ID in IDs:
        data[ID] = list(map(lambda x: str(x), data[ID]))
        newfile.write('\t'.join(data[ID]) + '\n')
    newfile.close()






def write_clinical_oncokb(outputfile, mapfile, cancer_code):
    '''
    (str, str, str) -> None
    
    Write clinical file for oncokb-annotator
    
    Parameters
    ----------    
    - outputfile (str): Path to the outputfile
    - mapfile (str): Mapping file (map.csv) that contains paths to maf, segments, expression and fusion files    
    - cancer_code (str): Cancer code from OncoTree
    '''
       
    infile = open(mapfile)    
    content = infile.read().rstrip().split('\n')
    infile.close()
    # create a list of samples
    L = [i.split(',')[1].strip() for i in content]
    
    if L:
        newfile = open(outputfile, 'w')
        header = ['SAMPLE_ID', 'ONCOTREE_CODE']
        newfile.write('\t'.join(header) + '\n')
        for i in L:
            newfile.write('\t'.join([i, cancer_code]) + '\n')
        newfile.close()


   
def concatenate_seg_files(data, outputfile):
    '''
    (str, str) -> None
    
    Concatenates all the segmentation files into a text file outputfile.
    Returns the path to the concatenated file or the empty string if there is no data files    
    
    Parameters
    ----------
    - data (str): Dictionary with data files extracted from the processed mapping file
    - outputfile (str): Path to the concatenated seg file
    '''
        
    # make a list of seg files
    segfiles = []
    # make a parallel list of samples
    samples = []
    
    for donor in data:
        for sample in data[donor]:
            file = data[donor][sample]['segments']
            if file != 'NA' and os.path.isfile(file):
                segfiles.append(file)
                samples.append(sample)
    
    # get the header of the seg file
    if segfiles:
        infile = open(segfiles[0])
        header = infile.readline()
        infile.close()
        
        concatenated_file = outputfile
        
        newfile = open(outputfile, 'w')
        newfile.write(header)
        # concatenate segmentation files
        for i in range(len(segfiles)):
            file = segfiles[i]
            sample = samples[i]
            # get content of seg file
            infile = open(file)
            content = infile.read().split('\n')
            while '' in content:
                content.remove('')
            infile.close()
            # remove header
            content.pop(0)
            # replace ID field with sample name
            for j in range(len(content)):
                content[j] = content[j].split('\t')
                content[j][0] = samples[i]
                content[j] = '\t'.join(content[j])
            # write content of seg file to concatenated file
            newfile.write('\n'.join(content))
        newfile.close()        
    else:
        concatenated_file = ''
        
    return concatenated_file


def write_metadata(outputfile, project_name, data_type, genome):
    '''
    (str, str, str, str) -> None

    Write CNA metadata for a given data type
        
    Parameters
    ----------
    - outputfile (str): Path to the output file
    - project_name (str): Name of project: field project_name in configuration file
    - data_type (str): Type of CNA metadata.
                       Accepted valued: discrete, log2-value, seg, expression, sv, zscore, maf
    - genome (str): Reference genome (hg19 or hg38)
    '''
    
    if data_type == 'discrete':
        stable_id = 'gistic'
        description = 'profile_description: Putative copy-number calls:  Values: -2=homozygous deletion; -1=hemizygous deletion; 0=neutral/no change; 1=gain; 2=high level amplification'
        name = 'Putative copy-number alterations from GISTIC'
        filename = 'data_CNA.txt'
        alteration = 'COPY_NUMBER_ALTERATION'
        data = data_type.upper()
        show_profile = 'true'
    elif data_type == 'log2-value':
        stable_id = 'log2CNA'
        description = 'profile_description: Log2 copy-number values'
        name = 'Log2 copy-number values' 
        filename = 'data_log2CNA.txt'
        alteration = 'COPY_NUMBER_ALTERATION'
        data = data_type.upper()
        show_profile = 'false'
    elif data_type == 'seg':
        filename = 'data_segments.txt'
        description = 'description: Segment data'
        alteration = 'COPY_NUMBER_ALTERATION'
        data = data_type.upper()
    elif data_type == 'sv':
        alteration = 'STRUCTURAL_VARIANT'
        stable_id = 'structural_variants' 
        name = 'Structural Variants'
        filename = 'data_sv.txt'
        description = 'profile_description: Structural variant data'
        data = 'SV'
        show_profile = 'true'
    elif data_type == 'expression':
        alteration = 'MRNA_EXPRESSION'
        data = 'CONTINUOUS'
        stable_id = 'rna_seq_mrna'
        filename = 'data_expression.txt'
        name = 'mRNA expression RNA-Seq'
        description = 'profile_description: Expression levels RNA-Seq'
        show_profile = 'false'
    elif data_type == 'zscore':
        alteration = 'MRNA_EXPRESSION'
        data = 'Z-SCORE'
        stable_id = 'rna_seq_mrna_median_Zscores'
        filename = 'data_expression_zscores.txt'
        description = 'profile_description: Expression levels z-scores'
        name = 'mRNA expression z-scores'
        show_profile = 'true'
    elif data_type == 'maf':
        data = data_type.upper()
        alteration = 'MUTATION_EXTENDED'
        description = 'profile_description: Mutations'
        name = 'Mutations'
        stable_id = 'mutations'
        filename = 'data_mutations_extended.txt'
        show_profile = 'true'

    # collect file text commun to all metadata data types
    L = ['cancer_study_identifier: {0}'.format(project_name),
         'data_filename: {0}'.format(filename),
         'datatype: {0}'.format(data),
         description,
         'genetic_alteration_type: {0}'.format(alteration)]

    # add specific data type text
    if data_type in ['discrete', 'log2-value', 'sv', 'expression', 'zscore', 'maf']:
        L.extend(['stable_id: {0}'.format(stable_id),
                  'show_profile_in_analysis_tab: {0}'.format(show_profile),
                  'profile_name: {0}'.format(name)])
    elif data_type == 'seg':
        L.append('reference_genome_id: {0}'.format(genome))

    newfile = open(outputfile, 'w')
    newfile.write('\n'.join(L))
    newfile.close()        
    

def get_fusfiles_header(fusfiles, samples):
    '''
    (list, list) -> dict
    
    Returns a dictionary with the input data type WT, WG or both  for each sample
    
    Parameters
    ----------
    - fusfiles (list): List of fusion files
    - samples (list): Parallel list of samples
    '''
    
    D = {}
    
    for i in range(len(fusfiles)):
        file = fusfiles[i]
        sample = samples[i]
        # check that fusion file has data
        if check_fusion_data(file):
            D[sample] = []
            infile = open(file)
            header = infile.readline().rstrip().split('\t')
            infile.close()
            for j in range(len(header)):
                if 'WT.' in header[j]:
                    D[sample].append('WT')
                if 'WG.' in header[j]:
                    D[sample].append('WG')
                D[sample].sort()
            
    return D


def extract_fusion(fusion_file):
    '''
    (str) -> list
    
    Returns a list of dictionaies each representing a line of data in fusion_file
    annotated with the column header 
    
    Parameters
    ----------
    - fusion_file (str): Path to the fusion file generated by the mavis workflow
    '''
    
    L = []
    
    infile = open(fusion_file)
    header = infile.readline().rstrip().split('\t')
    for line in infile:
        line = line.rstrip()
        if line:
            line = line.split('\t')
            d = {}
            for i in range(len(header)):
                if 'WT.' in  header[i]:
                    d['WT'] = line[i]
                elif 'WG.' in header[i]:
                    d['WG'] = line[i]
                else:
                    d[header[i]] = line[i]
            L.append(d)        
    infile.close()
    return L



def concatenate_fusion_files(data, outputfile):
    '''
    (str, str) -> None
    
    Concatenates all the fusion files into a text outputfile.
    Returns the path to the concatenated file or the empty string if there is no data files    
    
    Parameters
    ----------
    - data (dict): Dictionary with files extracted from the processed map file
    - outputfile (str): Path to the concatenated seg file
    '''

    # make a list of maf files
    fusionfiles = []
    # make a parallel list of samples
    samples = []
        
    for donor in data:
        for sample in data[donor]:
            file = data[donor][sample]['fusion']
            if file != 'NA' and os.path.isfile(file):
                # check if fusion file has data
                if check_fusion_data(file):
                    fusionfiles.append(file)
                    samples.append(sample)
        
    if fusionfiles:
        concatenated_file = outputfile

        # determine if the headers have WT, WG or both
        header_types = get_fusfiles_header(fusionfiles, samples)
        data_types = []
        for i in header_types:
            data_types.extend(header_types[i])
        data_types = sorted(list(set(data_types)))
    
        header = ['#tracking_id', 'library', 'annotation_id', 'product_id', 'event_type',
                  'gene1', 'gene1_direction', 'gene2', 'gene2_direction', 'gene1_aliases',
                  'gene2_aliases', 'gene_product_type', 'transcript1', 'transcript2',
                  'fusion_splicing_pattern', 'fusion_cdna_coding_start', 'fusion_cdna_coding_end',
                  'fusion_mapped_domains', 'fusion_protein_hgvs', 'annotation_figure',
                  'genes_encompassed', 'break1_chromosome', 'break1_position_start',
                  'break1_position_end', 'break1_orientation', 'exon_last_5prime',
                  'exon_first_3prime', 'break1_strand', 'break2_chromosome',
                  'break2_position_start', 'break2_position_end', 'break2_orientation',
                  'break2_strand', 'protocol', 'tools', 'call_method', 'break1_homologous_seq',
                  'break1_split_reads', 'break2_homologous_seq', 'break2_split_reads',
                  'contig_alignment_score', 'contig_remapped_reads', 'contig_seq',
                  'spanning_reads', 'flanking_pairs', 'linking_split_reads', 'untemplated_seq',
                  'cdna_synon', 'protein_synon', 'supplementary_call', 'net_size',
                  'assumed_untemplated', 'dgv']
    
        # add extra columns indicating the origin of the data
        for i in data_types:
            header.insert(-1, i)
        # add sample to header
        header.insert(0, 'Sample')
    
        # write header to outputfile
        newfile = open(outputfile, 'w')
        newfile.write('\t'.join(header) + '\n')
    
        for i in range(len(fusionfiles)):
            file = fusionfiles[i]
            sample = samples[i]
            # extract data from file
            datafile = extract_fusion(file)
            for d in datafile:
                newline = []
                for i in range(len(header)):
                    if header[i] == 'Sample':
                        newline.append(sample)
                    elif header[i] in d:
                        newline.append(d[header[i]])
                    elif header[i] not in d:
                        newline.append('')
                newfile.write('\t'.join(newline) + '\n')
        newfile.close()        

    else:
        concatenated_file = ''
        
    return concatenated_file


def list_samples_with_expression(data):
    '''
    (dict) -> samples
    
    Returns a list of samples with expression data
    
    Parameters
    ----------
    - data (dict): Dictionary with data files extracted from the processed mapping file
    '''

    # make a list of samples with expression data
    samples = []

    for donor in data:
        for sample in data[donor]:
            file = data[donor][sample]['expression']
            if file != 'NA' and os.path.isfile(file):
                samples.append(sample)

    return samples


def expression_samples_to_file(data, outputfile):
    '''
    (dict, str) -> None
    
    Write list of samples with rna data to outputfile
    
    Parameters
    ----------
    - data (dict): Dictionary with data files extracted from the processed mapping file
    - outputfile (str): path to the outputfile
    '''

    # make a list of samples with expression data
    samples = list_samples_with_expression(data)
    # write samples to file
    newfile = open(outputfile, 'w')
    newfile.write('\n'.join(samples) + '\n')
    newfile.close()     





def extract_expression(gepfile, count):
    '''
    (str, str) -> dict
    
    Returns a dictionary of gene, fpkm or tpm key, value pairs
    
    Parameters
    ----------
    - gepfile (str): Path to the rsem file with fkpm and tpm counts
    - count (str): fpkm or tpm
    '''
    
    # create a dict to store fpkm for each gene
    D = {}
    infile = open(gepfile)
    header = list(map(lambda x: x.lower(), infile.readline().split('\t')))
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split('\t')
            gene = line[0]
            expression = line[header.index(count)]
            assert gene not in D
            D[gene] = expression
    infile.close()
    return D
 
   
 
def extract_expression_isofox(gepfile):
    '''
    (str) -> dict
    
    Returns a dictionary of gene, adjusted tpm key, value pairs
    
    Parameters
    ----------
    - gepfile (str): Path to the gene expression file from isofox
    '''
    
    # create a dict to store fpkm for each gene
    D = {}
    infile = open(gepfile)
    header = infile.readline().rstrip().split(',')
    for line in infile:
        line = line.rstrip()
        if line:
            line = line.split(',')
            gene = line[header.index('GeneId')]
            expression = str(float(line[header.index('AdjTPM')]))
            assert gene not in D
            D[gene] = expression
    infile.close()
    return D
    
 
    

def collect_expression(expression_files, samples, count):
    '''
    (list, list, str) -> dict
    
    Returns a dictionary with fpkm or tpm values for all genes and samples with rna data in directory gepdir
    
    Parameters
    ----------
    - expression_files (list): List of rsem files
    - samples (list): Parallel list of samples. 
    - count (str): fpkm or tpm
    '''
    
    # create a dict to store fpkm or tpm for each gene and sample
    D = {}
    for i in range(len(expression_files)):
        file = expression_files[i]
        sample = samples[i]
        # get the fkpm or tpm for each gene
        expression = extract_expression(file, count)
        assert sample not in D
        D[sample] = expression
    return D



def collect_expression_isofox(expression_files, samples):
    '''
    (list, str) -> dict
    
    Returns a dictionary with fpkm or tpm values for all genes and samples with rna data in directory gepdir
    
    Parameters
    ----------
    - expression_files (list): List of rsem files
    - samples (list): Parallel list of samples. 
    '''
    
    D = {}
    for i in range(len(expression_files)):
        file = expression_files[i]
        sample = samples[i]
        # collect adjusted tpm for each gene
        expression = extract_expression_isofox(file)
        assert sample not in D
        D[sample] = expression
    return D



def write_expression_to_file(D, outputfile):
    '''
    (dict, str) -> None
    
    Write the fpkm or tpm values for all samples and genes in dictionary D to outputfile
    
    Parameters
    ----------
    - D (dict): Dictionary with fpkm or tpm for all samples and gene {sample: {gene: fpkm or tpm}}
    - outputfile (str): Path to the outputfile
    '''
    
    # list all samples
    samples = list(D.keys())
    # list all genes
    genes = []
    for i in samples:
        genes.extend(list(D[i].keys()))
        genes = list(set(genes))
        
    # write fpkm or tpm to file
    newfile = open(outputfile, 'w')    
    header = ['gene_id'] + samples
    newfile.write('\t'.join(header) + '\n')    
    for gene in genes:
        line = [gene]
        for sample in samples:
            line.append(D[sample][gene])
        newfile.write('\t'.join(line) + '\n')    
    newfile.close()
    
   
def concatenate_expression_from_gep_files(data, count, outputfile):
    '''
    (str, str, str) -> None
    
    Write fpkm or tpm for all genes and samples with rna data to outputfile
    Returns the path to the concatenated file or the empty string if there is no data files    
        
    Parameters
    ----------
    - data (dict): Dictionary with files extracted from the processed map file
    - count (str): fpkm or tpm
    - outputfile (str): Path to the outputfile with fpkm or tpm for each sample and gene
    '''
    
    # make a list of maf files
    expressionfiles = []
    # make a parallel list of samples
    samples = []
    
    for donor in data:
        for sample in data[donor]:
            file = data[donor][sample]['expression']
            if file != 'NA' and os.path.isfile(file):
                expressionfiles.append(file)
                samples.append(sample)
    
    if expressionfiles:
        concatenated_file = outputfile
        # extract fpkm or tmp from each file
        expression = collect_expression(expressionfiles, samples, count)   
        # write fpkm to outputfile
        write_expression_to_file(expression, outputfile)
    else:
        concatenated_file = ''

    return concatenated_file
    

def concatenate_expression_from_isofox_files(data, outputfile):
    '''
    (str, str) -> None
    
    Write adjusted tpm for all genes and samples with rna data to outputfile
    Returns the path to the concatenated file or the empty string if there is no data files    
        
    Parameters
    ----------
    - data (dict): Dictionary with files extracted from the processed map file
    - outputfile (str): Path to the outputfile with fpkm or tpm for each sample and gene
    '''
    
    # make a list of maf files
    expressionfiles = []
    # make a parallel list of samples
    samples = []
    
    for donor in data:
        for sample in data[donor]:
            file = data[donor][sample]['expression']
            if file != 'NA' and os.path.isfile(file):
                expressionfiles.append(file)
                samples.append(sample)
    
    if expressionfiles:
        concatenated_file = outputfile
        # extract adjusted tmp from each gene and each sample
        expression = collect_expression_isofox(expressionfiles, samples)
        # write fpkm to outputfile
        write_expression_to_file(expression, outputfile)
    else:
        concatenated_file = ''

    return concatenated_file


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
            break
    infile.close()
    return header
    

def concatenate_maf_files(data, outputfile):
    '''
    (dict, str) -> str
    
    Concatenates all the gzipped maf files into a text file outputfile
    with column Tumor_Sample_Barcode replaced by the sample name.
    Returns the path to the concatenated file or the empty string if there is no data files    
    
    Parameters
    ----------
    - data (dict): Dictionary with files extracted from the processed map file
    - outputfile (str): Path to the concatenated maf file (unzipped)
    '''
     
    # make a list of maf files
    maffiles = []
    # make a parallel list of samples
    samples = []
    
    for donor in data:
        for sample in data[donor]:
            if data[donor][sample]['snv'] != 'NA' and os.path.isfile(data[donor][sample]['snv']):
                maffiles.append(data[donor][sample]['snv'])
                samples.append(sample)
    
    # get the header of the maf file
    if maffiles:
        header = get_maf_header(maffiles[0])
        concatenated_file = outputfile
        newfile = open(outputfile, 'w')
        newfile.write(header + '\n')
        for i in range(len(maffiles)):
            sample = samples[i]
            file = maffiles[i]
            # read the content of the file
            infile = gzip.open(file, 'rt')
            # skip header
            for line in infile:
                if 'version' not in line and 'Hugo_Symbol' not in line:
                    if line.rstrip() != '':
                        line = line.split('\t')
                        # replace column 'Tumor_Sample_Barcode' with sample name
                        j = header.index('Tumor_Sample_Barcode')
                        line[j] = sample
                        # write to outputfile
                        newfile.write('\t'.join(line))
            infile.close() 
    else:
        concatenated_file = ''
    
    return concatenated_file
    
    
                
def removed_filtered_data(outputfile, removed_things, header):
    '''
    (str, str, str)
    
    Write a discription
    
    Parameters
    ----------
    - outputfile (str):
    - removed_things (str):
    - header (str):
    '''
    #open files
    newfile = open(outputfile, 'w')
    # get file header
    #write header to outputfile
    for file_content in header:
        newfile.write(file_content+'\t')
    #write mutations removed to file
    for line in removed_things:
        newfile.write('\n')
        for value in line.values():
            newfile.write(f"{value}\t")

    newfile.close() 



def filter_mutations(maffile, outputfile, depth_filter, alt_freq_filter, gnomAD_AF_filter, keep_variants, removedfile):
    '''
    (str, str, int, float, float, bool) -> (int, int)
    
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
    - keep_variants (bool): Keep variants with missing gnomAD_AF values when Matched_Norm_Sample_Barcode is unmatched if True 
    - removedfile (str): file name for removed file
    '''
    
    # count the number of mutations before and after filtering
    total, kept = 0, 0
    
    # open files
    newfile = open(outputfile, 'w')
    infile = open(maffile)
    
    # get file header
    header = infile.readline().rstrip('\n').split('\t')
    
    # write header to outputfile
    newfile.write('\t'.join(header) + '\n')
    
    # make a list of accepted mutations
    valid_mutations = ['Frame_Shift_Del',  'Frame_Shift_Ins', 'In_Frame_Del',
                          'In_Frame_Ins', 'Missense_Mutation', 'Nonsense_Mutation',
                          'Nonstop_Mutation', 'Silent', 'Splice_Site', 'Translation_Start_Site',
                          'Splice_Region', 'Targeted_Region', "5'Flank"]
    exclude = ['str_contraction', 't_lod_fstar']    

    #ADDED TO SPEED UP SEARCHING
    h_s=header.index('Hugo_Symbol')
    chr_n=header.index('Chromosome')
    s_p=header.index('Start_Position')
    e_p=header.index('End_Position')
    ref_a=header.index('Reference_Allele')
    var_c=header.index('Variant_Classification')
    var_t=header.index('Variant_Type')
    tum_bar=header.index('Tumor_Sample_Barcode')
    tum_al=header.index('Allele')

    filtered_out_list = []

    # apply maf filters to all mutations
    for line in infile:
        # count total mutations
        total += 1
        line = line.rstrip('\n')
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
                    # check that mutation has supporting read counts
                    if line[header.index('t_alt_count')] and line[header.index('t_depth')]:
                        # filter based on ratio t_alt_count / t_depth
                        if int(line[header.index('t_alt_count')]) / int(line[header.index('t_depth')]) >= alt_freq_filter:
                            #start
                            if not (line[var_c] == "5'Flank" and line[h_s] != "TERT"):
                                # filter based on Matched_Norm_Sample_Barcode
                                if line[header.index('Matched_Norm_Sample_Barcode')] == "unmatched":
                                    # check gnomAD_AF. field may be blank, check if value recorded
                                    try:
                                        float(line[header.index('gnomAD_AF')])
                                    except:
                                        # check if variants are kept or not
                                        if keep_variants:
                                            # variants are kept anyway when gnomAD_AF is not defined
                                            newline = line
                                            kept += 1
                                        else:
                                            # no value for gnomAD_AF, do not keep mutation
                                            newline = ''
                                            filtered_out_list.append({"hugo_symbol": line[h_s], "reason": "gnomAD_AF has no value & unmatched barcode", "variant_type": line[var_t], "chr": line[chr_n], "start_pos": line[s_p], "end_pos": line[e_p], "ref_allel": line[ref_a], "tumor_bar": line[tum_bar]})
                                    else:
                                        # compare gnomAD_AF to folder
                                        if float(line[header.index('gnomAD_AF')]) < gnomAD_AF_filter:
                                            newline = line
                                            kept += 1
                                        else:
                                            filtered_out_list.append({"hugo_symbol": line[h_s], "reason": "gnomAD_AF value to larger & unmatched barcode", "variant_type": line[var_t], "chr": line[chr_n], "start_pos": line[s_p], "end_pos": line[e_p], "ref_allel": line[ref_a], "tumor_bar": line[tum_bar]})
                                else:
                                    newline = line
                                    kept += 1
                            else:
                                filtered_out_list.append({"hugo_symbol": line[h_s], "reason": "Hugo symbol and variant class error", "variant_type": line[var_t], "chr": line[chr_n], "start_pos": line[s_p], "end_pos": line[e_p], "ref_allel": line[ref_a], "tumor_bar": line[tum_bar]})
                            #end
                        else:
                            filtered_out_list.append({"hugo_symbol": line[h_s], "reason": "t_alt_count / t_depth", "variant_type": line[var_t], "chr": line[chr_n], "start_pos": line[s_p], "end_pos": line[e_p], "ref_allel": line[ref_a], "tumor_bar": line[tum_bar]})
                    else:
                        # discard mutations without supporting read count
                        newline = ''
                        filtered_out_list.append({"hugo_symbol": line[h_s], "reason": "unsupported read count", "variant_type": line[var_t], "chr": line[chr_n], "start_pos": line[s_p], "end_pos": line[e_p], "ref_allel": line[ref_a], "tumor_bar": line[tum_bar]})
                else:
                    filtered_out_list.append({"hugo_symbol": line[h_s], "reason": "Depth to small", "variant_type": line[var_t], "chr": line[chr_n], "start_pos": line[s_p], "end_pos": line[e_p], "ref_allel": line[ref_a], "tumor_bar": line[tum_bar]})
            else:
                line = line.split('\t')
                removal_reason = 'non valid mutation; ' + line[var_c]
                filtered_out_list.append({"hugo_symbol": line[h_s], "reason": f"{removal_reason}", "variant_type": line[var_t], "chr": line[chr_n], "start_pos": line[s_p], "end_pos": line[e_p], "ref_allel": line[ref_a], "tumor_bar": line[tum_bar]})
            if newline:
                newfile.write('\t'.join(newline) + '\n')
    newfile.close()
    header_list = ["hugo_symbol", "reason", "variant_type", "chr", "start_pos", "end_pos", "ref_allel", "tumor_bar"]
    removed_filtered_data(removedfile, filtered_out_list, header_list)                
    return total, kept                
  

def remove_indels(maffile, outputfile, removedfile):
    '''
    (str, str) -> (int, int)
    
    Writes records from maffile to outputfile with indels removed.
    Returns the total number of mutations before and after filtering
    Precondition: the maf file is unzipped.
    
    Parameters
    ----------
    - maffile (str): Path to the maf file (unzipped)
    - outputfile (str): Path to the output file
    - removedfile (str): Name of removed file
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

    #ADDED TO SPEED UP SEARCHING
    h_s=header.index('Hugo_Symbol')
    chr_n=header.index('Chromosome')
    s_p=header.index('Start_Position')
    e_p=header.index('End_Position')
    ref_a=header.index('Reference_Allele')
    var_c=header.index('Variant_Classification')
    var_t=header.index('Variant_Type')
    tum_bar=header.index('Tumor_Sample_Barcode')


    filtered_out_list = []
    
    for line in infile:
        # count total mutations
        total += 1
        line = line.rstrip()
        if line != '':
            if line[header.index('Variant_Type')] not in ['DEL', 'INS']:
                # record mutations without indels and update counter
                newfile.write('\t'.join(line) + '\n')
                kept += 1
            else:
                line = line.split('\t') 
                removal_reason = 'INDEL; ' + line[var_c]
                filtered_out_list.append({"hugo_symbol": line[h_s], "reason": f"{removal_reason}", "variant_type": line[var_t], "chr": line[chr_n], "start_pos": line[s_p], "end_pos": line[e_p], "ref_allel": line[ref_a], "tumor_bar": line[tum_bar]})

    newfile.close()

    header_list = ["hugo_symbol", "reason", "variant_type", "chr", "start_pos", "end_pos", "ref_allel", "tumor_bar"]
    removed_filtered_data(removedfile, filtered_out_list, header_list)


    newfile.close()                
    return total, kept                


def parse_fusion(fusion_file):
    '''
    (str) -> dict
    
    Returns a dictionary with fusion information
    
    Parameters
    ----------
    - fusion_file (str): Path to the file with fusion information
    '''
    
    D = {}
    
    infile = open(fusion_file)
    header = infile.readline().rstrip().split('\t')
    content = infile.read().rstrip()
    infile.close()
    
    if content:
        content = content.split('\n')
        for i in content:
            i = i.rstrip().split('\t')
            hugo = i[0]
            entrez = i[1]
            center = i[2]
            sample = i[3]
            fusion = i[4]
            dna = i[5]
            rna = i[6]
            method = i[7]
            frame = i[8]
            status = i[9]
        
            # change gene separator in fusion name
            if 'None-' in fusion:
                fusion = fusion.replace('None-', 'None;')
            elif '-None' in fusion:
                fusion = fusion.replace('-None', ';None')
            else:
                if hugo + '-'  in fusion:
                    fusion = fusion.replace(hugo + '-', hugo + ';')
                elif '-' + hugo in fusion:
                    fusion = fusion.replace('-' + hugo, ';' + hugo)
            
            d = {'hugo': hugo,
                'entrez': entrez,
                 'center': center,
                 'fusion': fusion,
                 'dna': dna,
                 'rna': rna,
                 'method': method,
                 'frame': frame,
                 'status': status}
        
            if sample not in D:
                D[sample] = {}
            if fusion not in D[sample]:
                D[sample][fusion] = []
            D[sample][fusion].append(d) 
            
       
    return D   


def list_genes(fusion_file):
    '''
    (str, str) -> str
    
    Returns a list of valid Hugo genes with fusion events
        
    Parameters
    ----------
    - fusion_file (str): Path to the input fusion file
    '''
    
    infile = open(fusion_file)
    header = infile.readline().rstrip().split('\t')
    genes = []
    for line in infile:
        line = line.rstrip()
        if line:
            line = line.split('\t')
            genes.append(line[0])
    genes = list(set(genes))
    
    return genes

    
def convert_fusion_to_sv(fusion_file, sv_file):
    '''
    (str, str) -> str
    
    Convert the fusion file into the SV format 
    
    Parameters
    ----------
    - fusion_file (str): Path to the input fusion file
    - sv_file (str): Path to the output fusion file
    '''

    newfile = open(sv_file, 'w')
    header = ['Sample_Id', 'Site1_Hugo_Symbol', 'Site1_Entrez_Gene_Id',
              'Site2_Hugo_Symbol', 'Site2_Entrez_Gene_Id', 'SV_Status',
              'Center', 'Event_Info', 'DNA_support', 'RNA_support', 
              'Method', 'Site2_Effect_On_Frame', 'Fusion_Status']
    newfile.write('\t'.join(header) + '\n')

    # list all gene names with fusion events
    gene_names = list_genes(fusion_file)
    # parse fusion file
    genes = parse_fusion(fusion_file)
    
    for sample in genes:
        for fusion in genes[sample]:
            Site1_Hugo_Symbol, Site2_Hugo_Symbol, Site1_Entrez_Gene_Id, Site1_Entrez_Gene_Id = '', '', '', ''
            event = fusion.split(';')
            if 'None' in event:
                while 'None' in event:
                    event.remove('None')
                assert len(event) == 1
                Site1_Hugo_Symbol = event[0]
                Site2_Hugo_Symbol = ''
                #Site2_Entrez_Gene_Id = ''
            else:
                if event[0] not in gene_names:
                    Site1_Hugo_Symbol = event[1]
                    Site2_Hugo_Symbol = ''
                    #Site2_Entrez_Gene_Id = ''
                elif event[1] not in gene_names:
                    Site1_Hugo_Symbol = event[0]
                    Site2_Hugo_Symbol = ''
                    #Site2_Entrez_Gene_Id = ''
                else:
                    Site1_Hugo_Symbol, Site2_Hugo_Symbol = event[0], event[1]
                    assert Site1_Hugo_Symbol in gene_names and Site2_Hugo_Symbol in gene_names
                       
            center, dna, rna, method, frame, status = '', '', '', '', '', ''
            for d in genes[sample][fusion]:
                if d['hugo'] == Site1_Hugo_Symbol:
                    Site1_Entrez_Gene_Id = d['entrez']
                    center = d['center']
                    dna = d['dna']
                    rna = d['rna']
                    method = d['method']
                    frame = d['frame']
                    status = d['status']
                elif d['hugo'] == Site2_Hugo_Symbol:
                    Site2_Entrez_Gene_Id = d['entrez']            
                    center = d['center']
                    dna = d['dna']
                    rna = d['rna']
                    method = d['method']
                    frame = d['frame']
                    status = d['status']
                if not Site2_Hugo_Symbol:
                    Site2_Entrez_Gene_Id = ''
            
            assert center and dna and rna and method and frame and status
            L = [sample, Site1_Hugo_Symbol, Site1_Entrez_Gene_Id, Site2_Hugo_Symbol,
                 Site2_Entrez_Gene_Id, 'SOMATIC', center, fusion.replace(';', '-'),
                 dna, rna, method, frame, status]
            newfile.write('\t'.join(L) + '\n')
    
    newfile.close()
    


# def get_sample_info(clinical_samples, outputfile):
#     '''
#     (str, str) -> None
    
#     Write clinical data for samples compatible with IGV
    
#     Parameters
#     ----------
#     - clinical_samples (str): Path to the file with clinical samples information
#     - outputfile (str): Path to the outputfile with content re-formatted for IGV
#     '''
    
#     infile = open(clinical_samples)
#     # find the header, skip commented lines
#     samples = []
#     for line in infile:
#         if not line.startswith('#'):
#             if line.startswith('PATIENT_ID'):
#                 header = line
#             else:
#                 samples.append(line.rstrip())
#     infile.close()
#     while '' in samples:
#         samples.remove('')            
#     for i in range(len(samples)):
#         samples[i] = samples[i].split('\t')
#     # edit header
#     header = header.rstrip().split('\t')
#     header[0] = 'TRACK_ID'
#     header[1] = 'PARTICIPANT_ID'
        
#     newfile = open(outputfile, 'w')
#     newfile.write('\t'.join(header) + '\n')
#     for i in samples:
#         newfile.write('\t'.join([i[1], i[0]]) + '\n')
#     newfile.close()
    

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



def extract_config_params(config):
    '''
    (configparser.ConfigParser) -> dict
    
    Returns a dictionary with parameters for each section of the config file
        
    Parameters
    ----------
    - config (configparser.ConfigParser): Config file parsed with configparser
    '''
    
    D = {}
    
    # extract resource files
    D['Resources'] = {}
    for i in list(config['Resources'].keys()):
        D['Resources'][i] = config['Resources'][i]
    
    # extract options
    D['Options'] = {}
    for i in list(config['Options'].keys()):
        if i == 'keep_variants':
            D['Options'][i] = config['Options'].getboolean(i)
        elif i == 'gamma':
            D['Options'][i] = config['Options'].getint(i)
        else:
            D['Options'][i] = config['Options'][i]
            
    # extract  parameters
    D['Parameters'] = {}
    for i in list(config['Parameters'].keys()):
        if '.' in config['Parameters'][i]:
            D['Parameters'][i] = config['Parameters'].getfloat(i)
        else:
            D['Parameters'][i] = config['Parameters'].getint(i)
    
    # extract filtering parameters
    D['Filters'] = {}
    for i in list(config['Filters'].keys()):
        if i in ['depth_filter', 'alt_freq_filter', 'gnomad_af_filter']:
            D['Filters'][i] = config['Filters'].getfloat(i)
        elif i in ['tglpipe', 'filter_variants', 'filter_indels']:
            D['Filters'][i] = config['Filters'].getboolean(i)
    
    # extract workflows
    D['Workflows'] = {}
    for i in list(config['Workflows'].keys()):
        D['Workflows'][i] = config['Workflows'][i]

    return D



def check_configuration(params):
    '''
    (dict) -> None
    
    Check the content of the config file and raise a ValueError if required information is missing or not as expected
        
    Parameters
    ----------
    - params (dict): Dictionary with parameters extracted from the config 
    '''
    
    # raise an error if a section is omitted
    missing_sections = [i for i in ['Resources', 'Options', 'Parameters', 'Filters', 'Workflows'] if i not in params.keys()] 
    if missing_sections:
        raise ValueError('ERROR. Missing sections {0} from config'.format(', '.join(missing_sections)))
            
    # check paths from resources
    expected_resources = ['token', 'enscon_hg38', 'enscon_hg19', 'entcon', 'genebed_hg38', 'genebed_hg19', 'genelist', 'oncolist']
    missing_resources = [i for i in expected_resources if i not in params['Resources'].keys()]
    if missing_resources:
        raise ValueError('ERROR. Missing resources: {0}'.format(', '.join(missing_resources)))
    invalid_resource_files = []
    for i in expected_resources:
        # genelist is optional and may be None
        if i == 'genelist' and params['Resources'][i] and os.path.isfile(params['Resources'][i]) == False:
            invalid_resource_files.append(i)                
        elif params['Resources'][i] is None or os.path.isfile(params['Resources'][i]) == False:
            invalid_resource_files.append(i)
    if invalid_resource_files:
        raise ValueError('ERROR. Provide valid path for {0}'.format(', '.join(invalid_resource_files)))
    
    # check options
    missing_options = [i for i in ['mapfile', 'outdir', 'study', 'center', 'cancer_code', 'keep_variants', 'count', 'gamma'] if params['Options'][i] is None]
    if missing_options:
        raise ValueError('ERROR. Provide values for {0} in the config'.format(', '.join(missing_options)))
    # check map file
    if os.path.isfile(params['Options']['mapfile']) == False:
        raise ValueError('ERROR: Provide valid path to mapfile in config')
    # check that bbolean filter parameter is provided
    if params['Options']['keep_variants'] not in [True, False]:
        raise ValueError('ERROR. {0} is not a boolean. Use true or false'.format(params['Options']['keep_variants']))
            
    # check parameters
    expected_parameters = ['gain', 'amplification', 'heterozygous_deletion', 'homozygous_deletion', 'minfusionreads']
    missing_parameters = [i for i in expected_parameters if i not in params['Parameters'].keys()]
    if missing_parameters:
        raise ValueError('ERROR. Missing Parameters: {0}'.format(', '.join(missing_parameters)))
    # check value types
    for i in ['gain', 'amplification', 'heterozygous_deletion', 'homozygous_deletion', 'minfusionreads']:
        try:
            float(params['Parameters'][i])
        except:
            raise ValueError('ERROR. {0} is not a number'.format(i))
      
    # check filters
    expected_filters = ['tglpipe', 'filter_variants', 'depth_filter', 'alt_freq_filter', 'gnomad_af_filter', 'filter_indels']
    missing_filters = [i for i in expected_filters if i not in params['Filters'].keys()]
    if missing_filters:
        raise ValueError('ERROR. Missing Filters: {0}'.format(', '.join(missing_filters)))
    # check value types
    for i in ['tglpipe', 'filter_variants', 'filter_indels']:
        try:
            bool(params['Filters'][i])
        except:
            raise ValueError('ERROR. {0} is not a boolean'.format(i))
    for i in ['depth_filter', 'alt_freq_filter', 'gnomAD_AF_filter']:
        try:
            float(params['Filters'][i])
        except:
            raise ValueError('ERROR. {0} is not a number'.format(i))

    # check data types
    if params['Workflows']['snp'] not in ['vep', 'pave']:
        raise ValueError('ERROR. Expecting vep or pave for snp')
    if params['Workflows']['expression'] not in ['rsem', 'isofox']:
        raise ValueError('ERROR. Expecting rsem or isofox for expression')
    if params['Workflows']['cna'] not in ['purple', 'sequenza']:
        raise ValueError('ERROR. Expecting purple or sequenza for cna')
    if params['Workflows']['fusion'] not in ['mavis', 'linx']:
        raise ValueError('ERROR. Expecting mavis or linx for fusion')
        
    


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
   
    
    
# def copy_segmentation_data(cbiodir, suppdir):
#     '''
#     (str, str) -> None
    
#     Copy file processed segmentation file data_segments.txt from the cbioportal
#     import folder cbiodir to the supplementary folder suppdir and rename it to data_segments.seg
    
#     Parameters
#     ---------
#     - cbiodir (str): Path to the cbioportal import folder
#     - suupdir (str): Path the supplementary folder 
#     '''
    
#     if os.path.isfile(os.path.join(cbiodir, 'data_segments.txt')):
#         shutil.copy(os.path.join(cbiodir, 'data_segments.txt'), os.path.join(suppdir, 'data_segments.seg'))


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


    
def check_input_mafs(mapfile):
    '''
    (str) -> None
    
    Exits if the input maf files have different headers
        
    Parameters
    ----------
    - mapfile (str): Mapping file (map.csv) that contains paths to maf, segment, expression and fusion files    
    '''

    # make a list of input maf files 
    mafs = extract_files_from_map(mapfile, 'snv')
    
    # make a list of file headers
    headers = []
    if mafs:
        for file in mafs:
            infile = gzip.open(file, 'rt')
            for line in infile:
                if 'Hugo_Symbol' in line:
                    headers.append(line.rstrip())
                    break
            infile.close()
    
    # check if multiple headers
    if len(list(set(headers))) > 1:
        sys.exit('Input MAF files have different headers')


def extract_samples_case_file(case_file):
    '''
    (str) -> list
    
    Returns a list of case samples extracted from the case_file of a previous import folder
        
    - case_file (str): Path to the case file located in folder cbioportal_import_data/case_lists/
                       Possible files are: cases_sequenced.txt, cases_rna_seq_mrna.txt,
                       cases_cna.txt, cases_cnaseq.txt, cases_3way_complete.txt, cases_sv.txt
    '''
    
    infile = open(case_file)
    content = infile.read().rstrip().split('\n')
    infile.close()
    samples = content[-1].split(':')[1]
    samples = list(map(lambda x: x.strip(), samples.split('\t')))
    
    return samples                    



def check_fusion_data(fusfile):
    '''
    (str) -> bool
    
    Returns True if fusfile contains SV variants and False otherwise
    
    Parameters
    ----------
    - fusfile (str): Path to the concatenated fusion file
    '''
    
    infile = open(fusfile)
    header = infile.readline()
    data = infile.read().strip()
    infile.close()
    
    if data:
        return True
    else:
        return False


def copy_resource(file, destination_dir):
    '''
    (str, str) -> None
    
    Copy file in the destination directory
    
    Parameters
    ----------
    - file (str): Path to the file to be copied
    - destination_dir (str): Path to the destination where file is copied
    '''
    
    filename = os.path.basename(file)
    copied_file = os.path.join(destination_dir, filename)
    shutil.copyfile(file, copied_file)    
    




def get_paths(paths):
    '''
    (str) -> list
    
    Returns a list of paths to files or directories contained in file paths
    
    Parameters
    ----------
    - paths (str): Path to the file containing file or folder paths
    '''

    infile = open(paths)
    L = infile.read().rstrip().split('\n')
    infile.close()

    return L


def check_configs_import_folders(config_files, import_folders):
    '''
    (list, list) -> None
    
    Exits if config files and import folders are not conform
    
    Parameters
    ----------
    - config_files (list): List of paths to config files of import folders to merge
    - import_folders (list): List of paths to import folders to merge
    '''

    if len(config_files) == 0:
        sys.exit('No specified config files. Use -configs or -configPaths')
    if len(import_folders) == 0:
        sys.exit('No specified import folders. Use -importFolders or -importPaths')

    # check that numbers of import folders and config files match        
    if len(config_files) > len(import_folders):
        missing = len(config_files) - len(import_folders) 
        sys.exit('Missing {0} import folders'.format(missing))
    elif len(config_files) < len(import_folders):
        missing = len(import_folders) - len(config_files) 
        sys.exit('Missing {0} config files'.format(missing))

    # check that config files are valid
    if all(map(lambda x: os.path.isfile(x), config_files)) == False:
        invalid = list(map(lambda x: os.path.isfile(x), config_files)).count(False)
        sys.exit('{0} config files have invalid paths'.format(invalid))    
    # check that import folders are valid
    if all(map(lambda x: os.path.isdir(x), import_folders)) == False:
        invalid = list(map(lambda x: os.path.isdir(x), import_folders)).count(False)
        sys.exit('{0} import folders have invalid paths'.format(invalid))    




def collect_config_parameters(config_files):
    '''
    (list) -> dict
     
    Returns a dictionary with lists of parameter values across config files
    
    Parameters
    ----------
    - config_files (list): List of config files of import folders to merge
    '''
    
    # collect all parameters
    D = {}
    for file in config_files:
        # parse config file
        config = configparser.ConfigParser(allow_no_value=True)
        config.read(file)
        # extract variables from config
        params = extract_config_params(config)
        # check config
        check_configuration(params)
        for section in params:
            if section not in D:
                D[section] = {}
            for k in params[section]:
                if k not in D[section]:
                    D[section][k] = []
                D[section][k].append(params[section][k])    
                D[section][k] = list(set(D[section][k]))             
    
    return D
    



def check_config_parameters(D):
    '''
    (dict) -> (bool, list)    
    
    Returns True if parameters are the same across the config files and False otherwise,
    and also returns a list of parameters with differences
    
    Parameters
    ----------
    - D (dict): Dictionaries with parameters extracted from config files
    '''

    # set parameter to evaluate that all parameters are the same across configs
    same_parameters = True
    
    # report parameter with differences     
    
    differences = []
    
    for section in D:
        # skip Workflows. can only merge identical data types
        if section != 'Workflows':
            for k in D[section]:
                # skip genome. can only merge data with the same reference
                if k != 'genome':
                    if len(D[section][k]) > 1:
                        # update boolean, configs have different parameters
                        same_parameters = False
                        differences.append(k)
    
    return same_parameters, differences




def check_genome(D, genome):
    '''
    (dict, str) -> None
    
    Exits if the reference genome does not match the genome from config files
    or if multiple genome references are specified in the config files.
    
    Parameters
    ----------
    - D (dict): Dictionary with parameters extracted from the config files of the import folders to be merged
    - genome (str): Reference genome
    '''
    
    
    if len(D['Options']['genome']) > 1:
        sys.exit('The config files have different genomes: {0}'.format(';'.join(D['Options']['genome'])))
    elif len(D['Options']['genome']) == 0:
        sys.exit('The config files do not have genomes')
    elif D['Options']['genome'][0] != genome:
        sys.exit('''Attempting to merge data with different genome references.
                 ({0} from command, {1} from configs)'''.format(genome, D['Options']['genome'][0]))


def check_datatypes(D):
    '''
    (dict) -> None
    
    Exits if the data types specified are different
    
    Parameters
    ----------
    - D (dict): Dictionary with parameters extracted from the config files of the import folders to be merged
    '''
    
    L = []
    for k in D['Workflows']:
        if len(D['Workflows'][k]) > 1:
            L.append('{0}: {1}'.format(k, ','.join(D['Workflows'][k])))
    if len(D['Options']['count']) > 1:
        L.append('count: {0}'.format(','.join(D['Workflows']['count'])))
    if len(D['Options']['gamma']) > 1:
        L.append('gamma: {0}'.format(','.join(D['Workflows']['gamma'])))
    if L:
        sys.exit('Attempting to merge different data types: {0}'.format(';'.join(L)))
    

def check_metadata(D, metadata, metadata_type):
    '''
    (dict, str, str) -> None
    
    Print a warning if the metadata specified in the command does not match the
    metadata from config files if multiple values are specified in the config files.
    
    Parameters
    ----------
    - D (dict): Dictionary with parameters extracted from the config files of the import folders to be merged
    - metadata (str): Cancer code, project, study name or description 
    - metadata_type (str): The type of metadata to check (cancer_code, project_name, study, description)
    '''
    
    name = ' '.join(metadata_type.split('_'))
    vals = ';'.join(D[metadata_type])
    
    if len(D['Options'][metadata_type]) > 1:
        print('WARNING: The config files have different {0}s: {1}'.format(name, vals))
    elif len(D['Options'][metadata_type]) == 0:
        print('WARNING: The config files do not have {0}s'.format(name))
    elif D['Options'][metadata_type][0] != metadata:
        print('''WARNING: Attempting to merge data with different {0}s. 
                 ({1} from command, {2} from configs)'''.format(name, metadata, vals))
       

def group_files(import_folders, case_list = False):
    '''
    (list, bool) -> dict
    
    Returns a dictionary with list of file paths for each expected file name
    in the import folder or in the case_lists folder if case_list is True
        
    Parameters
    ----------
    - import_folders (list): List of import folders to merge
    - case_list (bool): Get the files in the case_lists folder is True 
                        and in the import folder if False
    '''

    D = {}

    for i in import_folders:
        if case_list:
            # get the files in case_lists folder
            folder = os.path.join(i, 'case_lists')
        else:
            folder = i
        files = [os.path.join(folder, j) for j in os.listdir(folder)]
        for file in files:
            filename = os.path.basename(file)
            if filename in D:
                D[filename].append(file)
            else:
                D[filename] = [file]
    
    return D  






def merge_case_files(case_dir, import_folders, project_name, excluded_samples):
    '''
    (str, list, str, list, list)
    
    Extract samples from each case file in case_lists directory of the import folders
    and write a new case file in the case_dir merging all samples from previous import folders
    
    Parameters
    ==========
    - case_dir (str): Path to the case_lists directory where the merged case file is written
    - import_folders (list): List of paths to import folders to merge
    - project_name (str): Name of the project
    - excluded_samples (list): List of samples to remove 
    '''
    
    # get the list of files to merge in the case_lists folder
    case_lists_files = group_files(import_folders, True)

    # collect samples for each type of case file and write a new case file to case_dir
    for filename in case_lists_files:
        # extract samples from each file corresponding to the file name
        samples = []
        for case_file in case_lists_files[filename]:
            samples.extend(extract_samples_case_file(case_file))

        if filename == 'cases_sequenced.txt':
            name = 'Samples profiled for mutations'
            description = 'This is this case list that contains all samples that are profiled for mutations.'
            stable_id = '{0}_sequenced'.format(project_name)
        elif filename == 'cases_rna_seq_mrna.txt':
            name = 'Samples profiled for rnaseq'
            description = 'This is this case list that contains all samples that are profiled for rnaseq.'
            stable_id = '{0}_rna_seq_mrna'.format(project_name)    
        elif filename == 'cases_cna.txt':
            name = 'Samples profiled for cnas'
            description = 'This is this case list that contains all samples that are profiled for cnas.'
            stable_id = '{0}_cna'.format(project_name)
        elif filename == 'cases_cnaseq.txt':
            name = 'Samples profiled for cnas and sequencing'
            description = 'This is this case list that contains all samples that are profiled for mutations and cnas.'
            stable_id = '{0}_cnaseq'.format(project_name)
        elif filename == 'cases_3way_complete.txt':
            name = 'Samples profiled for all of mutation, cnas, and rnaseq'
            description = 'This is this case list that contains all samples that are profiled for mutations, cnas, and rnaseq.'
            stable_id = '{0}_3way_complete'.format(project_name)
        elif filename == 'cases_sv.txt':
            name = 'Samples profiled for structural variants'
            description = 'This is this case list that contains all samples that are profiled for structural variants.'
            stable_id = '{0}_sv'.format(project_name)
   
    
        # remove samples to exclude
        if excluded_samples:
            for i in excluded_samples:
                if i in samples:
                    samples.remove(i)
          
        if samples:
            # write outputfile if samples exist
            outputfile = os.path.join(case_dir, filename)
            newfile = open(outputfile, 'w')
            L = ['cancer_study_identifier: {0}'.format(project_name),
                 'stable_id: {0}'.format(stable_id),
                 'case_list_name: {0}'.format(name),
                 'case_list_description: {0}'.format(description),
                 'case_list_ids: {0}'.format('\t'.join(samples))]
            newfile.write('\n'.join(L))
            newfile.close()



def merge_meta_files(cbiodir, import_folders, project_name, study, description, genome, cancer_code):
    '''
    (str, list, str)
    
    Extract samples from each case file in case_lists directory of the import folders
    and write a new case file in the case_dir merging all samples from previous import folders
    
    Parameters
    ----------
    - cbiodir (str): Path to the cbioportal_import_data folder 
    - import_folders (list): List of paths to import folders to merge
    - project_name (str): Name of the project
    - study (str): Name of the study
    - description (str): Study description
    - genome (str): Reference genome    
    - cancer_code (str): Oncotree cancer code
    '''
    
    # get the list of files to merge in the import folder 
    cbiofiles = group_files(import_folders, False)
    
    # # keep only the metadata files
    # to_remove = [i for i in cbiofiles if not i.startswith('meta_')]
    # for i in to_remove:
    #     del cbiofiles[i]
    
    if 'meta_study.txt' in cbiofiles:
        write_meta_study(os.path.join(cbiodir, 'meta_study.txt'), study, project_name, description, genome, cancer_code)
    if 'meta_clinical_samples.txt' in cbiofiles:
        write_meta_clinical(cbiodir, project_name, 'sample')
    if 'meta_clinical_patients.txt' in cbiofiles:
        write_meta_clinical(cbiodir, project_name, 'patient')
    
    metafiles =  {'meta_mutations_extended.txt': 'maf',
                  'meta_CNA.txt': 'discrete',
                  'meta_log2CNA.txt': 'log2-value',
                  'meta_segments.txt': 'seg',
                  'meta_expression.txt': 'expression',
                  'meta_expression_zscores.txt': 'zscore',
                  'meta_sv.txt': 'sv'}
    
    for i in metafiles:
        if i in cbiofiles:
            write_metadata(os.path.join(cbiodir, i), project_name, metafiles[i], genome)
    
    print('wrote metadata files')
    
    
    
def remove_merged_metafiles(cbiodir):
    '''
    (str) -> None

    Remove merged metadata files if corresponding data files do not exist
    because samples have been excluded
    
    Parameters
    ----------
    - cbiodir (str): Path to the cbioportal_import_data folder 
    '''    
    
    # make a list of metada files in cbiodir
    # exclude meta_study.txt
    metafiles = [os.path.join(cbiodir, i) for i in os.listdir(cbiodir) if 'meta_' in i and 'study' not in i]
    
    if metafiles:
        for i in metafiles:
            filename = os.path.basename(i)
            datafilename = filename.replace('meta', 'data')
            datafile = os.path.join(cbiodir, datafilename)
            if os.path.isfile(datafile) == False:
                os.remove(i)
                print('Removed metadafile {0} because {1} does not exist'.format(filename, datafilename)) 
        
        
    
def get_simple_header(infile):
    '''
    (_io.TextIOWrapper) -> str
    
    Returns the header of the open file
    
    Parameters
    ----------
    infile : (_io.TextIOWrapper): File open for reading
    '''    
    
    return infile.readline().rstrip()
    

def get_multiline_header(infile):
    '''
    (_io.TextIOWrapper) -> list
    
    Returns a list with all the lines of the file header
    
    Parameters
    ----------
    infile : (_io.TextIOWrapper): File open for reading
    '''

    header = []
    for line in infile:
        if line.startswith('#'):
            header.append(line)
        else:
            break
    header.append(line)
     
    return header
    
    
def check_file_headers(cbiofiles):
    '''
    (dict) -> None
    
    Exits if data file headers are different (ie, cannot be merged)
    
    Parameters
    ----------
    - cbiofiles (dict): Dictionary with list of file paths for each expected
                        file found in the import folders to merge
    '''    
    
    # check header for each data file
    D = {}
        
    for i in cbiofiles:
        if i in ['data_mutations_extended.txt', 'data_segments.txt', 'data_sv.txt']:
            D[i] = True
            # get the header of the first file
            infile = open(cbiofiles[i][0])
            header = get_simple_header(infile)
            infile.close()
            for file in cbiofiles[i]:
                infile = open(file)
                h = get_simple_header(infile)
                infile.close()
                if header != h:
                    D[i] = False
        elif i in ['data_clinical_patients.txt', 'data_clinical_samples.txt']:
            D[i] = True
            infile = open(cbiofiles[i][0])
            header = get_multiline_header(infile)
            infile.close()
            for file in cbiofiles[i]:
                infile = open(file)
                h = get_multiline_header(infile)
                infile.close()
                if header != h:
                    D[i] = False        
    
    different_headers = [i for i in D if D[i] == False]
    
    if different_headers:
        sys.exit('WARNING: Some data files have different headers and cannopt be merged: {0}'.format(';'.join(different_headers)))    
    
    
def merge_expression_CNA_files(cbiodir, cbiofiles, excluded_samples):
    '''
    (str, dict, list) -> None
    
    Merge expression and CNA data files into the cbiodir directory
        
    Parameters
    ----------
    - cbiodir (str): Path to the cbioportal_import_data folder 
    - cbiofiles (dict): Dictionary with list of file paths for each expected
                        file found in the import folders to merge
    - excluded_samples (list): List of samples to remove
    '''
    
    for i in cbiofiles:
        if i in ['data_expression.txt', 'data_expression_zscores.txt',
                 'data_log2CNA.txt', 'data_CNA.txt']:
            D = {}
            for file in cbiofiles[i]:
                infile = open(file)
                header = infile.readline().rstrip().split('\t')
                for line in infile:
                    line = line.rstrip()
                    if line:
                        line = line.split('\t')
                        gene = line[0]
                        if gene not in D:
                            D[gene] = {}
                        for j in range(1, len(line)):
                            sample = header[j]
                            expression = line[j]
                            if sample not in excluded_samples:
                                D[gene][sample] = expression
                infile.close()            
                
            if D:
                h = ['Hugo_Symbol']
                samples = []
                for gene in D:
                    samples.extend(list(D[gene].keys()))
                    samples = list(set(samples))
                samples.sort()
                h.extend(samples)
                genes = sorted(D.keys())
                newfile = open(os.path.join(cbiodir, i), 'w')
                newfile.write('\t'.join(h) + '\n')
                for gene in genes:
                    L = [gene]
                    for sample in samples:
                        assert sample in D[gene]
                        L.append(D[gene][sample])
                    newfile.write('\t'.join(L) + '\n')
                newfile.close()
        
    

def merge_data_files(cbiodir, cbiofiles, excluded_samples):
    '''
    (str, dict, list) -> None
    
    Merge data files into the cbiodir directory and exclude samples if needed
        
    Parameters
    ----------
    - cbiodir (str): Path to the cbioportal_import_data folder 
    - cbiofiles (dict): Dictionary with list of file paths for each expected
                        file found in the import folders to merge
    - excluded_samples (list): List of samples to remove
    '''

    for i in cbiofiles:
        if i in ['data_mutations_extended.txt', 'data_segments.txt', 'data_sv.txt']:
            data = []
            for file in cbiofiles[i]:
                infile = open(file)
                header = infile.readline().split('\t')
                for line in infile:
                    if excluded_samples:
                        if line.strip():
                            line = line.split('\t')
                            if i  == 'data_mutations_extended.txt':
                                tumor_sample = line[header.index('Tumor_Sample_Barcode')]
                                matched_normal = line[header.index('Matched_Norm_Sample_Barcode')]
                                if tumor_sample not in excluded_samples and matched_normal not in excluded_samples:
                                    data.append('\t'.join(line))
                            elif i in ['data_segments.txt', 'data_sv.txt']:
                                sample = line[0]
                                if sample not in excluded_samples:
                                    data.append('\t'.join(line))
                    else:
                        data.append(line)
                    
            newfile = open(os.path.join(cbiodir, i), 'w')
            newfile.write('\t'.join(header))
            newfile.write(''.join(data))
            newfile.close()
        
        
        
def convert_samples_to_donors(samples):
    '''
    (list) -> list
    
    Returns a list of donor Ids  corresponding to the list of samples
    
    Parameters
    ----------
    samples (list): List of sample Ids
    '''        
        
    donors = ['_'.join(i.split('_')[:2]) for i in samples]
    donors = list(set(donors))
    
    return donors
    
        


def merge_clinical_files(cbiodir, cbiofiles, excluded_samples):
    '''
    (str, dict, list) -> None
    
    Merge clinical data files into the cbiodir directory
        
    Parameters
    ----------
    - cbiodir (str): Path to the cbioportal_import_data folder 
    - cbiofiles (dict): Dictionary with list of file paths for each expected
                        file found in the import folders to merge
    - excluded_samples (list): List of samples to remove
    '''

    if excluded_samples:
        # get the donor ids
        donors = convert_samples_to_donors(excluded_samples)

    for i in cbiofiles:
        if i in ['data_clinical_patients.txt', 'data_clinical_samples.txt']:
            data = []
            for file in cbiofiles[i]:
                infile = open(file)
                # collect multi line header
                header = get_multiline_header(infile)
                # keep trailing white space in lines (multiple tabs)
                for line in infile:
                    if excluded_samples:
                        line = line.split('\t')
                        # skip records if samples or donors are xcluded
                        if (line[0] in donors or line[0] in excluded_samples) \
                            or (line[1] in donors or line[1] in excluded_samples):
                                continue
                        else:
                            data.append('\t'.join(line))
                    else:
                        data.append(line)
                infile.close()
                          
            newfile = open(os.path.join(cbiodir, i), 'w')
            newfile.write(''.join(header))
            newfile.write(''.join(data))
            newfile.close()




def check_mapfile(params, mapfile):
    '''
    (dict, str) -> int
    
    Check the expected number of columns in the map file
    and raise a ValueError if the number of columns differ from expected based
    on the source of data. Returns the number of columns
    
    Parameters
    ----------
    - params (dict): Dictionary with options extracted from the config
    - mapfile (str): Path to the mapping file
    '''
    
    L = []
    infile = open(mapfile)
    for line in infile:
        line = line.rstrip()
        if line:
            line = line.split(',')
            L.append(len(line))
    infile.close()
    L = list(set(L))    
       
    if len(L) == 0:
        line_length = 0
        raise ValueError('No data in {0}'.format(mapfile))
    elif len(L) > 1:
        line_length = -1
        raise ValueError('Multiple line length detected in {0}. Replace missing values with NA'.format(mapfile))
    else:
        line_length = L[0]
        if params['Workfows']['cna'] == 'purple':
            # expecting donor,sample,maf,purity,somatic,expression,fusion
            # with purity and somatic being purple outputs used to compute the segments
            if line_length != 7:
                raise ValueError('Expecting 7 columns: donor,sample,maf,purity,somatic,expression,fusion in {0} with {1}'.format(mapfile, params['Workfows']['cna']))
        else:
            if line_length != 6:
                raise ValueError('Expecting 6 columns: donor,sample,maf,segments,expression,fusion in {0}'.format(mapfile))
    
    return line_length
    
    
    
def parse_mapfile(params, mapfile):
    '''
    (dict, str) -> dict
    
    Returns a dictionary with data files from the mapping file for each donor and sample
    
    Parameters
    ----------
    - params (dict): Dictionary with options extracted from the config
    - mapfile (str): Path to the mapping file
    '''
    
    infile = open(mapfile)
    content = infile.read().split('\n')
    content = list(map(lambda x: x.strip().split(','), content))
    infile.close()
    L = list(set(list(map(lambda x: len(x), content))))
    assert len(L) == 1
    cols = L[0]    
    
    if params['Workflows']['cna'] != 'purple':
        assert cols == 7
    else:
        assert cols == 6
        
    D = {}
    
    for i in content:
        donor = i[0]
        sample = i[1]
        maf = i[2]
        if params['Workflows']['cna'] != 'purple':
            segments = i[3]
            expression = i[4]
            fusion = i[5]
        else:
            segments = i[3:5]
            expression = i[5]
            fusion = i[6]
        if donor in D:
            D[donor] = {}
        assert sample not in D[donor]
        D[donor][sample] = {'snv': maf,
                            'segments': segments,
                            'expression': expression,
                            'fusion': fusion}
    
    return D



def parse_processed_mapfile(processed_mapfile):
    '''
    (str) -> dict
    
    Returns a dictionary with data files from the mapping file for each donor and sample
    
    Parameters
    ----------
    - processed_mapfile (str): Path to the processed mapping file with paths to data files
    '''
        
    D = {}
    
    infile = open(processed_mapfile)
    for line in infile:
        line = line.rstrip()
        if line:
            line = line.split(',')
            donor = line[0]
            sample = line[1]
            maf = line[2]
            segments = line[3]
            expression = line[4]
            fusion = line[5]
                    
            if donor not in D:
                D[donor] = {}
            assert sample not in D[donor]
            D[donor][sample] = {'snv': maf,
                                'segments': segments,
                                'expression': expression,
                                'fusion': fusion}
    infile.close()
   
    return D


def get_sequenza_segfile(sampledir, sequenza_zip, gamma):
    '''
    (str, str, int)
    
    Returns the sequenza segmentation file corresponding to the given gamma
    after unziping sequenza files to the sample directory
    
    Parameters
    ----------
    - sampledir (str): Directory where sequenza output is unzipped
    - sequenza_zip (str): Path to the sequenza output zip file
    - gamma (int): Particular gamma solution for segmnents
    '''
    
    segfile = 'NA'
    
    # unzip the sequenza file
    if os.path.isfile(sequenza_zip):
        with zipfile.ZipFile(sequenza_zip, 'r') as zipref:
            zipref.extractall(sampledir)
        # get the seg file of interest
        segfile = glob.glob(sampledir +'/gammas/{0}/*.seg'.format(gamma))
        assert len(segfile) == 1
        segfile = segfile[0]
        assert os.path.isfile(segfile)
        
    return segfile           
    
    
def preprocess_mapfile(params, mapfile, gamma, outdir):
    '''
    (dict, str, int, str) -> str
    
    Parse the mapfile and 
    
    Parameters
    ----------
    - params (dict): Dictionary with options extracted from the config
    - mapfile (str): Path to the mapping file
    - gamma (str): Gamma value of the sequenza workflow
    - outdir (str): Path to the output directory 
    '''

    processed_map = os.path.join(outdir, 'preprocessed_map.csv')
    newfile = open(processed_map, 'w')

    data = parse_mapfile(params, mapfile)

    # link the files to the data folders

    for donor in data:
        for sample in data[donor]:
            L = [donor, sample, 'NA', 'NA', 'NA', 'NA']
            segdir = os.path.join(outdir, 'segdir')   
            if params['Workflows']['cna'] == 'sequenza':
                # extract the seg file corresponding to gamma
                # create a sequenza directory
                sequenzadir =  os.path.join(segdir, 'sequenza')
                donordir = os.path.join(sequenzadir, donor)
                sampledir = os.path.join(donordir, sample)
                os.makedirs(sampledir, exist_ok=True)
                # get the sequenza seg file
                segfile = get_sequenza_segfile(sampledir, data[donor][sample]['segments'], gamma)
                L[3] = segfile       
            elif params['Workflows']['cna'] == 'purple':
                # create a purple directory
                # generate segmentation file for purple
                purpledir =  os.path.join(segdir, 'purple')
                donordir = os.path.join(purpledir, donor)
                sampledir = os.path.join(donordir, sample)
                os.makedirs(sampledir, exist_ok=True)
                #get the cnv and purity files
                cnv_file = data[donor][sample]['segments'][0]
                purity_file = data[donor][sample]['segments'][1]
                # generate segmentation file from purple
                outputfile = os.path.join(sampledir, '{0}.purple.seg'.format(sample))
                generate_purple_segmentation(cnv_file, purity_file, sample, outputfile)
                L[3] = outputfile
                
            # replace files for other data type
            L[2] = data[donor][sample]['snv']
            L[4] = data[donor][sample]['expression']
            L[5] = data[donor][sample]['fusion']
            newfile.write(','.join(L) + '\n')

    newfile.close()
    return processed_map


def make_import_folder(args):
    '''
    (list) -> None
    
    Generate folder cbioportal_import_folder with metadata and processed data files
    to upload to cBioPortal
        
    Parameters
    ----------
    - config (str): Path to the config file
    - clinical (str): Path to the sample clinical file
    '''
    
    if os.path.isfile(args.config) == False:
        sys.exit('Provide a valid config file')
    if args.clinical:
        if os.path.isfile(args.clinical) == False:
            sys.exit('Provide a valid clinical file')
        
    # parse config file
    config = configparser.ConfigParser(allow_no_value=True)
    config.read(args.config)
    # extract variables from config
    params = extract_config_params(config)
    # check config
    check_configuration(params)
    print('read and checked config')
    
    # check mapfile
    mapfile = params['Options']['mapfile']
    num_columns = check_mapfile(params, mapfile)
    print('checked map file: {0} columns'.format(num_columns))
        
    # create output directory and output sub-folders. remove output directory if it exists
    outdir = params['Options']['outdir']
    cbiodir, casedir, suppdir = create_output_directories(outdir)
    print('created output directories')

    # copy config file to out directory
    copy_resource(args.config, outdir)
    copy_resource(mapfile, outdir)
    print('copied {0} to {1}'.format(os.path.basename(args.config), outdir))
    print('copied {0} to {1}'.format(os.path.basename(mapfile), outdir))
     
    # create input directories for each data type 
    mafdir, segdir, gepdir, fusdir = create_input_directories(outdir)
    print('created input directories')
    
    # preprocess the map file; extract sequenza gamma segments or generate purple segment files
    processed_mapfile = preprocess_mapfile(params, mapfile, params['Options']['gamma'], outdir)
    print('preprocessed map file: {0}'.format(processed_mapfile))

    # check that input maf files, if any, have the same format and the same header
    check_input_mafs(processed_mapfile)
    print('checked maf headers')
    
    # check genome version in the maf files, if provided
    genome = params['Options']['genome']
    check_genome_version(mapfile, genome)
    print('checked genome version')
    
    # get genome specific variables
    if genome == 'hg38':
        enscon, genebed = params['Resources']['enscon_hg38'], params['Resources']['genebed_hg38']
    elif genome == 'hg19':
        enscon, genebed = params['Resources']['enscon_hg19'], params['Resources']['genebed_hg19']
    print('determined genome specific resources')
        
    # check that cancer type is correctly defined
    cancer_code = params['Options']['cancer_code']
    check_cancer_type(cancer_code)
    print('checked cancer code')

    # write meta study and clinical files
    project_name = params['Options']['project_name']
    write_meta_study(os.path.join(cbiodir, 'meta_study.txt') ,
                     params['Options']['study'],
                     project_name,
                     params['Options']['description'],
                     genome, cancer_code)
    write_meta_clinical(cbiodir, project_name, 'sample')
    write_meta_clinical(cbiodir, project_name, 'patient')
    print('wrote study and clinical metadata')    
    
    # write cases
    write_cases(os.path.join(casedir, 'cases_sequenced.txt'), project_name, processed_mapfile, 'seq')
    write_cases(os.path.join(casedir, 'cases_rna_seq_mrna.txt'), project_name, processed_mapfile, 'rna')
    write_cases(os.path.join(casedir, 'cases_cna.txt'), project_name, processed_mapfile, 'cna')
    write_cases(os.path.join(casedir, 'cases_cnaseq.txt'), project_name, processed_mapfile, 'cna_seq')
    write_cases(os.path.join(casedir, 'cases_3way_complete.txt'), project_name, processed_mapfile, 'cna_seq_rna')
    write_cases(os.path.join(casedir, 'cases_sv.txt'), project_name, processed_mapfile, 'sv')
    print('wrote cases')

    # write patient clinical information
    clinical_outputfile = os.path.join(cbiodir, 'data_clinical_patients.txt')
    center = params['Options']['center']
    write_patient_minimal_clinical_information(clinical_outputfile, processed_mapfile, center)
    
    # get the user defined sample clinical information
    clinical_info = get_clinical_data(args.clinical) if args.clinical else None
    clinical_outputfile = os.path.join(cbiodir, 'data_clinical_samples.txt')
    write_sample_minimal_clinical_information(clinical_outputfile, processed_mapfile, center, sample_info = clinical_info)
    print('wrote clinical information')
    
    # write clinical input file for oncokb-annotator
    clinical_oncokb = os.path.join(suppdir, 'oncokb_clinical_info.txt')
    write_clinical_oncokb(clinical_oncokb, processed_mapfile, cancer_code)
    print('wrote oncoKb clinical information') 
    
    # concatenate data files if they exist  
    # parse processed map file
    data = parse_processed_mapfile(processed_mapfile)
    # concatenate VEP mafiles
    mutation_file = concatenate_maf_files(data, os.path.join(mafdir, 'all_mutations.maf.txt'))
    # concatenate segmentation files
    segfile = concatenate_seg_files(data, os.path.join(segdir, 'input.seg.txt'))
    # concatenate fusion files
    fusfile = concatenate_fusion_files(data, os.path.join(fusdir, 'input.fus.txt'))
    # concatenate expression files
    # extract and concatenate fpkm from gep files
    if params['Workflows']['expression'] == 'rsem':
        count = params['Options']['count']
        assert count in ['fpkm', 'tpm']
        gepfile = concatenate_expression_from_gep_files(data, count, os.path.join(gepdir, 'input.expression.txt'))
    elif params['Workflows']['expression'] == 'isofox':
        gepfile = concatenate_expression_from_isofox_files(data, os.path.join(gepdir, 'input.expression.txt'))
    print('concatenated data files')

    
    
    #### to do:
    
    ### concantenate hmf mafs
    ### concatenate hmf fusion
    ### concatenate hmf segments
    
       
    # filter maf files and write metadata
    # define maffile, output of MafAnnotator
    filter_variants = params['Filters']['filter_variants']
    depth_filter = params['Filters']['depth_filter']
    alt_freq_filter = params['Filters']['alt_freq_filter']
    gnomAD_AF_filter = params['Filters']['gnomad_af_filter'] 
    filter_indels = params['Filters']['filter_indels']
    keep_variants = params['Options']['keep_variants']
    tglpipe = params['Filters']['tglpipe']    
    
    
    if mutation_file:
        maffile = os.path.join(mafdir, 'input.maf.txt')
        # filter mutations and indels if option is activated
        if filter_variants:
            total, kept = filter_mutations(mutation_file, os.path.join(mafdir, 'filtered.mutations.txt'), depth_filter, alt_freq_filter, gnomAD_AF_filter, keep_variants, os.path.join(mafdir, 'filtered.mutations.removed.txt'))
            print("before mutations filtering: ", total)
            print("after mutations filtering: ", kept)
            print('filtered variants')
        # filter indels if option is activated
        if filter_indels:
            # check if variants are filtered
            if filter_variants:
                maf_to_filter = os.path.join(mafdir, 'filtered.mutations.txt')
                maf_filtered = os.path.join(mafdir, 'filtered.mutations.indels.txt')
                maf_rem = os.path.join(mafdir, 'filtered.mutations.indels.remove.txt')
            else:
                maf_to_filter = os.path.join(mafdir, 'all_mutations.maf.txt')
                maf_filtered = os.path.join(mafdir, 'filtered.indels.txt')
                maf_rem = os.path.join(mafdir, 'filtered.indels.remove.txt')
            # output file for MafAnnotator
            maf_input_annotation = maf_filtered
            total, kept = remove_indels(maf_to_filter, maf_filtered, maf_rem)
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
        oncokb_token = get_token(params['Resources']['token'])
    
        # annotate mafs with oncokb-annotate        
        maf_annotation = subprocess.call('MafAnnotator -i {0} -o {1} -c {2} -b {3}'.format(maf_input_annotation, maffile, os.path.join(suppdir, 'oncokb_clinical_info.txt'), oncokb_token), shell=True)
        # check exit code
        if maf_annotation:
            sys.exit('Error when running MafAnnotator.')
        else:
            print('Annotated variants with MafAnnotator')
    
        # generate mutations data
        print('Processing Mutation data from {0}'.format(maffile))
        # only do the filtering steps if tglpipe is set to TRUE
        if tglpipe:
            print('tglpipe is set to true, filtering data according to tgl specifications')
            df_cbio_anno = procVEP(maffile)
            df_cbio_filt = df_cbio_anno[df_cbio_anno['TGL_FILTER_VERDICT'] == 'PASS']
            # get snvs for dcsigs
            df_snv = df_cbio_filt[df_cbio_filt['Variant_Type'] == 'SNP']
            # for cbioportal input
            df_cbio_filt.to_csv(os.path.join(cbiodir, 'data_mutations_extended.txt'), sep="\t", index=False, na_rep='NA')
            # unfiltered data
            df_cbio_anno.to_csv(os.path.join(suppdir, 'unfiltered_data_mutations_extended.txt'), sep="\t", index=False, na_rep='NA')
        else:
            df_cbio_filt = pd.read_csv(maffile, sep="\t", header=0)
            df_snv = df_cbio_filt[df_cbio_filt['Variant_Type'] == 'SNP']
            df_cbio_filt.to_csv(os.path.join(cbiodir, 'data_mutations_extended.txt'), sep="\t", index=False, na_rep='NA')

        # write metadata file
        write_metadata(os.path.join(cbiodir, 'meta_mutations_extended.txt'), project_name, 'maf', genome)
        print('wrote mutations metadata')
    
    
    gain = params['Parameters']['gain']
    amplification = params['Parameters']['amplification']
    heterozygous_deletion = params['Parameters']['heterozygous_deletion']
    homozygous_deletion = params['Parameters']['homozygous_deletion']
    oncolist = params['Resources']['oncolist']
    genelist = params['Resources']['genelist']
       
    # generate CNA data and metadata files if input segmentation file exists
    if segfile:
        # generate metadata files
        # write cna metadata
        write_metadata(os.path.join(cbiodir, 'meta_CNA.txt'), project_name, 'discrete', genome)    
        write_metadata(os.path.join(cbiodir, 'meta_log2CNA.txt'), project_name, 'log2-value', genome)    
        write_metadata(os.path.join(cbiodir, 'meta_segments.txt'), project_name, 'seg', genome) 
        print('wrote CNA metadata files')
        if genelist:
            print('Restricting CNAs to the list of genes provided in {0}'.format(genelist))
          
        # function returns list of 3 objects
        print('Processing CNA data from {0}'.format(segfile))
        segData, df_cna, df_cna_thresh = preProcCNA(segfile, genebed, gain, amplification, heterozygous_deletion, homozygous_deletion, oncolist, genelist)
        # write cbio files
        print('writing seg file')
        segData.to_csv(os.path.join(cbiodir, 'data_segments.txt'), sep='\t', index=False)
        df_cna.to_csv(os.path.join(cbiodir, 'data_log2CNA.txt'), sep='\t', index=True)
        df_cna_thresh.to_csv(os.path.join(cbiodir, 'data_CNA.txt'), sep='\t', index=True)
        # write the truncated data_CNA file (remove genes which are all zero) for oncoKB annotator
        df_CNA = df_cna_thresh.loc[~(df_cna_thresh == 0).all(axis=1)]
        df_CNA.to_csv(os.path.join(suppdir, 'data_CNA_short.txt'), sep='\t', index=True)
    
    
    # generate expression data and metadata file if input file exists
    if gepfile:
        # write metadata files
        write_metadata(os.path.join(cbiodir, 'meta_expression.txt'), project_name, 'expression', genome)
        write_metadata(os.path.join(cbiodir, 'meta_expression_zscores.txt'), project_name, 'zscore', genome)
        print('wrote expression metadata files')
        # write all samples with rna data to file 
        gep_study_file = os.path.join(outdir, 'gep_study.list')
        expression_samples_to_file(data, gep_study_file)
        # generate expression data files
        print('Processing RNASEQ data from {0}'.format(gepfile))
        # get list of samples in study with expression data
        study_samples = list_samples_with_expression(data)
        # preprocess the full data frame
        df = preProcRNA(gepfile, enscon, genelist)
        print('getting STUDY-level data')
        # subset data to STUDY level data for cbioportal
        df_study = df[study_samples]
        # write the raw STUDY data
        df_study.to_csv(os.path.join(cbiodir, "data_expression.txt"), sep="\t", header=True, index=True)
        # z-scores STUDY
        df_zscore = compZ(df_study)
        df_zscore.to_csv(os.path.join(cbiodir, "data_expression_zscores.txt"), sep="\t", header=True, index=True)
        print('wrote expression data files')        
    
    # generate fusion data and metadata if input file exists
    entcon = params['Resources']['entcon']
    minfusionreads = params['Parameters']['minfusionreads']
    
    if fusfile:
        # write SV metadata
        write_metadata(os.path.join(cbiodir, 'meta_sv.txt'), project_name, 'sv', genome)
        print('wrote SV metadata')
        # process the fusion file
        fusion_cbio = preProcFus(fusfile, minfusionreads, entcon)
        # write fusion file
        data_fusion = os.path.join(cbiodir, "data_fusions.txt")
        fusion_cbio.to_csv(data_fusion, sep="\t", index=False)
        print('wrote fusion data files')
        # convert fusion file to SV format
        # get the path to sv file
        data_sv = os.path.join(cbiodir, "data_sv.txt")
        # convert to sv file
        convert_fusion_to_sv(data_fusion, data_sv)
        # move fusion file to supplementary directory
        if os.path.isfile(data_fusion):
            new_data_fusion = os.path.join(suppdir, os.path.basename(data_fusion))
            os.rename(data_fusion, new_data_fusion)     
                    
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


def merge_import_folder(args):
    '''
    (list) -> None
    
        
    Parameters
    ----------
    - configs (list | None): List of paths to config files
    - import_folders (list | None): List of paths to impot folders
    - config_paths (str | None): File with paths to the config files of the import folders to merge
    - import_paths (str): File with paths to the import folders to merge
    - cancer_code (str): Cancer code. See http://oncotree.mskcc.org
    - genome (str): Reference genome. hg19 or hg38
    - merge_with_differences (bool): Merges import folders with parameter differences if True. Default is False
    - outdir (str): Path to the output directory
    - center (str): Genomic center
    - project (str): Name of the project
    - description (str): Project description
    - study (str): Format of the study name: Project: Top-level-OncoTree, Concept (PI, Centre) AcronymProject description
    - exclude_samples (list | None): list of excluded samples
    - exclude_file (str | None): Path to file with samples to exclude
    '''
    
    
    # check that config files and import folders are specified
    if (args.configs and args.config_paths) or (args.configs is None and args.config_paths is None):
        sys.exit('Use -configs or -configPaths to specify the paths to the config files')    
    if (args.import_folders and args.import_paths) or (args.import_folders is None and args.import_paths is None):
        sys.exit('Use -importFolders or -importPaths to specify the paths to the import folders')    
       
    if  args.exclude_samples and args.exclude_file:
        sys.exit('Use -excludeSamples or -excludeFile to remove samples from the final import folder')    
       
    # get import folders and config files
    if args.configs:
        config_files = args.configs
    elif args.config_paths:
        config_files = get_paths(args.config_paths)
    if args.import_folders:
        import_folders = args.import_folders
    elif args.import_paths:
        import_folders = get_paths(args.import_paths)

    # check that config files and import folders are valid
    check_configs_import_folders(config_files, import_folders)
      
    # collect parameters
    config_params = collect_config_parameters(config_files)
    
    # check the parameters in the configs
    same_parameters, differences = check_config_parameters(config_params)
    print('collected config parameters')

    # check genome
    check_genome(config_params, args.genome)
    print('checked genome')
    
    # check workflows
    check_datatypes(config_params)
    print('checked data types')
    
    # check resources
    check_metadata(config_params, args.cancer_code, 'cancer_code')
    check_metadata(config_params, args.project, 'project_name')
    check_metadata(config_params, args.study, 'study')
    check_metadata(config_params, args.description, 'description')
    print('checked metadata')
    check_cancer_type(args.cancer_code)
    print('checked cancer code')


    # check if merging is allowed if differences exist
    if same_parameters == False and args.merge_with_differences == False:
        sys.exit('''WARNING: Attempting to merge import folders with parameter differences.
                 Parameters with differences in config files: {0}
                 Use flag --merge_with_differences to force merging'''.format(';'.join(differences)))
    else:
        if same_parameters:
            print('Merging import folders with identical parameters')
        elif args.merge_with_differences:
            print('''WARNING: Merging import folders with parameter differences.
                  Parameters with differences in config files: {0}'''.format(';'.join(differences)))
                 
        # create output folders
        cbiodir = os.path.join(args.outdir, 'cbioportal_import_data')
        casedir = os.path.join(cbiodir, 'case_lists')
        for i in [cbiodir, casedir]:
            os.makedirs(i, exist_ok=True)
        print('created output sub-folders:', cbiodir, casedir, sep = '\n')
                
        # get the list of samples to exclude
        if  args.exclude_samples:
            excluded_samples = args.exclude_samples
        elif args.exclude_file:
            infile = open(args.exclude_file)
            excluded_samples = infile.read().rstrip().split('\n')
            infile.close()
        else:
            excluded_samples = []
        if excluded_samples:
            print('Excluding {0} samples from the final data files'.format(len(excluded_samples)))
        
            
        # get the list of files to merge in the import folder 
        cbiofiles = group_files(import_folders, False)
            
        # check if the data files headers are the same (ie, files can be merged) 
        check_file_headers(cbiofiles) 
        
        # write the case files in the cbioportal_import_data/case_lists folder
        # remove samples from the case files if samples in excluded samples
        merge_case_files(casedir, import_folders, args.project, excluded_samples)
        print('Wrote case files')
        
        # write the metadata files in the cbioportal_import_data folder
        merge_meta_files(cbiodir, import_folders, args.project, args.study, args.description, args.genome, args.cancer_code)
        print('Wrote metadata files')

        # merge mutation, sv and segment data files
        # remove samples from the data files if samples in excluded samples
        merge_data_files(cbiodir, cbiofiles, excluded_samples)
        print('Merged mutation and SV data files')
        
        # merge clinical data files
        # remove samples from the clinical files if samples in excluded samples
        merge_clinical_files(cbiodir, cbiofiles, excluded_samples)
        print('Merged clinical data files')
        
        # merge expression and CNA data files
        # remove samples from the expression and CNA files if samples in excluded samples
        merge_expression_CNA_files(cbiodir, cbiofiles, excluded_samples)
        print('Merged expression and CNA data files')

        
        # remove case files if datafiles are not present because samples were excluded
        remove_merged_metafiles(cbiodir)
        
            
        # THINGS TO DO
        # allow clinical merging of clinical files with different headers:
            # clinical_samples and clinical_patients may have custom fields: 
            # enable merging of clinical files with different hears and fields
            # do not include clinical data files when checking for header


    print('Success! {0} import folders have been merged into {1}'.format(len(import_folders), cbiodir))

 

def generate_mapfile(args):
    '''
    (str, str)
    
    Generates a map file with data required for creating the cbioportal import folder
       
    Parameters
    ----------
    - outdir (str): Path to the output directory
    - input_data (str): Path to the json file with input data
    - gamma (str): Gamma value of the sequenza workflow. Default is 500
    '''
    
    infile = open(args.input_data)
    data = json.load(infile)
    infile.close()
    
    # create directories
    os.makedirs(args.outdir, exist_ok=True)
        
    mapfile = os.path.join(args.outdir, 'map.csv')
    newfile = open(mapfile, 'w')
        
    for donor in data:
        for sample in data[donor]:
            L = [donor,sample,'NA','NA','NA','NA']
            for workflow in data[donor][sample]:
                if 'sequenza' in workflow.lower():
                    # create a sequenza directory
                    sequenzadir =  os.path.join(args.outdir, 'sequenza')
                    donordir = os.path.join(sequenzadir, donor)
                    sampledir = os.path.join(donordir, sample)
                    os.makedirs(sampledir, exist_ok=True)
                    if os.path.isfile(data[donor][sample][workflow]):
                        with zipfile.ZipFile(data[donor][sample][workflow], 'r') as zipref:
                            zipref.extractall(sampledir)
                    # get the seg file of interest
                    segfile = glob.glob(sampledir +'/gammas/{0}/*.seg'.format(args.gamma))
                    assert len(segfile) == 1
                    segfile = segfile[0]
                    assert os.path.isfile(segfile)
                    L[3] = segfile       
                elif 'purple' in workflow.lower():
                    # create a purple directory
                    purpledir =  os.path.join(args.outdir, 'purple')
                    donordir = os.path.join(purpledir, donor)
                    sampledir = os.path.join(donordir, sample)
                    os.makedirs(sampledir, exist_ok=True)
                    #get the cnv and purity files
                    cnv_file = data[donor][sample][workflow]['cnv']
                    purity_file = data[donor][sample][workflow]['purity']
                    # generate segmentation file from purple
                    outputfile = os.path.join(sampledir, '{0}.purple.seg'.format(sample))
                    generate_purple_segmentation(cnv_file, purity_file, sample, outputfile)
                    L[3] = outputfile
                elif 'varianteffectpredictor' in workflow.lower():
                    L[2] = data[donor][sample][workflow]
                elif 'rsem' in workflow.lower():
                    L[4] = data[donor][sample][workflow]
                elif 'mavis' in workflow.lower():
                    L[5] = data[donor][sample][workflow]
            newfile.write(','.join(L) + '\n')


    newfile.close()
 



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
    
    # merge import folders
    mg_parser = subparsers.add_parser('merge', help="Merge previous import folders")
    mg_parser.add_argument('-configs', dest='configs', nargs = '*', help="White space separated list of paths to config files")
    mg_parser.add_argument('-importFolders', dest='import_folders', nargs = '*', help="White space separated list of paths to impot folders")
    mg_parser.add_argument('-configPaths', dest='config_paths', help="File with paths to the config files of the import folders to merge")
    mg_parser.add_argument('-importPaths', dest='import_paths', help="File with paths to the import folders to merge")
    mg_parser.add_argument('-cancerCode', dest='cancer_code', help="Cancer code. See http://oncotree.mskcc.org", required=True)
    mg_parser.add_argument('-genome', dest='genome', choices = ['hg19', 'hg38'], help="Reference genome", required=True)
    mg_parser.add_argument('--merge_with_differences', dest='merge_with_differences', action='store_true', help='Merges import folders with parameter differences if activated. Excludes genome differences. Default is False')
    mg_parser.add_argument('-outdir', dest='outdir', help='Path to the output directory', required = True)
    mg_parser.add_argument('-center', dest='center', help='genomic center', required = True)
    mg_parser.add_argument('-project', dest='project', help='Name of the project', required = True)
    mg_parser.add_argument('-description', dest='description', help='Project description', required = True)
    mg_parser.add_argument('-study', dest='study', help='Format of the study name: Project: Top-level-OncoTree, Concept (PI, Centre) AcronymProject description', required = True)
    mg_parser.add_argument('-excludeSamples', dest='exclude_samples', nargs= '*',  help='White space sparated list of excluded samples')
    mg_parser.add_argument('-excludeFile', dest='exclude_file', help='Path to file with samples to exclude')
    mg_parser.set_defaults(func=merge_import_folder)
        
    # generate maping file
    m_parser = subparsers.add_parser('map', help="Generate map file")
    m_parser.add_argument('-o', '--outdir', dest='outdir', help='Path to the output directory', required = True)
    m_parser.add_argument('-i', '--input_data', dest='input_data', help='Path to the json file with input data', required = True)
    m_parser.add_argument('-g', '--gamma', dest='gamma', default='500', help='Gamma value of the sequenza output. Default is 500')
    m_parser.set_defaults(func=generate_mapfile)
    
    # get arguments from the command line
    args = parser.parse_args()
    # pass the args to the default function
    args.func(args)
