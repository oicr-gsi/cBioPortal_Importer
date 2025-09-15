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
#import pandas as pd 
import gc 
#from rpy2 import robjects 
#from rpy2.robjects import pandas2ri
#from rpy2.robjects.packages import importr
#from scipy.stats import zscore

import warnings


#pandas2ri.activate()
#cntools = importr('CNTools')



def extract_files_from_map(mapfile, data_type):
    '''
    (str, str) -> list
    
    Returns a list of input files from the mapping file
        
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


# def extract_purple_ploidy_purity(purple_purity_file):
#     '''
#     (str) -> (float, float)
    
#     Returns a tuple with purity and ploidy extracted from file .purple.purity.tsv
 
#     Parameters
#     ----------
#     - purple_purity_file (str): Path to the Purple .purple.purity.tsv file
#     '''
    
#     infile = open(purple_purity_file)
#     header = infile.readline().rstrip().split('\t')
#     content = infile.readline().rstrip().split('\t')
#     infile.close()
    
#     purity = float(content[header.index('purity')])
#     ploidy = float(content[header.index('ploidy')])

#     return purity, ploidy
    

# def convert_purple_to_seg(somatic_file, purity, ploidy, sample_name):
#     '''
#     (str, str) -> list
    
#     Returns a list of fields extracted from the purple somatic file
#     corresponding to the sequenza segmentation output file 
    
#     Fields are mapped and renamed as:
#     chromosome -> chrom
#     start -> loc.start
#     end -> loc.end 
#     bafCount -> num.mark
#     copyNumber -> seg.mean adjusted with purity and ploidy  
    
#     Parameters
#     ----------
#     - somatic_file (str): Path to the purple "purple.cnv.somatic.tsv" file
#     - sample_name (str): Name of the sample
#     - purity (float): Purity estimated from purple
#     - ploidy (float): Ploidy estimated from purple
#     '''
    
    # L = [['chrom', 'loc.start', 'loc.end', 'num.mark', 'seg.mean']]
    
    # infile = open(somatic_file)
    # header = infile.readline().rstrip().split('\t')
    # for line in infile:
    #     line = line.rstrip()
    #     if line:
    #         line = line.split('\t')
    #         # compute seg.mean       
    #         copyNumber = float(line[header.index('copyNumber')])
            
    #         with warnings.catch_warnings():
    #             warnings.simplefilter("error", RuntimeWarning)
    #             try:
    #                 adjusted_copy_number = np.log2(1 + (purity *(copyNumber - ploidy)/ploidy))
    #             except:
    #                 adjusted_copy_number = np.nan
    #             finally:
    #                 L.append([sample_name, line[header.index('chromosome')],
    #                           line[header.index('start')],
    #                           line[header.index('end')],
    #                           line[header.index('bafCount')],
    #                           str(adjusted_copy_number)])
    # infile.close()
    
    # return L


# def generate_purple_segmentation(somatic_file, purple_purity_file, sample_name, outputfile):
#     '''
#     (str, str, str) -> None
    
#     Convert the purple somatic file to a segmentation file and write to outputfile
    
#     Parameters
#     ----------
#     - somatic_file (str): Path to the purple "purple.cnv.somatic.tsv" file
#     - purple_purity_file (str): Path to the Purple .purple.purity.tsv file
#     - sample_name (str): Name of the sample
#     - outputfile (str): Path to the outputfile
#     '''

#     purity, ploidy = extract_purple_ploidy_purity(purple_purity_file)
    
#     L = convert_purple_to_seg(somatic_file, purity, ploidy, sample_name)
#     newfile = open(outputfile, 'w')
#     for i in L:
#         newfile.write('\t'.join(i) + '\n')
#     newfile.close()


# def is_sequenza_segmentation(segfile):
#     '''
#     (str) -> bool
    
#     Returns True if the segfile is the segmentation file from sequenza
    
#     Parameters
#     ----------
#     - segfile (str):Path to a segmentation file
#     '''
    
#     infile = open(segfile)
#     header = infile.readline().rstrip().split('\t')
#     infile.close()
    
#     return header == ['ID', 'chrom', 'loc.start', 'loc.end', 'num.mark', 'seg.mean']


# def split_column_take_max(df, columns):
#     '''
#     Splits columns by ';' and takes the maximum value from a list of columns
#     Parameters
#     ----------
#     - df (pd.DataFrame) : The input data that contains the columns to be processed
#     - columns (list) : A list of columns names (str) to be processed
#     Returns
#     -------
#     - df (pd.DataFrame) : The modified dataframe 
#     '''
#     df[columns] = df[columns].fillna(0).replace('', 0)

#     for column in columns:
#         df[column] = df[column].apply(
#             lambda x: max(map(int, x.split(';'))) if isinstance(x, str) and ';' in x else int(x)
#         )
    
#     for column in columns:
#         df[column] = pd.to_numeric(df[column], errors='coerce')  

#     return df


# def preProcFus(datafile, readfilt, entrfile):
#     '''
#     Preprocesses the fusion input data and formats it for CBioPortal output. 
#     Parameters
#     -----------
#     - datafile (str) : Path to the input fusion data file
#     - readfilt (int) : Minimum number of reads for fusion calls
#     - entrfilw (str) : Path to the Entrez gene ID file
#     Returns
#     --------
#     - df_cbio (pd.DataFrame) : Processed and formatted data
#     '''
#     data = pd.read_csv(datafile, sep="\t")
#     entr = pd.read_csv(entrfile, sep="\t")

#     # reformat the filtering columns to split and take the max value within cell
#     columns = ['contig_remapped_reads', 'flanking_pairs', 'break1_split_reads', 'break2_split_reads', 'linking_split_reads']
#     data = split_column_take_max(data, columns)

#     # add a column which pulls the correct read support column    
#     data['read_support'] = data.apply(
#     lambda row: row['contig_remapped_reads'] if 'contig' in row['call_method'] 
#     else row['flanking_pairs'] if 'flanking reads' in row['call_method'] 
#     else (row['break1_split_reads'] + row['break2_split_reads'] + row['linking_split_reads']) if 'split reads' in row['call_method'] 
#     else 0, 
#     axis=1)

#     # filter by minimum read support
#     data = data[data['read_support'] > readfilt]

#     # sort descending read support 
#     data = data.sort_values(by='read_support', ascending=False)

#     # get unique fusions for each sample
#     data['gene1_aliases'] = data['gene1_aliases'].fillna('')
#     data['gene2_aliases'] = data['gene2_aliases'].fillna('')

#     data['fusion_tuples'] = data.apply(
#     lambda row: '-'.join(sorted([str(x) if x else 'None' for x in [row['gene1_aliases'], row['gene2_aliases']]])), axis=1)

#     # add index which is sample, tuple
#     data['index'] = data['Sample'] + data['fusion_tuples']

#     # deduplicate
#     data_dedup = data.drop_duplicates(subset='index', keep='first')

#     # gene1 should not equal gene2
#     data_dedup = data_dedup[data_dedup['gene1_aliases'] != data_dedup['gene2_aliases']]

#     # merge in entrez gene ids
#     data_dedup = pd.merge(data_dedup, entr, how='left', left_on='gene1_aliases', right_on='Hugo_Symbol')
#     data_dedup = pd.merge(data_dedup, entr, how='left', left_on='gene2_aliases', right_on='Hugo_Symbol', suffixes=('.x', '.y'))

#     # add some missing columns
#     data_dedup['DNA_support'] = 'no'
#     data_dedup['RNA_support'] = 'yes'
#     data_dedup['Center'] = 'TGL'
#     data_dedup['Frame'] = 'frameshift'
#     data_dedup['Fusion_Status'] = 'unknown'

#     # write out the nice header 
#     header = ['Hugo_Symbol', 'Entrez_Gene_Id', 'Center', 'Tumor_Sample_Barcode', 'Fusion', 'DNA_support', 'RNA_support', 'Method', 'Frame', 'Fusion_Status']

#     # get left gene data
#     col_left = ['gene1_aliases', 'Entrez_Gene_Id.x', 'Center', 'Sample', 'fusion_tuples', 'DNA_support', 'RNA_support', 'tools', 'Frame', 'Fusion_Status']
#     data_left = data_dedup[col_left]
#     data_left.columns = header

#     # get right gene data
#     col_right = ['gene2_aliases', 'Entrez_Gene_Id.y', 'Center', 'Sample', 'fusion_tuples', 'DNA_support', 'RNA_support', 'tools', 'Frame', 'Fusion_Status']
#     data_right = data_dedup[col_right]
#     data_right.columns = header

#     data_left.loc[:, 'Hugo_Symbol'] = data_left['Hugo_Symbol'].replace("", float('nan')).fillna(data_right['Hugo_Symbol'])
#     data_left.loc[:, 'Entrez_Gene_Id'] = data_left['Entrez_Gene_Id'].fillna(data_right['Entrez_Gene_Id'])
#     data_right.loc[:, 'Hugo_Symbol'] = data_right['Hugo_Symbol'].replace("", float('nan')).fillna(data_left['Hugo_Symbol'])
#     data_right.loc[:, 'Entrez_Gene_Id'] = data_right['Entrez_Gene_Id'].fillna(data_left['Entrez_Gene_Id'])

#     # append it all together
#     df_cbio = pd.concat([data_left, data_right])
#     df_cbio['Entrez_Gene_Id'] = df_cbio['Entrez_Gene_Id'].astype(int)

#     # remove rows where gene is not known (this still keeps the side of the gene which is known)
#     df_cbio = df_cbio.dropna()

#     return df_cbio



# def preProcCNA(segfile, genebed, gain, amp, htz, hmz, oncolist, genelist=None):
#     '''
#     Processes CNA data by applying thresholds and performing gene-level segmentation.
#     Parameters
#     ----------
#     - segfile (str) : Path to the input concatenated segmentation file from Sequenza
#     - genebed (str) : Path to a tab-delimited bed file defining the genomic positions of canonical genes
#     - gain (float) : Threshold for gains in CNA data
#     - amp (float) : Threshold for amplification in CNA data
#     - htz (float) : Threshold for heterozygous deletion in CNA data
#     - hmz (float) : Threshold for homozygous deletion in CNA data
#     - oncolist (str) : Path to a tab-delimited file containing the list of cancer genes
#     - genelist (str or None) : (Optional) Path to a file containing a list of Hugo Symbols to filter the genes, or None if no filtering is required
#     Returns
#     -------
#     - segData (pd.DataFrame) : Dataframe containing the processed segmentation data
#     - df_cna (pd.DataFrame) : Dataframe with gene-level log2 copy number alterations
#     - df_cna_thresh (pd.DataFrame) : Dataframe with thresholded CNA values (5-state matrix)
#     '''
#     # check if genelist is None when preProcCNA.py is called by pycbio.py and genelist is omitted
#     if genelist:
#         print('genelist is used during CNA processing')
#     else:
#         print('genelist is not used during CNA processing')
        
    
#     # read oncogenes
#     oncogenes = pd.read_csv(oncolist, sep='\t')

#     # small fix segmentation data
#     segData = pd.read_csv(segfile, sep='\t')
#     segData['chrom'] = segData['chrom'].str.replace('chr', '')

#     # remove NA
#     segData.dropna(inplace=True)

#     # convert pandas dataframe to R dataframe 
#     segData_r = pandas2ri.py2rpy(segData)

#     # set thresholds
#     print('setting thresholds')
#     gain = float(gain)
#     amp = float(amp)
#     htz = float(htz)
#     hmz = float(hmz)

#     # get the gene info
#     print('getting gene info')
#     geneInfo = pd.read_csv(genebed, sep='\t')

#     # convert pandas dataframe to R dataframe 
#     geneInfo_r = pandas2ri.py2rpy(geneInfo)

#     # make CN matrix gene level
#     print('converting seg')
#     cnseg = cntools.CNSeg(segData_r)
#     print('get segmentation by gene')
#     rdByGene = cntools.getRS(cnseg, by='gene', imput=False, XY=False, geneMap=geneInfo_r, what='median', mapChrom='chrom', mapStart='start', mapEnd='end')
#     print('get reduced segmentation data')
#     reducedseg_df = cntools.rs(rdByGene)

#     # convert data from R back to Py
#     reducedseg = pandas2ri.rpy2py(reducedseg_df)

#     # some reformatting and return log2cna data
#     df_cna = reducedseg.iloc[:, [4] + list(range(5, reducedseg.shape[1]))]
#     df_cna = df_cna.drop_duplicates(subset=[df_cna.columns[0]])
#     df_cna.columns = ['Hugo_Symbol'] + list(df_cna.columns[1:])

#     # set thresholds and return 5-state matrix
#     print('thresholding cnas')
#     df_cna_thresh = df_cna.copy()
#     df_cna_thresh.iloc[:, 1:] = df_cna_thresh.iloc[:, 1:].apply(pd.to_numeric)

#     # threshold data
#     for col in df_cna_thresh.columns[1:]:
#         df_cna_thresh[col] = df_cna_thresh[col].apply(lambda x: 2 if x > amp
#                                                 else (-2 if x < hmz
#                                                       else (1 if gain < x <= amp
#                                                             else (-1 if hmz <= x < htz
#                                                                   else 0)
#                                                         )
#                                                     )
#                                                 )
    
#     # fix rownames of log2cna data
#     df_cna.set_index('Hugo_Symbol', inplace=True)
#     df_cna = df_cna.round(4)
#     df_cna_thresh.set_index(df_cna_thresh.columns[0], inplace=True)

#     # subset of oncoKB genes
#     df_cna_thresh_onco = df_cna_thresh.loc[df_cna_thresh.index.isin(oncogenes.index)]

#     # subset if gene list is given
#     if genelist is not None:
#         with open(genelist, 'r') as file:
#             keep_genes = [line.strip() for line in file]
            
#         df_cna = df_cna.loc[df_cna.index.isin(keep_genes)]
#         df_cna_thresh = df_cna_thresh.loc[df_cna_thresh.index.isin(keep_genes)]
      
#     return segData, df_cna, df_cna_thresh



# def addVAFtoMAF(maf_df, alt_col, dep_col, vaf_header):
#     '''
#     Adds a VAF column to a MAF DataFrame.
#     Parameters
#     ----------
#     - maf_df (pd.DataFrame): The MAF DataFrame containing mutation data.
#     - alt_col (str): The name of the column representing the alternate allele count.
#     - dep_col (str): The name of the column representing the depth. 
#     - vaf_header (str): The name for the new column that will hold the VAF values.
#     Returns
#     -------
#     - maf_df (pd.DataFrame) : The modified dataframe with the VAF column. 
#     '''
#     # print a warning if any values are missing (shouldn't happen), but change them to 0
#     if maf_df[alt_col].isnull().any() or maf_df[dep_col].isnull().any():
#         print('Warning! Missing values found in one of the count columns')
#         maf_df[alt_col] = maf_df[alt_col].fillna(0)
#         maf_df[dep_col] = maf_df[dep_col].fillna(0)
    
#     # ensure factors end up as numeric
#     maf_df[alt_col] = pd.to_numeric(maf_df[alt_col], errors='coerce')
#     maf_df[dep_col] = pd.to_numeric(maf_df[dep_col], errors='coerce')

#     # ensure position comes after alternate count field 
#     bspot = maf_df.columns.get_loc(alt_col)
#     maf_df.insert(bspot + 1, vaf_header, maf_df[alt_col] / maf_df[dep_col])

#     # check for any NAs
#     if maf_df[vaf_header].isnull().any():
#         print('Warning! There are missing values in the new vaf column')
#         maf_df[vaf_header] = maf_df[vaf_header].fillna(0)
    
#     return maf_df
    

# def procVEP(datafile):
#     '''
#     Processes the input MAF file, adding various computed columns, and applying multiple filters to prepare the data for further analysis.
#     Parameters
#     ----------
#     - datafile (str): The file path to the input in tab-separated format.
#     Returns
#     -------
#     - df_anno (pd.DataFrame) : The modified dataframe. 
#     '''
#     print("--- reading data ---")
#     data = pd.read_csv(datafile, sep="\t")

#     print('--- doing some formatting ---')

#     # add vaf columns 
#     print('add tumor_vaf')
#     data = addVAFtoMAF(data, 't_alt_count', 't_depth', 'tumor_vaf')
#     print('add normal_vaf')
#     data = addVAFtoMAF(data, 'n_alt_count', 'n_depth', 'normal_vaf')

#     # clear memory (important when the mafs are huge with millions and millions of lines)
#     df_anno = data
#     del data
#     gc.collect()

#     # add oncogenic yes or no columns 
#     print('add oncogenic status')
#     df_anno['oncogenic_binary'] = df_anno['oncogenic'].apply(lambda x: 'YES' if x in ['Oncogenic', 'Likely Oncogenic'] else 'NO')

#     # add common_variant yes or no columns
#     df_anno['ExAC_common'] = df_anno['FILTER'].apply(lambda x: 'YES' if 'common_variant' in x else 'NO')

#     # add POPMAX yes or no columns 
#     print('add population level frequency')
#     gnomad_cols = ['gnomAD_AFR_AF', 'gnomAD_AMR_AF', 'gnomAD_ASJ_AF', 'gnomAD_EAS_AF', 'gnomAD_FIN_AF', 'gnomAD_NFE_AF', 'gnomAD_OTH_AF', 'gnomAD_SAS_AF']
#     df_anno[gnomad_cols] = df_anno[gnomad_cols].fillna(0)
#     df_anno['gnomAD_AF_POPMAX'] = df_anno[gnomad_cols].max(axis=1)

#     # caller artifact filters 
#     print('apply filters')
#     df_anno['FILTER'] = df_anno['FILTER'].replace('clustered_events', 'PASS')
#     df_anno['FILTER'] = df_anno['FILTER'].replace('common_variant', 'PASS')

#     # some specific filter flags should be rescued if oncogenic (i.e. EGFR had issues here)
#     print("rescue filter flags if oncogenic")
#     df_anno['FILTER'] = df_anno.apply(
#         lambda row: 'PASS' if row['oncogenic_binary'] == 'YES' and row['FILTER'] in ['triallelic_site', 'clustered_events;triallelic_site', 'clustered_events;homologous_mapping_event'] else row['FILTER'],
#         axis=1
#     )

#     # Artifact Filter
#     print('artifact filter')
#     df_anno['TGL_FILTER_ARTIFACT'] = df_anno['FILTER'].apply(lambda x: 'PASS' if x == 'PASS' else 'Artifact')

#     # ExAC Filter
#     print('exac filter')
#     df_anno['TGL_FILTER_ExAC'] = df_anno.apply(
#         lambda row: 'ExAC_common' if row['ExAC_common'] == 'YES' and row['Matched_Norm_Sample_Barcode'] == 'unmatched' else 'PASS',
#         axis=1
#     )

#     # gnomAD_AF_POPMAX Filter
#     print('population frequency filter')
#     df_anno['TGL_FILTER_gnomAD'] = df_anno.apply(
#         lambda row: 'gnomAD_common' if row['gnomAD_AF_POPMAX'] > 0.001 and row['Matched_Norm_Sample_Barcode'] == 'unmatched' else 'PASS',
#         axis=1
#     )

#     # VAF Filter
#     print('VAF Filter')
#     df_anno['TGL_FILTER_VAF'] = df_anno.apply(
#         lambda row: 'PASS' if row['tumor_vaf'] >= 0.1 or 
#         (row['tumor_vaf'] < 0.1 and row['oncogenic_binary'] == 'YES' and 
#          ((row['Variant_Classification'] in ['In_Frame_Del', 'In_Frame_Ins']) or 
#           row['Variant_Type'] == 'SNP')) 
#         else 'low_VAF', axis=1
#     )

#     # Mark filters 
#     print('Mark filters')
#     df_anno['TGL_FILTER_VERDICT'] = df_anno.apply(
#         lambda row: 'PASS' if row['TGL_FILTER_ARTIFACT'] == 'PASS' and 
#                           row['TGL_FILTER_ExAC'] == 'PASS' and 
#                           row['TGL_FILTER_gnomAD'] == 'PASS' and 
#                           row['TGL_FILTER_VAF'] == 'PASS' 
#                   else ';'.join([row['TGL_FILTER_ARTIFACT'], row['TGL_FILTER_ExAC'], 
#                                  row['TGL_FILTER_gnomAD'], row['TGL_FILTER_VAF']]), axis=1
#     )

#     return df_anno




# def readGep(gepfile):
#     '''
#     Reads in a sample study file and returns a list of sample names or IDs.
#     Parameters
#     ----------
#     - gepfile (str): Path to the sample study file.
#     Returns
#     -------
#     - study_samples (list): List of sample names or IDs.
#     '''
#     study_samples = []
#     with open(gepfile, 'r') as file:
#         study_samples = [line.strip() for line in file]

#     return study_samples


# def preProcRNA(fpkmfile, enscon, genelist=None):
#     '''
#     Preprocesses RNA expression data by merging with gene symbol annotations, optionally subsetting by a gene list.
#     Parameters
#     ----------
#     - fpkmfile (str): Path to the input concatenated FPKM data from RSEM workflow.
#     - enscon (str): Path to a tab-delimited file with ENSEMBLE gene ID and Hugo_Symbol.
#     - genelist (str, optional): Path to a file with list of Hugo Symbols to report in the final results. Default is None.
#     Returns
#     -------
#     - df (pd.DataFrame): Preprocessed DataFrame with gene expression data.
#     '''
#     # check if genelist is None when preProcRNA.py is called by pycbio.py and genelist is omitted
#     if genelist:
#         print('genelist is used during RNA processing')
#     else:
#         print('genelist is not used during RNA processing')
        
#     # read in data
#     gep_data = pd.read_csv(fpkmfile, sep="\t")
#     ens_conv = pd.read_csv(enscon, sep="\t", header=None)

#     # rename columns
#     ens_conv.columns = ["gene_id", "Hugo_Symbol"]

#     # merge in Hugo's, re-order columns, deduplicate
#     df = pd.merge(gep_data, ens_conv, on="gene_id", how="left")
#     df = df.iloc[:, [-1] + list(range(1, df.shape[1] - 1))]
#     df = df[~df.duplicated(subset=[df.columns[0]])]

#     df.set_index("Hugo_Symbol", inplace=True)

#     # subset if gene list is given
#     if genelist is not None:
#         with open(genelist, 'r') as file:
#             keep_genes = [line.strip() for line in file]

#         df = df[df.index.isin(keep_genes)]
    
#     # return the data frame
#     return df


# def compZ(df):
#     '''
#     Computes row-wise z-scores for a given DataFrame and fills NaN values with zero.
#     Parameters
#     ----------
#     - df (pd.DataFrame): Input DataFrame with gene expression data.
#     Returns
#     -------
#     - df_zscore (pd.DataFrame): DataFrame with z-scores for each gene.
#     '''
#     # scale row-wise
#     df_zscore = df.apply(lambda x: (x - x.mean()) / x.std(), axis=1)
#     df_zscore = df_zscore.fillna(0)

#     return df_zscore




def create_input_directories(outdir, mapfile, merge_maf, merge_seg, merge_fus, merge_gep):
    '''
    (str, str, str, str, str, str) -> None
    
    Create sub-directories in outdir for each type of file listed in map file if at least 1 such file exists.
    Also create sub-directories if data from a previous impor folder exists and need to be merged even if no such data 
    exist in the current map file.
        
    Parameters
    ----------
    - outdir (str): Path to the output directory where mafdir, sgedir, fusdir and gepdir folders are located 
    - mapfile (str): Mapping file (map.csv) that contains paths to maf, seg, gep and mavis files
    - merge_maf (str): Path the maf file that need to be merged or empty string  
    - merge_seg (str): Path the sequenza file that need to be merged or empty string
    - merge_fus (str): Path the mavis file that need to be merged or empty string
    - merge_gep (str): Path the rsem file that need to be merged or empty string 
    '''

    # create input directories based on mapping file   
    for i in ['maf', 'seg', 'gep', 'fus']:
        files = extract_files_from_map(mapfile, i)
        if files:
            filedir = os.path.join(outdir, '{0}dir'.format(i))
            os.makedirs(filedir, exist_ok=True)

    # create input directories if data from previous import folder needs to be merged
    data_files = [merge_maf, merge_seg, merge_fus, merge_gep]
    data_dirs = ['mafdir', 'segdir', 'gepdir', 'fusdir']
    for i in range(len(data_files)):
        if data_files[i] and os.path.isfile(data_files[i]):
            filedir = os.path.join(outdir, data_dirs[i])
            os.makedirs(filedir, exist_ok=True)



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
    

def check_genome_version(mapfile, genome, merge_maf=None):
    '''
    (str, str, str | None) -> None
    
    Check if the genome version (hg19 or hg38) found in each MAF file listed in the map file,
    and in the merge maf, if they exist, are the same as the expected genome ref from the config file
        
    Parameters
    ----------
    - mapfile (str): Mapping file (map.csv) that contains paths to maf, seg, gep and mavis files    
    - genome (str): Genome reference from the config file
    - merge_maf (str | None): Path the maf file that need to be merged or None 
    '''
    
    # get mafs listed in mapfile
    mafs = extract_files_from_map(mapfile, 'maf')
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
    
    # check if merging mafs
    if merge_maf:
        infile = open(merge_maf)
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
    



def write_cases(outputfile, project_name, mapfile, data_type, merge_samples = None, discarded_samples = None):
    '''
    (str, str, str, str, str|None, str|None) -> None
    
    Write list of samples profiled for data type 
            
    Parameters
    ----------
    - outputfile (str): Path to the outputfile
    - project_name (str): Field project_name in configuration file
    - mapfile (str): Mapping file (map.csv) that contains paths to maf, seg, gep and mavis files    
    - data_type (str): The type of data considered.
                       Values accepted are: seq, rna, cna, cna_seq, cna_seq_rna, sv
    - merge_samples (list | None): List of case samples for a given data type from a previous import folder that needs to be merged
    - discarded_samples (list | None): List of samples and/or donors to remove
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
    # merge samples from previous import folder if they exist
    samples.extend(merge_samples)
    
    # remove samples
    for i in discarded_samples:
        if i in samples:
            samples.remove(i)
       
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


def update_clinical_sample_header_with_merging_data(merge_sample_clinical_info, header):
    '''
    (dict, list) -> list

    Returns an updated header including the clinical fields from the clinical sample file
    of a previous import folder that needs to be merged
           
    Parameters
    ----------
    - merge_sample_clinical_info (dict): Dictionary with the sample clinical information of a previous import folder
    - header (list): Lists of column names, data_types and priority from the header of the data_clinical_samples.txt file
    '''
    
    for i in merge_sample_clinical_info:
        for j in merge_sample_clinical_info[i]:
            if j not in header[-1]:
                # add column names to header
                for k in [0,1,4]:
                    header[k].append(j.upper())
                # add data type
                header[2].append(merge_sample_clinical_info[i][j]['datatype'])
                # add priority
                header[3].append('1')
    return header



def parse_clinical_patients(append_data, merge_import_folder, clinical_patient = 'data_clinical_patients.txt'):
    '''
    (bool, str, str) -> list
    
    Returns a list with patient clinical information from a previous import folder if 
    needs to be merged with the current patient information
    
    Parameters
    ----------
    
    - append_data (bool): Create an import folder by merging data from an existing import folder if True
    - merge_import_folder (str): Path to the previous import folder in which data should be merged
    - clinical_patient (str): path to the file with patient clinical information
    '''
    
    if append_data:
        # get expected folders in the import folder
        merge_cbiodir, merge_casedir, merge_suppdir, merge_mafdir, merge_segdir, merge_gepdir, merge_fusdir = get_directories(merge_import_folder)
        filepath = os.path.join(merge_cbiodir, clinical_patient)
        if os.path.isfile(filepath):
            infile = open(filepath)
            content = infile.read().rstrip().split('\n')
            infile.close()
            # get rid of the header
            while any(map(lambda x: x.startswith('#'), content)):
                positions = list(map(lambda x: x.startswith('#'), content))
                pos = [i for i in range(len(positions)) if positions[i]]
                if pos:
                    content.pop(pos[0])
            if 'PATIENT_ID' in content[0]:
                content.pop(0)
            for i in range(len(content)):
                content[i] = content[i].split('\t')
    else:
        content = []
    return content
    


def parse_clinical_samples(append_data, merge_import_folder, clinical_sample = 'data_clinical_samples.txt'):
    '''
    (str, str, str) -> dict
    
    Returns a dictionary with sample clinical information extracted from the 
    clinical sample file of a previous import folder which data needs to be merged
        
    Parameters
    ----------
    - append_data (bool): Create an import folder by merging data from an existing import folder if True
    - merge_import_folder (str): Path to the previous import folder in which data should be merged
    - clinical_sample (str): path to the file with sample clinical information
    '''
    
    # create a dict to store clinical info {'patient;sample': 'field': value}
    D = {}
    if append_data:
        # get expected folders in the import folder
        merge_cbiodir, merge_casedir, merge_suppdir, merge_mafdir, merge_segdir, merge_gepdir, merge_fusdir = get_directories(merge_import_folder)
        filepath = os.path.join(merge_cbiodir, clinical_sample)
        if os.path.isfile(filepath):
            infile = open(filepath)
            content = infile.read().rstrip().split('\n')
            infile.close()

            datatype = content[2].split('\t')
            header = content[4].split('\t')
    
            for i in content[5:]:
                i = i.split('\t')
                patient, sample = i[0], i[1]
                ID = patient + ';' + sample
                D[ID] = {}
                for j in range(len(i)):
                    if j >= 2:
                        D[ID][header[j]] = {'value': i[j], 'datatype': datatype[j]}
                
    return D



def write_patient_minimal_clinical_information(outputfile, mapfile, centre, merge_patient_clinical_info=None, discarded_donors=None):
    '''
    (str, str, str, list | None, list | None) -> None
    
    Write clinical files with minimal clinical information
    
    Parameters
    ----------    
    - outputfile (str): Path to the outputfile
    - mapfile (str): Mapping file (map.csv) that contains paths to maf, seg, gep and mavis files    
    - centre (str): Genomic centre (eg TGL, OICR)
    - merge_patient_clinical_info (list | None): List with patient clinical information to be merged
    - discarded_donors (list | None): List of donors to exclude
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
    # add clinical information to be merged if it exists
    U.extend(merge_patient_clinical_info)
    

    # remove donors 
    if discarded_donors:
        remove = [i for i in U if i[0] in discarded_donors]
        for i in remove:
            U.remove(i)
              
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



def write_sample_minimal_clinical_information(outputfile, mapfile, centre, sample_info = None, merge_sample_clinical_info = None, discarded_samples=None):
    '''
    (str, str, str, str, dict | None, dict | None, list | None) -> None
    
    Write clinical files with minimal clinical information
    
    Parameters
    ----------    
    - outputfile (str): Path to the outputfile
    - mapfile (str): Mapping file (map.csv) that contains paths to maf, seg, gep and mavis files    
    - centre (str): Genomic centre (eg TGL, OICR)
    - sample_info (dict | None): Dictionary with patient and sample information. Populates the sample clinical file
    - merge_sample_clinical_info (dict, | None): Dictionary with clinical sample information to be merged 
    - discarded_samples (list | None): List of samples to exclude
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

    # update header with data type to be merged
    if merge_sample_clinical_info:
        T = update_clinical_sample_header_with_merging_data(merge_sample_clinical_info, T)
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
    if merge_sample_clinical_info:
        for ID in merge_sample_clinical_info:
            data[ID] = [ID.split(';')[0], ID.split(';')[1]] + ['' for j in range((len(T[0]) - 2))]
            # add values to clinical fields            
            for field in merge_sample_clinical_info[ID]:
                data[ID][positions[field]] = merge_sample_clinical_info[ID][field]['value']
    if sample_info:
        for ID in sample_info:
            # add values to clinical fields
            for field in sample_info[ID]:
                data[ID][positions[field]] = sample_info[ID][field]      

    
    for ID in data:
        assert len(data[ID]) == len(T[0])
    
    # remove samples
    if discarded_samples:
        remove = []
        for i in discarded_samples:
            for j in data:
                if i in j:
                    remove.append(j)
        for i in remove:
            del data[i]
    
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



def parse_clinical_oncokb(append_data, merge_import_folder, clinical_oncokb = 'oncokb_clinical_info.txt'):
    '''
    (bool, str, str) -> list
    
    Returns a list of samples extracted from the clinical oncoKB file from a previous import folder if it exists
        
    Parameters
    ----------
    - append_data (bool): Create an import folder by merging data from an existing import folder if True
    - merge_import_folder (str): Path to the previous import folder in which data should be merged
    - clinical_oncokb (str): Path to the clinical file for oncoKB annotation
    '''
    
    L = []
    
    if append_data:
        # get expected folders in the import folder
        merge_cbiodir, merge_casedir, merge_suppdir, merge_mafdir, merge_segdir, merge_gepdir, merge_fusdir = get_directories(merge_import_folder)
        filepath = os.path.join(merge_suppdir, clinical_oncokb)
        if os.path.isfile(filepath):
            infile = open(filepath)
            infile.readline()
            L = infile.read().rstrip().split('\n')
            for i in range(len(L)):
                L[i] = L[i].split('\t')[0]
            infile.close()
            
    return L

   

def write_clinical_oncokb(outputfile, mapfile, cancer_code, merge_clinical_oncokb=None, discarded_samples=None):
    '''
    (str, str, str, list | None, list | None) -> None
    
    Write clinical file for oncokb-annotator
    
    Parameters
    ----------    
    - outputfile (str): Path to the outputfile
    - mapfile (str): Mapping file (map.csv) that contains paths to maf, seg, gep and mavis files    
    - cancer_code (str): Cancer code from OncoTree
    - merge_clinical_oncokb (list, None): List of samples from the clinical oncokb file of a previous import folder that needs to be merged
    - discarded_samples (list | None): List of samples to exclude
    '''
       
    # get library and cancer type from the mapfile
    infile = open(mapfile)    
    content = infile.read().rstrip().split('\n')
    infile.close()
    # create a list of samples
    L = [i.split(',')[1].strip() for i in content]
    
    # add samples to merge if they exist
    if merge_clinical_oncokb:
        L.extend(merge_clinical_oncokb)
    
    # remove samples
    if discarded_samples:
        remove = [i for i in L if i in discarded_samples]
        for i in remove:
            L.remove(i)
    
    if L:
        newfile = open(outputfile, 'w')
        header = ['SAMPLE_ID', 'ONCOTREE_CODE']
        newfile.write('\t'.join(header) + '\n')
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
           os.symlink(files[i], target)
    else:
        print('Cannot link {0} files. No files exist in mapping file {1}'.format(data_type, mapfile))



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
            
    
def concatenate_seg_files(segdir, outputfile, merge_seg=None):
    '''
    (str, str, str | None) -> None
    
    Concatenate seg files located in segdir into outputfile
    
    Parameters
    ----------
    - segdir (str): Directory with seg files
    - outputfile (str): Path to the concatenated seg file
    - merge_seg (str | None): Path the sequenza file that need to be merged if it exists
    '''

    # make a list of seg files
    segfiles = [os.path.join(segdir, i) for i in os.listdir(segdir) if '.seg' in i]
    # get the header of the seg file
    if segfiles:
        infile = open(segfiles[0])
        header = infile.readline()
        infile.close()
    elif merge_seg:
        infile = open(merge_seg)
        header = infile.readline()
        infile.close()
    
    newfile = open(outputfile, 'w')
    newfile.write(header)
    
    # concatenate sequencza files from the map file
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
    
    # add sequenza files from the previous import folder if it exists
    if merge_seg:
        infile = open(merge_seg)
        #skip header
        infile.readline()
        merge_data = infile.read().rstrip()
        infile.close()
        newfile.write(merge_data + '\n')
    
    newfile.close()        


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
    


def select_fusion_file_for_header(fusfiles):
    '''
    (list) -> str
    
    Returns the first fusion file with calls from the list of fusion files or 
    the first file if there no calls in each of the files
    '''
    
    if fusfiles:
        L = []
        for i in fusfiles:
            infile = open(i)
            content = infile.read().strip().split('\n')
            infile.close()
            if len(content) > 1:
                L.append(i)  
                break
        if L:
            return L[0]
        else:
            return fusfiles[0]
    else:
        return []


def get_fusfiles_header(fusfiles, merge_fus=None):
    '''
    (list) -> dict
    
    Returns a dictionary with the input data type WT, WG or both  for each sample
    
    Parameters
    ----------
    - fusfiles (list): List of fusion files
    - merge_fus (str): Path to the concantenated fusion fileto be merged if it exists
    '''
    
    D = {}
    
    for i in fusfiles:
        # check that fusion file has data
        if check_fusion_data(i):
            sample = get_sample_from_filename(i)
            D[sample] = []
            infile = open(i)
            header = infile.readline().rstrip().split('\t')
            infile.close()
            for j in range(len(header)):
                if 'WT.' in header[j]:
                    D[sample].append('WT')
                elif 'WG.' in header[j]:
                    D[sample].append('WG')
                D[sample].sort()
    
    if merge_fus:
        infile = open(merge_fus)
        header = infile.readline().rstrip().split('\t')
        samples = infile.read().rstrip().split('\n')
        if samples:
            for i in range(len(samples)):
                samples[i] = samples[i].split('\t')[0]
        L = []
        for j in range(len(header)):
            if 'WT.' in header[j] or header[j] == 'WT':
                L.append('WT')
            elif 'WG.' in header[j] or header[j] == 'WG':
                L.append('WG')
        for i in samples:
            D[i] = L
            
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




def extract_merged_fusion(merge_fusion):
    '''
    (str) -> list
    
    Returns a list of dictionaies each representing a line of data in merge_fusion
    annotated with the column header 
    
    Parameters
    ----------
    - merge_fusion (str): Path to the merged fusion file
    '''
    
    L = []
    
    infile = open(merge_fusion)
    header = infile.readline().rstrip().split('\t')
    for line in infile:
        line = line.rstrip()
        if line:
            line = line.split('\t')
            d = {}
            for i in range(len(header)):
                if 'WT.' in  header[i] or header[i] == 'WT':
                    d['WT'] = line[i]
                elif 'WG.' in header[i] or header[i] == 'WG':
                    d['WG'] = line[i]
                else:
                    d[header[i]] = line[i]
            L.append(d)        
    infile.close()
    return L


def concatenate_fusion_files(fusdir, outputfile, merge_fus=None):
    '''
    (str, str, str | None) -> None
    
    Concatenates fusion files located in fusdir into outputfile. Also adds fusion data
    from previous import folder that needs to be merged, if it exists
    
    Parameters
    ----------
    - fusdir (str): Directory with fusion files
    - outputfile (str): Path to the concatenated seg file
    - merge_fus (str | None): Path the mavis file that need to be merged if it exists
    '''

    # make a list of fusion files
    fusfiles = [os.path.join(fusdir, i) for i in os.listdir(fusdir) if '.fus' in i]
    
    # determine if the headers have WT, WG or both
    header_types = get_fusfiles_header(fusfiles, merge_fus)
      
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
    
    for file in fusfiles:
        # check the content of fusion file
        if check_fusion_data(file):
            # extract sample name from file name
            sample = get_sample_from_filename(file)
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
        
    # add data from pervious import folder to merge if it exists
    if merge_fus:
        # extract data from the fusion file
        merge_data = extract_merged_fusion(merge_fus)
        for d in merge_data:
            newline = []
            for i in range(len(header)):
                if header[i] in d:
                    newline.append(d[header[i]])
                else:
                    newline.append('')
            newfile.write('\t'.join(newline) + '\n')
        
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


def parse_merge_gep(merge_gep):
    '''
    (str) -> dict
    
    Returns a dictionary with fpkm for each gene and sample in the merge_gep file
    
    Parameters
    ----------
    - merge_gep (str): Path to the fpkm file that need to be merged or empty string
    '''

    D = {}
    infile = open(merge_gep)
    header = infile.readline().rstrip().split('\t')
    samples = header[1:]
    for line in infile:
        line = line.rstrip()
        if line:
            line = line.split('\t')
            gene = line[0]
            fpkm = line[1:]
            assert len(samples) == len(fpkm)
            for i in range(len(samples)):
                if samples[i] not in D:
                    D[samples[i]] = {}
                D[samples[i]][gene] = fpkm[i]
    infile.close()
    
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
    

def concatenate_fpkm_from_gep_files(gepdir, outputfile, merge_gep):
    '''
    (str, str, str) -> None
    
    Write fpkm for all genes and samples with rna data to outputfile
    
    Parameters
    ----------
    - gepdir (str): Directory with rsem files
    - outputfile (str): Path to the outputfile with fpkm for each sample and gene
    - merge_gep (str): Path to the fpkm file that need to be merged or empty string
    '''
    
    # collect all fpkm for each gene and sample 
    fpkm = collect_fpkm(gepdir) if gepdir else {}
    merge_fpkm = parse_merge_gep(merge_gep) if merge_gep else {}
    # merge both dictionaries
    D = {}
    if fpkm:
        D.update(fpkm)
    if merge_fpkm:
        D.update(merge_fpkm)
        
    # write fpkm to outputfile
    write_fpkm_to_file(D, outputfile)


   
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
    


def concatenate_maf_files(mafdir, outputfile, merge_maf=None):
    '''
    (str, str, str | None) -> None
    
    Concatenates all the gzipped maf files in mafdir into a text file outputfile
    with column Tumor_Sample_Barcode replaced by the sample name.
    Also merge the concatenated mafs from a previous import folder if merge_maf is provided
    
    Parameters
    ----------
    - mafdir (str): Directory containing the maf files
    - outputfile (str): Path to the concatenated maf file (unzipped)
    - merge_maf (str | None): Path the maf file that need to be merged or None
    '''
     
    # make a list of maf files
    maffiles = [os.path.join(mafdir, i) for i in os.listdir(mafdir) if '.maf.gz' in i]
    # get the header of the maf file
    if maffiles:
        header = get_maf_header(maffiles[0])
    elif merge_maf:
        infile = open(merge_maf)
        header = infile.readline().rstrip()
        infile.close()
    
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
        infile.close() 
    
    # add concatenated mafs from previous import folder
    if merge_maf:
        infile = open(merge_maf)
        # skip header
        infile.readline()
        # get all the mutations
        content = infile.read()
        infile.close()
        newfile.write(content)
       
    newfile.close()

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
    - keep_variants(bool): Keep variants with missing gnomAD_AF values when Matched_Norm_Sample_Barcode is unmatched if True 
    - removedfile (str):
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
                          "5'Flank", "Splie_Region", "Targeted_Region"]
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
    #must match with the keys
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
    - removedfile (str):
    '''
     # count the number of mutations before and after filtering
    total, kept = 0, 0
    
    # open files
    newfile = open(outputfile, 'w')
    infile = open(maffile)
    
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


    return total, kept 
              


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
    expected_resources = ['token', 'enscon_hg38', 'enscon_hg19', 'entcon', 'genebed_hg38', 'genebed_hg19', 'genelist', 'oncolist']
    missing_resources = [i for i in expected_resources if i not in list(config['Resources'].keys())]
    if missing_resources:
        raise ValueError('ERROR. Missing resources: {0}'.format(', '.join(missing_resources)))
    invalid_resource_files = [i for i in ['token', 'enscon_hg38', 'enscon_hg19', 'entcon', 'genebed_hg38', 'genebed_hg19', 'genelist', 'oncolist'] if config['Resources'][i] and os.path.isfile(config['Resources'][i]) == False]
    if invalid_resource_files:
        raise ValueError('ERROR. Provide valid path for {0}'.format(', '.join(invalid_resource_files)))
    
    # check options
    missing_options = [i for i in ['mapfile', 'outdir', 'study', 'center', 'cancer_code', 'keep_variants'] if config['Options'][i] is None]
    if missing_options:
        raise ValueError('ERROR. Provide values for {0} in the config'.format(', '.join(missing_options)))
    # check map file
    if os.path.isfile(config['Options']['mapfile']) == False:
        raise ValueError('ERROR: Provide valid path to mapfile in config')
    # check that bbolean filter parameter is provided
    try:
        keep_variants = config['Options'].getboolean('keep_variants')
    except:
        raise ValueError('ERROR. {0} is not a boolean. Use true or false'.format(keep_variants))
        
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
    
    resources = ['token', 'enscon_hg38', 'enscon_hg19', 'entcon', 'genebed_hg38', 'genebed_hg19', 'genelist', 'oncolist']
    L = [config['Resources'][i] for i in resources]
    token, enscon_hg38, enscon_hg19, entcon, genebed_hg38, genebed_hg19, genelist, oncolist = L
    return token, enscon_hg38, enscon_hg19, entcon, genebed_hg38, genebed_hg19, genelist, oncolist
    

def extract_options_from_config(config):
    '''
    (configparser.ConfigParser) -> (str, str, str)
    
    Returns the variables listed in the Options section of the config
    
    Parameters
    ----------
    - config (configparser.ConfigParser): Config file parsed with configparser
    '''    
    
    options = ['mapfile', 'outdir', 'project_name', 'description', 'study', 'center', 'cancer_code', 'genome']
    L = [config['Options'][i] for i in options]
    # add boolean filter
    L.append(config['Options'].getboolean('keep_variants'))
    mapfile, outdir, project_name, description, study, center, cancer_code, genome, keep_variants = L
    return mapfile, outdir, project_name, description, study, center, cancer_code, genome, keep_variants


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


def check_input_mafs(mapfile, merge_maf = None):
    '''
    (str, str | None) -> None
    
    Exits if the input maf files have different headers
        
    Parameters
    ----------
    - mapfile (str): Mapping file (map.csv) that contains paths to maf, seg, gep and mavis files    
    - merge_maf (str | None): Path the maf file that need to be merged or None 
    '''

    # make a list of input maf files 
    mafs = extract_files_from_map(mapfile, 'maf')
    
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
    
    if merge_maf:
        infile = open(merge_maf)
        for line in infile:
            if 'Hugo_Symbol' in line:
                headers.append(line.rstrip())
                break
        infile.close()
               
    # check if multiple headers
    if len(list(set(headers))) > 1:
        sys.exit('Input MAF files have different headers')
    


def get_directories(outdir):
    '''
    (str) -> list
    
    Retrieves the list of expected directories in the outdir directory
    
    Parameters
    ----------
    - outdir (str): Path to the output directory containing the import folder and associated directories
    '''
    
    # expected folders in the output directory
    cbiodir = os.path.join(outdir, 'cbioportal_import_data')
    casedir = os.path.join(cbiodir, 'case_lists')
    suppdir = os.path.join(outdir, 'supplementary_data')
    mafdir =  os.path.join(outdir, 'mafdir')
    segdir = os.path.join(outdir, 'segdir')
    gepdir = os.path.join(outdir, 'gepdir')
    fusdir = os.path.join(outdir, 'fusdir')
    
    # replace with empty string if folder doesn't exist
    # folder is created only if there are input files
    expected = [cbiodir, casedir, suppdir, mafdir, segdir, gepdir, fusdir]
    for i in range(len(expected)):
        if os.path.isdir(expected[i]) == False:
            expected[i] == ''
    
    return expected


def get_concatenated_input_files(filedir, filetype):
    '''
    (str, str) -> str
    
    Returns the path of the concatenated file in the filedir or the empty string
    if the file does not exist
    
    Parameters
    ----------
    - filedir (str): Path to the directory where the concatenated file is expected
    - filetype (str): Type of the expected file. Accepted values are maf, fus, seg and gep    
    '''
    
    if filedir:
        if filetype == 'maf':
            inputfile = os.path.join(filedir, 'all_mutations.maf.txt')
        elif filetype == 'fus':
            inputfile = os.path.join(filedir, 'input.fus.txt')
        elif filetype == 'gep':
            inputfile = os.path.join(filedir, 'input.fpkm.txt')
        elif filetype == 'seg':
            inputfile = os.path.join(filedir, 'input.seg.txt')
 
    if os.path.isfile(inputfile):
        return inputfile
    else:
        return '' 
    


def check_merging_option(append_data, merge_import_folder, outdir):    
    '''
    (bool, str, str)
    
    Raise an Error if the combitions of parameters is not not valid
    
    Parameters
    ----------
    - append_data (bool): Create an import folder by merging data from an existing import folder if True
    - merge_import_folder (str): Path to the previous import folder in which data should be merged
    - outdir (str): Output directory from the config file where the import folder is created
    '''
    
    # check that outdir is different than import folder dir
    # check that correct options are used when merging data
    if append_data:
        if outdir == merge_import_folder:
            raise ValueError('''WARNING. The output directory is the same as the importer folder you want to merge.
                          Review the config file and/or provide a different path to the merging import folder''') 
        if merge_import_folder is None or os.path.isdir(merge_import_folder) == False:
            raise ValueError('Provide the path to the previous import folder')
    else:
        if merge_import_folder:
            raise ValueError('Paths to the previous import folder can only be used when merging data')    


def get_data_to_merge(append_data, merge_import_folder):
    '''
    (bool, str) -> (str, str, str, str)
    
    Returns paths to maf, sequenza, mavis and rsem concatenated data from a previous import folder
    if they exist or an empty string, if data from a previous import folder is to be merged to 
    create a new import folder
    
    Parameters
    ----------
    - append_data (bool): Create an import folder by merging data from an existing import folder if True
    - merge_import_folder (str): Path to the previous import folder in which data should be merged
    '''
    
    # get raw data to merge (concatenated files before any filtering and processing) 
    if append_data:
        # get expected folders in the import folder
        merge_cbiodir, merge_casedir, merge_suppdir, merge_mafdir, merge_segdir, merge_gepdir, merge_fusdir = get_directories(merge_import_folder)
        # get maf input file
        merge_maf = get_concatenated_input_files(merge_mafdir, 'maf')
        # get seg input file
        merge_seg = get_concatenated_input_files(merge_segdir, 'seg')
        # get fus input file
        merge_fus = get_concatenated_input_files(merge_fusdir, 'fus')
        # get gep input file
        merge_gep = get_concatenated_input_files(merge_gepdir, 'gep')
    else:
        merge_maf, merge_seg, merge_fus, merge_gep = '', '', '', ''
    
    return merge_maf, merge_seg, merge_fus, merge_gep




def get_samples_merge(append_data, merge_import_folder, case_file):
    '''
    (bool, str, str) -> list
    
    Returns a list of case samples extracted from the case_file of a previous import folder
    
    - append_data (bool): Create an import folder by merging data from an existing import folder if True
    - merge_import_folder (str): Path to the previous import folder in which data should be merged
    - case_file (str): Name of the case file.The type of data considered.
                       Values accepted are: cases_sequenced.txt, cases_rna_seq_mrna.txt,
                       cases_cna.txt, cases_cnaseq.txt, cases_3way_complete.txt, cases_sv.txt
    '''
    
    samples = []
    
    # get raw data to merge (concatenated files before any filtering and processing) 
    if append_data:
        # get expected folders in the import folder
        merge_cbiodir, merge_casedir, merge_suppdir, merge_mafdir, merge_segdir, merge_gepdir, merge_fusdir = get_directories(merge_import_folder)
        filepath = os.path.join(merge_casedir, case_file)
        if os.path.isfile(filepath):
            infile = open(filepath)
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


def check_exclude_option(removed_samples):
    '''
    (str) -> None
    
    Raise an Error if the file path to the removed file is incorrect when this option is used
        
    Parameters
    ----------
    - removed_samples (str): Path to the files with samples to be excluded
    '''
    
    if removed_samples:
        if os.path.isfile(removed_samples) == False:
            raise ValueError('Provide the path to the file with excluded samples')
    

def parse_excluded_samples(file):
    '''
    (str) -> list
    
    Returns a list of samples to be excluded from the import folder
    
    Parameters
    ----------
    - file (str): Path to the file with the samples to exclude
    '''
    
    infile = open(file)
    samples = infile.read().strip().split('\n')
    infile.close()
    
    return samples
    
       
def map_donors_samples(mapfile, append_data, merge_import_folder):
    '''
    (str, bool, str) -> dict
    
    Returns a dictionary mapping all donors to their samples 
    from the current mapping file and a previous import folder if data needs to be merged
    
    Parameters
    ----------
    - mapfile (str): Path to the current mapping file
    - append_data (bool): Merge data if True
    - merge_import_folder (str): Path to a previous import folder
    '''
    
    # read the map file
    infile = open(mapfile)
    content = infile.read().rstrip().split('\n')
    infile.close()
    
    
    # get the cases from the previous import folder if merging folders
    merge_samples = [get_samples_merge(append_data, merge_import_folder, i) for i in   
                     ['cases_sequenced.txt', 'cases_rna_seq_mrna.txt', 'cases_cna.txt',
                      'cases_cnaseq.txt', 'cases_3way_complete.txt', 'cases_sv.txt']] 
    
    # collect all samples
    samples = [i.split(',')[1] for i in content]
    for i in merge_samples:
        samples.extend(i)
    samples = list(set(samples))
        
    # map donors to samples
    D = {}
    for  i in samples:
        donor = i.split('_')[:2]
        if donor not in D:
            D[donor] = []
        D[donor].append(i)
        
    return D

    
def exclude_donors(mapfile, append_data, merge_import_folder, removed_samples):
    '''
    (str, bool, str, str) -> list
    
    Returns a list of donors that need to be excluded because all their corresponding 
    samples are in the removed samples file
    
    Parameters
    ----------
    - mapfile (str): Path to the current mapping file
    - append_data (bool): Merge data if True
    - merge_import_folder (str): Path to a previous import folder
    - removed_samples (str): Path to the files with excluded samples
    '''
    
    D = map_donors_samples(mapfile, append_data, merge_import_folder)
    excluded = parse_excluded_samples(removed_samples)
    
    for sample in excluded:
        for donor in D:
            if sample in D[donor]:
                D[donor].removed(sample)
    # get all the donors that do not have any samples
    excluded_donors = [i for i in D if len(D[i]) == 0]
        
    return excluded_donors




def remove_samples_from_maf(maffile, discarded_samples):
    '''
    (str, list) -> int
    
    Delete mutations from the maf file for samples in discarded samples
    and returns the number of removed mutations
    
    Parameters
    ----------
    - maffile (str): Path to the concatenated maf file
    - discarded_samples (list): List of samples to exclude
    '''

    # create a list of data excluding the discarded samples
    content = []
    # count the number of mutations removed
    removed = 0
    
    infile = open(maffile)
    # add header
    content.append(infile.readline().rstrip())
    # parse file content and skip lines containing the discarded samples
    for line in infile:
        line = line.rstrip()
        if line:
            line = line.split('\t')
            # keep mutations for samples not in discarded samples
            if line[15] not in discarded_samples:
                content.append('\t'.join(line))
            else:
                removed += 1
    infile.close()
    
    if removed:
        # open maffile for writing
        newfile = open(maffile, 'w')
        for i in content:
           newfile.write(i + '\n')
        newfile.close()        
    
    return removed


def remove_samples_from_fusion(fusfile, discarded_samples):
    '''
    (str, list) -> int
    
    Delete fusions from the fusion file for samples in discarded samples
    and returns the number of removed fusions
    
    Parameters
    ----------
    - fusfile (str): Path to the concatenated fusion file
    - discarded_samples (list): List of samples to exclude
    '''
    
    # create a list of data excluding the discarded samples
    content = []
    # count the number of fusions removed
    removed = 0
    
    infile = open(fusfile)
    # add header
    content.append(infile.readline().rstrip())
    # parse file content and skip lines containing the discarded samples
    for line in infile:
        line = line.rstrip()
        if line:
            line = line.split('\t')
            # keep fusions for samples not in discarded samples
            if line[0] not in discarded_samples:
                content.append('\t'.join(line))
            else:
                removed += 1
    infile.close()
    
    if removed:
        # open fusfile for writing
        newfile = open(fusfile, 'w')
        for i in content:
            newfile.write(i + '\n')
        newfile.close()        
    
    return removed




def find_position_discarded_sample(gepfile, discarded_samples):
    '''
    (str, list) -> list
    
    Returns a list of indices corresponding to the positions of each discarded sample
    in the header of the gep file in decreasing order
        
    Parameters
    ----------
    - gepfile (str): Path to the concatenated gep file
    - discarded_samples (list): List of samples to exclude
    '''
    
    L = []
    
    infile = open(gepfile)
    header = infile.readline().rstrip().split('\t')
    infile.close()
    for i in discarded_samples:
        if i in header:
            L.append(header.index(i))    
    # sort the list in decreasing order 
    L = list(reversed(sorted(L)))
    
    return L


def remove_samples_from_gepfile(gepfile, discarded_samples):
    '''
    (str, list) -> int
    
    Delete fusions from the fusion file for samples in discarded samples
    and returns the number of removed fusions
    
    Parameters
    ----------
    - fusfile (str): Path to the concatenated fusion file
    - discarded_samples (list): List of samples to exclude
    '''

    # get the indices of the discarded samples in the header of the gepfile
    excluded_positions = find_position_discarded_sample(gepfile, discarded_samples)

    # create a list of data excluding the discarded samples
    content = []
    # count the number of genes/transcripts for which data is removed
    removed = 0    
    
    if excluded_positions:
        infile = open(gepfile)
        header = infile.readline().rstrip().split('\t')
        # remove samples from header
        for i in excluded_positions:
            header = header[:i] + header[i+1:]
        # add header
        content.append('\t'.join(header))
    
        # parse file content and skip lines containing the discarded samples
        for line in infile:
            line = line.rstrip()
            if line:
                line = line.split('\t')
                # exclude each column (sample) in line
                for i in excluded_positions:
                    line = line[:i] + line[i+1:]
                content.append('\t'.join(line))
                removed += 1
        infile.close()    
            
        # open fusfile for writing
        newfile = open(gepfile, 'w')
        for i in content:
            newfile.write(i + '\n')
        newfile.close()        
    
    return removed



def remove_samples_from_segfile(segfile, discarded_samples):
    '''
    (str, list) -> int
    
    Delete events from the concatenated seg file for samples in discarded samples
    and returns the number of removed events
    
    Parameters
    ----------
    - segfile (str): Path to the concatenated segfile
    - discarded_samples (list): List of samples to exclude
    '''
    
    # create a list of data excluding the discarded samples
    content = []
    # count the number of fusions removed
    removed = 0
    
    infile = open(segfile)
    # add header
    content.append(infile.readline().rstrip())
    # parse file content and skip lines containing the discarded samples
    for line in infile:
        line = line.rstrip()
        if line:
            line = line.split('\t')
            # keep fusions for samples not in discarded samples
            if line[0] not in discarded_samples:
                content.append('\t'.join(line))
            else:
                removed += 1
    infile.close()
    
    if removed:
        # open fusfile for writing
        newfile = open(segfile, 'w')
        for i in content:
            newfile.write(i + '\n')
        newfile.close()        
    
    return removed





def make_import_folder(args):
    '''
    (list) -> None
    
    Generate folder cbioportal_import_folder with metadata and processed data files
    to upload to cBioPortal
        
    Parameters
    ----------
    - config (str): Path to the config file
    - clinical (str): Path to the sample clinical file
    - append_data (bool): Create an import folder by merging data from an existing import folder if True
    - merge_import_folder (str): Path to the previous import folder in which data should be merged
    - removed_samples (str): Path to the files with samples to be excluded
    '''
    
    # parse config file
    config = configparser.ConfigParser(allow_no_value=True)
    config.read(args.config)
    # check config
    check_configuration(config)
    print('read and checked config')
    
    # extract variables from config
    token, enscon_hg38, enscon_hg19, entcon, genebed_hg38, genebed_hg19, genelist, oncolist = extract_resources_from_config(config)
    mapfile, outdir, project_name, description, study, center, cancer_code, genome, keep_variants = extract_options_from_config(config)
    gain, amplification, heterozygous_deletion, homozygous_deletion, minfusionreads = extract_parameters_from_config(config)
    depth_filter, alt_freq_filter, gnomAD_AF_filter, tglpipe, filter_variants, filter_indels = extract_filters_from_config(config)
    print('extracted variables from config')
    
    # check that outdir is different than import folder dir
    # check that correct options are used when merging data
    check_merging_option(args.append_data, args.merge_import_folder, outdir)    
    
    # get raw data to merge (concatenated files before any filtering and processing) 
    # if they exist, or empty strings
    merge_maf, merge_seg, merge_fus, merge_gep = get_data_to_merge(args.append_data, args.merge_import_folder)
    
    # check that input maf files, if any, have the same format and the same header
    check_input_mafs(mapfile, merge_maf)
    
    # check genome version in the maf files, if provided
    check_genome_version(mapfile, genome, merge_maf)
    # get genome specific variables
    if genome == 'hg38':
        enscon, genebed = enscon_hg38, genebed_hg38
    elif genome == 'hg19':
        enscon, genebed = enscon_hg19, genebed_hg19
    print('determined genome specific resources')
        
    # check that cancer type is correctly defined
    check_cancer_type(cancer_code)
    print('checked cancer code')

    # create output directory and output sub-folders. remove output directory if it exists
    cbiodir, casedir, suppdir = create_output_directories(outdir)
    print('created output directories')

    # create input directories for each file type from map file    
    create_input_directories(outdir, mapfile, merge_maf, merge_seg, merge_fus, merge_gep)
    print('created input directories')
        
    # get the samples to exclude
    check_exclude_option(args.removed_samples)
    if args.removed_samples:
        excluded_samples = parse_excluded_samples(args.removed_samples)
        excluded_donors =  exclude_donors(mapfile, args.append_data, args.merge_import_folder, args.removed_samples)
        print('Excluding {0} donors and {1} samples'.format(len(excluded_donors), len(excluded_samples)))
    else:
        print('All specified donors and samples are used.')
        excluded_samples, excluded_donors = [], []
        
    # write meta study and clinical files
    write_meta_study(os.path.join(cbiodir, 'meta_study.txt') , study, project_name, description, genome, cancer_code)
    write_meta_clinical(cbiodir, project_name, 'sample')
    write_meta_clinical(cbiodir, project_name, 'patient')
    print('wrote study and clinical metadata')    
    
    # write cases
    merge_samples = [get_samples_merge(args.append_data, args.merge_import_folder, i) for i in   
                     ['cases_sequenced.txt', 'cases_rna_seq_mrna.txt', 'cases_cna.txt',
                      'cases_cnaseq.txt', 'cases_3way_complete.txt', 'cases_sv.txt']] 
                        
    write_cases(os.path.join(casedir, 'cases_sequenced.txt'), project_name, mapfile, 'seq',
                merge_samples = merge_samples[0], discarded_samples = excluded_samples)
    write_cases(os.path.join(casedir, 'cases_rna_seq_mrna.txt'), project_name, mapfile, 'rna',
                merge_samples = merge_samples[1], discarded_samples = excluded_samples)
    write_cases(os.path.join(casedir, 'cases_cna.txt'), project_name, mapfile, 'cna',
                merge_samples = merge_samples[2], discarded_samples = excluded_samples)
    write_cases(os.path.join(casedir, 'cases_cnaseq.txt'), project_name, mapfile, 'cna_seq',
                merge_samples = merge_samples[3], discarded_samples = excluded_samples)
    write_cases(os.path.join(casedir, 'cases_3way_complete.txt'), project_name, mapfile, 'cna_seq_rna',
                merge_samples = merge_samples[4], discarded_samples = excluded_samples)
    write_cases(os.path.join(casedir, 'cases_sv.txt'), project_name, mapfile, 'sv',
                merge_samples = merge_samples[5], discarded_samples = excluded_samples) 
    print('wrote cases')

    
    # write patient clinical information
    clinical_outputfile = os.path.join(cbiodir, 'data_clinical_patients.txt')
    # get clinical file from previous import folder if merging
    merge_patient_clinical_info = parse_clinical_patients(args.append_data, args.merge_import_folder, 'data_clinical_patients.txt')
    write_patient_minimal_clinical_information(clinical_outputfile, mapfile, center,
                                               merge_patient_clinical_info=merge_patient_clinical_info,
                                               discarded_donors=excluded_donors)
        
    # write sample clinical information
    # get clinical sample file from previous import folder if merging
    merge_sample_clinical_info =  parse_clinical_samples(args.append_data, args.merge_import_folder, 'data_clinical_samples.txt')
    # get the user defined sample clinical information
    clinical_info = get_clinical_data(args.clinical) if args.clinical else None
    clinical_outputfile = os.path.join(cbiodir, 'data_clinical_samples.txt')
    write_sample_minimal_clinical_information(clinical_outputfile, mapfile, center,
                                              sample_info = clinical_info,
                                              merge_sample_clinical_info = merge_sample_clinical_info,
                                              discarded_samples = excluded_samples)
    print('wrote clinical information')
    
       
    # write clinical input file for oncokb-annotator
    clinical_oncokb = os.path.join(suppdir, 'oncokb_clinical_info.txt')
    merge_clinical_oncokb = parse_clinical_oncokb(args.append_data, args.merge_import_folder, 'oncokb_clinical_info.txt')
    write_clinical_oncokb(clinical_oncokb, mapfile, cancer_code,
                          merge_clinical_oncokb=merge_clinical_oncokb,
                          discarded_samples=excluded_samples)
    print('wrote oncoKb clinical information') 
    
    
    # link files
    for i in ['maf', 'seg', 'gep', 'fus']:
        link_files(outdir, mapfile, i)
    print('linked files to input directories')
     
        
    # concatenate input files
    # concatenate mafs    
    mafs = extract_files_from_map(mapfile, 'maf')
    mafdir = os.path.join(outdir, 'mafdir')
    if mafs or merge_maf:
        assert os.path.isdir(mafdir)
        mutation_file = os.path.join(mafdir, 'all_mutations.maf.txt')
        # concatenate all maf files into a plain text file
        concatenate_maf_files(mafdir, mutation_file, merge_maf)
        print('concatenated maf files')
    else:
        mutation_file = ''
    
    # check if seg dir and seg files exist
    segs = extract_files_from_map(mapfile, 'seg')
    segdir = os.path.join(outdir, 'segdir')
    if segs or merge_seg:
        assert os.path.isdir(segdir)
        segfile = os.path.join(segdir, 'input.seg.txt')
        # concatenate seg files if they exist
        concatenate_seg_files(segdir, segfile, merge_seg)
        print('concatenated seg files')
    else:
        segfile = ''

    # check if fusion dir and fusion files exist
    fusions = extract_files_from_map(mapfile, 'fus')
    fusdir = os.path.join(outdir, 'fusdir')
    if fusions or merge_fus:
        assert os.path.isdir(fusdir)
        fusfile = os.path.join(fusdir, 'input.fus.txt')
        # concatenate fusion files, if they exist
        concatenate_fusion_files(fusdir, fusfile, merge_fus)
        print('concatenated fusion files')
        # check if fusion data
        if check_fusion_data(fusfile) == False:
            fusfile = ''
            print('no fusion data in concatenated fusion file')
        
    else:
        fusfile = ''
    
    # extract and concatenate fpkm from gep files
    geps = extract_files_from_map(mapfile, 'gep')
    gepdir = os.path.join(outdir, 'gepdir')
    if geps or merge_gep:
        assert os.path.isdir(gepdir)
        gepfile = os.path.join(gepdir, 'input.fpkm.txt')
        # write fpkm to file
        concatenate_fpkm_from_gep_files(gepdir, gepfile, merge_gep)
        print('concatenated fpkm from gep files')
    else:
        gepfile = ''
    
    
    # remove samples from the data files
    if excluded_samples:
        if mutation_file:
            removed_mutations = remove_samples_from_maf(mutation_file, excluded_samples)
            print('Removed {0} mutations for {1} samples in mutation file {2}'.format(removed_mutations, len(excluded_samples), mutation_file))
        if segfile:
            removed_events = remove_samples_from_segfile(segfile, excluded_samples)
            print('Removed {0} events for {1} samples in seg file {2}'.format(removed_events, len(excluded_samples), segfile))
        if fusfile:
            removed_fusions = remove_samples_from_fusion(fusfile, excluded_samples)
            print('Removed {0} fusions for {1} samples in fusion file {2}'.format(removed_fusions, len(excluded_samples), fusfile))
        if gepfile:
            removed_expression = remove_samples_from_gepfile(gepfile, excluded_samples)
            print('Removed {0} genes for {1} samples in gep file {2}'.format(removed_expression, len(excluded_samples), gepfile))
    
           
    # filter maf files and write metadata
    # define maffile, output of MafAnnotator
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
        oncokb_token = get_token(token)
    
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
            df_cbio_filt = pd.read_csv(args.maffile, sep="\t", header=0)
            df_snv = df_cbio_filt[df_cbio_filt['Variant_Type'] == 'SNP']
            df_cbio_filt.to_csv(os.path.join(cbiodir, 'data_mutations_extended.txt'), sep="\t", index=False, na_rep='NA')

        # write metadata file
        write_metadata(os.path.join(cbiodir, 'meta_mutations_extended.txt'), project_name, 'maf', genome)
        print('wrote mutations metadata')
    
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
        list_gep_samples(gepdir, gep_study_file)
        # generate expression data files
        print('Processing RNASEQ data from {0}'.format(gepfile))
        # get list of samples in study
        study_samples = readGep(gep_study_file)
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
    g_parser.add_argument('--append', dest='append_data', action='store_true', help='Create an import folder by merging data from an existing import folder if True')
    g_parser.add_argument('-mid', '--MergeImportDirectory', dest='merge_import_folder', help='Path to the previous import folder in which data should be merged')
    g_parser.add_argument('-rs', '--RemovedSamples', dest='removed_samples', help='Path to the files with samples to be excluded')
    g_parser.set_defaults(func=make_import_folder)
    
    # generate import folder
    m_parser = subparsers.add_parser('map', help="Generate map file")
    m_parser.add_argument('-o', '--outdir', dest='outdir', help='Path to the output directory', required = True)
    m_parser.add_argument('-i', '--input_data', dest='input_data', help='Path to the json file with input data', required = True)
    m_parser.add_argument('-g', '--gamma', dest='gamma', default='500', help='Gamma value of the sequenza output. Default is 500')
    m_parser.set_defaults(func=generate_mapfile)
    
    # get arguments from the command line
    args = parser.parse_args()
    # pass the args to the default function
    args.func(args)
