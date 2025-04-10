import os
import argparse
import pandas as pd 
from rpy2 import robjects 
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr

pandas2ri.activate()
cntools = importr('CNTools')

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
    if genelist == "None":
        genelist = None 
        print('genelist is not used during CNA processing')
    else:
        print('genelist is used during CNA processing')
    
    # read oncogenes
    oncogenes = pd.read_csv(oncolist, sep='\t')

    # small fix segmentation data
    segData = pd.read_csv(segfile, sep='\t')
    segData['chrom'] = segData['chrom'].str.replace('chr', '')

    # convert pandas dataframe to R dataframe 
    segData_r = pandas2ri.py2rpy(segData)

    # set thresholds
    print('setting thresholds')
    gain = float(gain)
    amp = float(amp)
    htz = float(htz)
    hmz = float(hmz)

    # get the gene info
    print('getting gene info')
    geneInfo = pd.read_csv(genebed, sep='\t')

    # convert pandas dataframe to R dataframe 
    geneInfo_r = pandas2ri.py2rpy(geneInfo)

    # make CN matrix gene level
    print('converting seg')
    cnseg = cntools.CNSeg(segData_r)
    print('get segmentation by gene')
    rdByGene = cntools.getRS(cnseg, by='gene', imput=False, XY=False, geneMap=geneInfo_r, what='median', mapChrom='chrom', mapStart='start', mapEnd='end')
    print('get reduced segmentation data')
    reducedseg_df = cntools.rs(rdByGene)

    # convert data from R back to Py
    reducedseg = pandas2ri.rpy2py(reducedseg_df)

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
        keep_genes = [line.strip() for line in open(genelist)]
        df_cna = df_cna.loc[df_cna.index.isin(keep_genes)]
        df_cna_thresh = df_cna_thresh.loc[df_cna_thresh.index.isin(keep_genes)]
    
    return segData, df_cna, df_cna_thresh


if __name__ == "__main__":
    # create a parser for command line arguments 
    parser = argparse.ArgumentParser(
        description='Filter and process CNA data and generate outputs for CBioPortal'
    )

    parser.add_argument(
        'segfile',
        type=str,
        help='Path to the input concantenated segmentation file from sequenza'
    )
    parser.add_argument(
        'genebed',
        type=str,
        help='Path to a tab-delimted bed file which defines the genomic positions of the canonical genes.'
    )
    parser.add_argument(
        'genelist',
        type=str,
        help='Path to a file with list of Hugo Symbols to report in the final results.'
    )
    parser.add_argument(
        'oncolist',
        type=str,
        help='Path to a tab-delimted file with list of cancer genes'
    )
    parser.add_argument(
        'gain',
        type=float,
        help='CNA parameter: Threshold for gains'
    )
    parser.add_argument(
        'amp',
        type=float,
        help='CNA parameter: Threshold for amplification'
    )
    parser.add_argument(
        'htz',
        type=float,
        help='CNA parameter: Threshold for heterozygous deletion'
    )
    parser.add_argument(
        'hmz',
        type=float,
        help='CNA parameter: Threshold for homozygous deletion'
    )
    parser.add_argument(
        'outdir',
        type=str,
        help='Path to the output directory'
    )

    args = parser.parse_args()

    # make subdirectories 
    cbiodir = os.path.join(args.outdir, 'cbioportal_import_data')
    suppdir = os.path.join(args.outdir, 'supplementary_data')
    os.makedirs(cbiodir, exist_ok=True)
    os.makedirs(suppdir, exist_ok=True)

    # function returns list of 3 objects
    print(f'Processing CNA data from {args.segfile}')
    segData, df_cna, df_cna_thresh = preProcCNA(args.segfile, args.genebed, args.gain, args.amp, args.htz, args.hmz, args.oncolist, args.genelist)

    # write cbio files
    print(f'writing seg file')
    segData.to_csv(os.path.join(cbiodir, 'data_segments.txt'), sep='\t', index=False)
    df_cna.to_csv(os.path.join(cbiodir, 'data_log2CNA.txt'), sep='\t', index=True)
    df_cna_thresh.to_csv(os.path.join(cbiodir, 'data_CNA.txt'), sep='\t', index=True)

    # write the truncated data_CNA file (remove genes which are all zero) for oncoKB annotator
    df_CNA = df_cna_thresh.loc[~(df_cna_thresh == 0).all(axis=1)]
    df_CNA.to_csv(os.path.join(suppdir, 'data_CNA_short.txt'), sep='\t', index=True)
