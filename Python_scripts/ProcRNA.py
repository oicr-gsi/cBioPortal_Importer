import os
import argparse
import pandas as pd 
from scipy.stats import zscore

def readGep(gepfile):
    study_samples = []
    with open(gepfile, 'r') as file:
        study_samples = [line.strip() for line in file]

    return study_samples

# preprocess function
def preProcRNA(fpkmfile, enscon, genelist=None):
    # check if genelist is None when preProcRNA.py is called by pycbio.py and genelist is omitted
    if genelist == "None":
        genelist = None 
        print('genelist is not used during RNA processing')
    else:
        print('genelist is used during RNA processing')
    
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
    df.drop(df.columns[0], axis=1, inplace=True)

    # subset if gene list is given
    if genelist is not None:
        with open(genelist, 'r') as file:
            keep_genes = [line.strip() for line in file]

        df = df[df.index.isin(keep_genes)]
    
    # return the data frame
    return df

# simple zscore function
def compZ(df):
    # scale row-wise
    df_zscore = df.apply(zscore, axis=1)
    df_zscore = df_zscore.fillna(0)

    return df_zscore


if __name__ == "__main__":
    # create a parser for command line arguments 
    parser = argparse.ArgumentParser(
        description='Filter and process rsem data and generate outputs for CBioPortal'
    )

    parser.add_argument(
        'fpkmfile',
        type=str,
        help='Path to the input concantenated fpkm data from rsem workflow'
    )
    parser.add_argument(
        'enscon',
        type=str,
        help='Path to a tab-delimted file with ENSEMBLE gene ID and Hugo_Symbol'
    )
    parser.add_argument(
        'gepfile',
        type=str,
        help='Path to a list of samples in the study'
    )
    parser.add_argument(
        'genelist',
        type=str,
        nargs='?',
        default=None,
        help='(Optional) Path to a file with list of Hugo Symbols to report in the final results.'
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

    print(f'Processing RNASEQ data from {args.fpkmfile}')

    # get list of samples in study
    study_samples = readGep(args.gepfile)

    # preprocess the full data frame
    df = preProcRNA(args.fpkmfile, args.enscon, args.genelist)
    print('getting STUDY-level data')

    # subset data to STUDY level data for cbioportal
    col = [col for col in study_samples if col in df.columns]
    df_study = df[col]

    # write the raw STUDY data
    df_study.to_csv(os.path.join(cbiodir, "data_expression.txt"), sep="\t", header=True, index=True)

    # z-scores STUDY
    df_zscore = compZ(df_study)
    df_zscore.to_csv(os.path.join(cbiodir, "data_expression_zscores.txt"), sep="\t", header=True, index=True)