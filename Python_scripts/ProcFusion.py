import os
import argparse
import pandas as pd 

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
            lambda x: max(map(int, x.split(';'))) if isinstance(x, str) and ';' in x else int(x)
        )
    
    for column in columns:
        df[column] = pd.to_numeric(df[column], errors='coerce')  

    return df


def preProcFus(datafile, readfilt, entrfile):
    '''
    Preprocesses the input data and formats it for CBioPortal output. 
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
    df_cbio['Entrez_Gene_Id'] = df_cbio['Entrez_Gene_Id'].astype(int)

    # remove rows where gene is not known (this still keeps the side of the gene which is known)
    df_cbio = df_cbio.dropna()

    return df_cbio

if __name__ == "__main__":
    # create a parser for command line arguments 
    parser = argparse.ArgumentParser(
        description='Process gene fusion data from mavis and generate outputs for CBioPortal'
    )

    parser.add_argument(
        'fusfile',
        type=str,
        help='Path to the input- A concatenated mavis fusion data file'
    )
    parser.add_argument(
        'entcon',
        type=str,
        help='Path to a tab-delimted file with Entrez gene IDs and Hugo Symbol'
    )
    parser.add_argument(
        'minfusionreads',
        type=int,
        help='Minimum number of reads for fusion calls'
    )
    parser.add_argument(
        'outdir',
        type=str,
        help='Path to the output directory'
    )

    args = parser.parse_args()

    # make subdirectories 
    cbiodir = os.path.join(args.outdir, "cbioportal_import_data")
    suppdir = os.path.join(args.outdir, "supplementary_data")
    os.makedirs(cbiodir, exist_ok=True)
    os.makedirs(suppdir, exist_ok=True)

    # function returns list of 3 objects ### TO WRITE
    print(f"Reading fusion data from {args.fusfile}")
    fusion_cbio = preProcFus(args.fusfile, args.minfusionreads, args.entcon)

    # write FUS files
    print(f"writing fus file")
    fusion_cbio.to_csv(os.path.join(cbiodir, "data_fusions.txt"), sep="\t", index=False)