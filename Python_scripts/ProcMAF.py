import os
import argparse
import pandas as pd
import gc 

def addVAFtoMAF(maf_df, alt_col, dep_col, vaf_header):
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
    data = pd.read_csv(datafile, sep="\t")

    #print('--- doing some formatting ---')

    # add vaf columns 
    #print('add tumor_vaf')
    data = addVAFtoMAF(data, 't_alt_count', 't_depth', 'tumor_vaf')
    #print('add normal_vaf')
    data = addVAFtoMAF(data, 'n_alt_count', 'n_depth', 'normal_vaf')

    # clear memory (important when the mafs are huge with millions and millions of lines)
    df_anno = data
    del data
    gc.collect()

    # add oncogenic yes or no columns 
    #print('add oncogenic status')
    df_anno['oncogenic_binary'] = df_anno['oncogenic'].apply(lambda x: 'YES' if x in ['Oncogenic', 'Likely Oncogenic'] else 'NO')

    # add common_variant yes or no columns
    df_anno['ExAC_common'] = df_anno['FILTER'].apply(lambda x: 'YES' if 'common_variant' in x else 'NO')

    # add POPMAX yes or no columns 
    #print('add population level frequency')
    gnomad_cols = ['gnomAD_AFR_AF', 'gnomAD_AMR_AF', 'gnomAD_ASJ_AF', 'gnomAD_EAS_AF', 'gnomAD_FIN_AF', 'gnomAD_NFE_AF', 'gnomAD_OTH_AF', 'gnomAD_SAS_AF']
    df_anno[gnomad_cols] = df_anno[gnomad_cols].fillna(0)
    df_anno['gnomAD_AF_POPMAX'] = df_anno[gnomad_cols].max(axis=1)

    # caller artifact filters 
    #print('apply filters')
    df_anno['FILTER'] = df_anno['FILTER'].replace('clustered_events', 'PASS')
    df_anno['FILTER'] = df_anno['FILTER'].replace('common_variant', 'PASS')
    df_anno['FILTER'] = df_anno['FILTER'].apply(lambda x: 'PASS')

    # some specific filter flags should be rescued if oncogenic (i.e. EGFR had issues here)
    #print("rescue filter flags if oncogenic")
    df_anno['FILTER'] = df_anno.apply(
        lambda row: 'PASS' if row['oncogenic_binary'] == 'YES' and row['FILTER'] in ['triallelic_site', 'clustered_events;triallelic_site', 'clustered_events;homologous_mapping_event'] else row['FILTER'],
        axis=1
    )

    # Artifact Filter
    #print('artifact filter')
    df_anno['TGL_FILTER_ARTIFACT'] = df_anno['FILTER'].apply(lambda x: 'PASS' if x == 'PASS' else 'Artifact')

    # ExAC Filter
    #print('exac filter')
    df_anno['TGL_FILTER_ExAC'] = df_anno.apply(
        lambda row: 'ExAC_common' if row['ExAC_common'] == 'YES' and row['Matched_Norm_Sample_Barcode'] == 'unmatched' else 'PASS',
        axis=1
    )

    # gnomAD_AF_POPMAX Filter
    #print('population frequency filter')
    df_anno['TGL_FILTER_gnomAD'] = df_anno.apply(
        lambda row: 'gnomAD_common' if row['gnomAD_AF_POPMAX'] > 0.001 and row['Matched_Norm_Sample_Barcode'] == 'unmatched' else 'PASS',
        axis=1
    )

    # VAF Filter
    #print('VAF Filter')
    df_anno['TGL_FILTER_VAF'] = df_anno.apply(
        lambda row: 'PASS' if row['tumor_vaf'] >= 0.1 or 
        (row['tumor_vaf'] < 0.1 and row['oncogenic_binary'] == 'YES' and 
         ((row['Variant_Classification'] in ['In_Frame_Del', 'In_Frame_Ins']) or 
          row['Variant_Type'] == 'SNP')) 
        else 'low_VAF', axis=1
    )

    # Mark filters 
    #print('Mark filters')
    df_anno['TGL_FILTER_VERDICT'] = df_anno.apply(
        lambda row: 'PASS' if row['TGL_FILTER_ARTIFACT'] == 'PASS' and 
                          row['TGL_FILTER_ExAC'] == 'PASS' and 
                          row['TGL_FILTER_gnomAD'] == 'PASS' and 
                          row['TGL_FILTER_VAF'] == 'PASS' 
                  else ';'.join([row['TGL_FILTER_ARTIFACT'], row['TGL_FILTER_ExAC'], 
                                 row['TGL_FILTER_gnomAD'], row['TGL_FILTER_VAF']]), axis=1
    )

    return df_anno


if __name__ == "__main__":
    # create a parser for command line arguments 
    parser = argparse.ArgumentParser(
        description='Process and filter data from the Variant Effect Predictor (VEP) workflow and generate outputs for CBioPortal'
    )

    parser.add_argument(
        'maffile',
        type=str,
        help='Path to the input concatenated MAFs (Mutation Annotation Format) from VEP and annotated with OncoKb'
    )
    parser.add_argument(
        'tglpipe',
        type=bool,
        help='filters variants according to TGL specifications if True'
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


    print(f'Processing Mutation data from {args.maffile}')
    
    # only do the filtering steps if tglpipe is set to TRUE
    if args.tglpipe:
        print('tglpipe is set to true, filtering data according to tgl specifications')
        df_cbio_anno = procVEP(args.maffile)
        df_cbio_filt = df_cbio_anno[df_cbio_anno['TGL_FILTER_VERDICT'] == 'PASS']

        # get snvs for dcsigs
        df_snv = df_cbio_filt[df_cbio_filt['Variant_Type'] == 'SNP']

        # for cbioportal input
        df_cbio_filt.to_csv(os.path.join(cbiodir, 'data_mutations_extended.txt'), sep="\t", index=False)
        
        # unfiltered data
        df_cbio_anno.to_csv(os.path.join(suppdir, 'unfiltered_data_mutations_extended.txt'), sep="\t", index=False)
    
    else:
        df_cbio_filt = pd.read_csv(args.maffile, sep="\t", header=0)
        df_snv = df_cbio_filt[df_cbio_filt['Variant_Type'] == 'SNP']
        df_cbio_filt.to_csv(os.path.join(cbiodir, 'data_mutations_extended.txt'), sep="\t", index=False)