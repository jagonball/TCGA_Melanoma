'''
Cox's proportional hazard model
'''
import sys
from pathlib import Path
import pandas as pd
import numpy as np
from lifelines import CoxPHFitter
from file_sorter import create_folder


def main():
    ### Input parameters. ###
    project_name = 'TCGA_SKCM'
    output_folder = Path('C:/Repositories/Melanoma_TCGA/analysis/')
    project_folder = output_folder / project_name
    target_gene = 'BRD3OS'
    gene_folder = output_folder / 'TCGA_SKCM_analysis' / target_gene
    
    # Main data with all cases' information.
    main_data = Path('C:/Repositories/Melanoma_TCGA/data/clinical_patient_skcm.txt')
    # The columns we want.
    columns_main = ['bcr_patient_barcode', 'vital_status', 
                    'last_contact_days_to', 'death_days_to']
    gene_file = gene_folder / 'merge_BRD3OS.txt'

    ### Create survival folder if not already exist. ###
    survival_folder = create_folder('survival', gene_folder,
                                   verbose = True)

    # Read main_data into DataFrame.
    main_df = pd.read_table(main_data, sep = '\t', header = 0,
                            index_col = 'bcr_patient_barcode',
                            usecols = columns_main, skiprows = [1,2])
    #print(main_df.shape)
    #print(main_df.head())
    new_df = main_df[main_df['vital_status'] == 'Alive'].copy()
    #print(new_df.shape)
    # Remove rows without proper value.
    new_df = new_df[~new_df['last_contact_days_to'].str.startswith('[')]
    #print(new_df.shape)
    #print(new_df.head())
    temp_df = main_df[main_df['vital_status'] == 'Dead'].copy()
    #print(temp_df.shape)
    # Remove rows without proper value.
    temp_df = temp_df[~temp_df['death_days_to'].str.startswith('[')]
    #print(temp_df.shape)

    # Rename column to 'days' then concat DataFrame.
    new_df = new_df.rename(columns={'last_contact_days_to': 'days'})
    temp_df = temp_df.rename(columns={'death_days_to': 'days'})
    new_df = pd.concat([new_df, temp_df], axis=0, join='outer')
    # Drop the unwanted columns.
    new_df = new_df.drop(columns=['last_contact_days_to', 'death_days_to'])
    new_df = new_df.replace({'Alive': 0, 'Dead': 1})
    print(new_df.shape)
    #print(new_df.head(-5))
    # Save file to check.
    #new_df.to_csv(survival_folder / 'new_df.txt',
    #              sep = '\t')


    ### Annotate gene information. ###
    df_gene = pd.read_csv(gene_file, sep='\t', header=0,
                          index_col = 0)
    print(df_gene.shape)
    #print(df_gene.head())
    # Calculate mean and standard deviation, then setup a cutoff value.
    #gene_mean = df_gene['tpm_unstranded'].mean()
    #gene_std = df_gene['tpm_unstranded'].std()
    #cutoff = gene_mean + gene_std
    # Use quantile to set cutoff.
    set_quantile = 0.9
    gene_quantile = df_gene.quantile(set_quantile)
    cutoff = gene_quantile[0]
    print(f'cutoff is set to: {cutoff}')
    # Create new column with cutoff value.
    # if >= cutoff, vlaue = 1; otherwise value = 0.
    df_gene['BRD3OS'] = np.where(df_gene['tpm_unstranded'] >= cutoff, 1, 0)
    #df_gene['BRD3OS'] = df_gene['tpm_unstranded'].copy()
    #print(df_gene.head(-5))
    # Save file to check.
    #df_gene.to_csv(survival_folder / 'BRD3OS_cutoff.txt',
    #               sep = '\t')

    # Concat DataFrame with index.
    new_df = pd.concat([new_df, df_gene], axis=1, join='inner')
    # Drop the unwanted columns.
    new_df = new_df.drop(columns=['tpm_unstranded', 'fpkm_unstranded',
                                  'fpkm_uq_unstranded', #'tumor_stage',
                                  'stage_T', 'stage_N', 'stage_M'])
    print(new_df.shape)
    #print(new_df.head(-5))
    # Dictionary for stage rename.
    tumor_stage = {'Stage_0': 0,
                   'Stage_I': 1,
                   'Stage_II': 2,
                   'Stage_III': 3,
                   'Stage_IV': 4,
                   '_Not_Available_': np.nan,
                   'I_II_NOS': np.nan}
    # All stages: '[Not Available]', 'T0', 'T1~4', 'Tis', 'TX'.
    stage_T = {'T0': 0,
               'T1': 1,
               'T2': 2,
               'T3': 3,
               'T4': 4,
               '[Not Available]': np.nan,
               'TX': np.nan,
               'Tis': np.nan}

    # All stages: '[Not Available]', 'N0', 'N1~3', 'NX'.
    stage_N = {'N0': 0,
               'N1': 1,
               'N2': 2,
               'N3': 3,
               '[Not Available]': np.nan,
               'NX': np.nan}

    # All stages: '[Not Available]', 'M0', 'M1'.
    stage_M = {'M0': 0,
               'M1': 1,
               '[Not Available]': np.nan}

    new_df['tumor_stage'] = new_df['tumor_stage'].replace(tumor_stage)
    #new_df['stage_T'] = new_df['stage_T'].replace(stage_T)
    #new_df['stage_N'] = new_df['stage_N'].replace(stage_N)
    #new_df['stage_M'] = new_df['stage_M'].replace(stage_M)
    #print(new_df.head(-5))
    new_df = new_df.dropna()
    print(new_df.shape)
    #print(new_df.head(-5))

    # Save file to check.
    #new_df.to_csv(survival_folder / 'new_df_final.txt',
    #              sep = '\t')

    # Change dtype for CoxPHFitter.
    new_df['vital_status'] = new_df['vital_status'].astype('int64')
    new_df['days'] = new_df['days'].astype('int64')
    print(new_df['vital_status'].dtype)
    print(new_df['days'].dtype)
    print(new_df['BRD3OS'].dtype)

    ### Perform CoxPHFitter. ###
    cph = CoxPHFitter()
    cph.fit(new_df, duration_col='days', event_col='vital_status')
    
    cph.print_summary()  # access the individual results using cph.summary.
    # Save the output to file.
    original_stdout = sys.stdout # Save a reference to the original standard output.
    with open(survival_folder / 'CoxPHFitter.txt', 'w') as f:
        sys.stdout = f # Change the standard output to the file we created.
        cph.print_summary()
        sys.stdout = original_stdout # Reset the standard output to its original value.


if __name__ == '__main__':
    main()