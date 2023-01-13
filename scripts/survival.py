'''
Cox’s proportional hazard model
'''
from pathlib import Path
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.datasets import load_rossi


def main():
    ### Input parameters. ###
    project_name = 'TCGA_SKCM'
    output_folder = Path('C:/Repositories/Melanoma_TCGA/analysis/')
    project_folder = output_folder / project_name
    
    # Main data with all cases' information.
    main_data = Path('C:/Repositories/Melanoma_TCGA/data/clinical_patient_skcm.txt')
    columns_main = ['vital_status', 'last_contact_days_to', 'death_days_to']

    ### EXAMPLE ###
    #rossi = load_rossi()
    ##print(rossi.shape)
    ##print(rossi.head())
    #print(rossi['week'].dtype)
    #print(rossi['arrest'].dtype)
#
    #cph = CoxPHFitter()
    #cph.fit(rossi, duration_col='week', event_col='arrest')
    #
    #cph.print_summary()  # access the individual results using cph.summary

    # Read main_data into DataFrame.
    main_df = pd.read_table(main_data, sep = '\t', header = 0, 
                            usecols = columns_main, skiprows = [1,2])
    #print(main_df.shape)
    #print(main_df.head())
    new_df = main_df[main_df['vital_status'] == 'Alive'].copy()
    #print(new_df.shape)
    new_df = new_df[~new_df['last_contact_days_to'].str.startswith('[')]
    #print(new_df.shape)
    #print(new_df.head())
    temp_df = main_df[main_df['vital_status'] == 'Dead'].copy()
    #print(temp_df.shape)
    temp_df = temp_df[~temp_df['death_days_to'].str.startswith('[')]
    #print(temp_df.shape)
    new_df = new_df.rename(columns={'last_contact_days_to': 'days'})
    temp_df = temp_df.rename(columns={'death_days_to': 'days'})
    new_df = pd.concat([new_df, temp_df], axis=0, join='outer')
    new_df = new_df.drop(columns=['last_contact_days_to', 'death_days_to'])
    new_df = new_df.replace({'Alive': 0, 'Dead': 1})
    #print(new_df.shape)
    #print(new_df.head(-5))
    new_df['vital_status'] = new_df['vital_status'].astype('int64')
    new_df['days'] = new_df['days'].astype('int64')
    print(new_df['vital_status'].dtype)
    print(new_df['days'].dtype)


    cph = CoxPHFitter()
    cph.fit(new_df, duration_col='days', event_col='vital_status')
    
    cph.print_summary()  # access the individual results using cph.summary


if __name__ == '__main__':
    main()