'''
To do list:
1. search within star_gene_counts.tsv for LINC00094, MIR1270,
calculate mean value for each stage folder.
'''

import os
from pathlib import Path
import pandas as pd
from glob import glob
import matplotlib


def main():
    ### Input parameters. ###
    project_name = 'TCGA_SKCM'
    output_folder = Path('C:/Repositories/Melanoma_TCGA/analysis/')
    project_folder = output_folder / project_name
    files_re_list = ['*.rna_seq.augmented_star_gene_counts*.tsv']
    columns_to_check = ['gene_id', 'unstranded', 'stranded_first', 'stranded_second']
    columns_we_want = ['gene_name', 'tpm_unstranded', 'fpkm_unstranded', 'fpkm_uq_unstranded']

#    test_file = Path('C:/Repositories/Melanoma_TCGA/data/dupe_files_check/TCGA-D3-A1QA.rna_seq.augmented_star_gene_counts.tsv')
#    df = pd.read_table(test_file, usecols = columns_we_want,
#                       index_col = 0, skiprows = 1)
#    print(df.shape)
#    print(df.head())

    #print(df[df.index == 'MIR1270'])

    for file_re in files_re_list:
        # Path for glob search.
        search_path = project_folder / '*' / '*' / file_re
        # List of matching files in project folder.
        files_list = glob(str(search_path))
        print(f'"{len(files_list)}" matching files found for "{search_path}"')
        # Turn all "path" into "parent folder name".
        case_folder_list = [Path(f).parents[1].name + '/' +\
                            Path(f).parents[0].name\
                            for f in files_list]
        #print(len(case_folder_list))
        # Count the occurances of folder name.
        folder_count = pd.Series(case_folder_list).value_counts()
        target_folder = folder_count[folder_count > 1].index
        #print(target_folder)
        # Iterate "target_folder" for files.
        for folder in target_folder:
            #print(folder)
            target_path = project_folder / folder / file_re
            dupe_files = glob(str(target_path))
            file_num = 0    # Give number to each file.
            file_compare = {}   # Dictionary for file compare.
            for file in dupe_files:
                df = pd.read_table(file, usecols = columns_to_check,
                                   index_col = 0, skiprows = 1, nrows = 4)
                #print(df)
                #print(df.iloc[0, 0])
                # Add file_num and value to dictionary for comparison.
                file_compare[file_num] = df.iloc[0, 0]
                file_num += 1
            #print(file_compare)
            # Find the "key"(s) other than the smallest value.
            file_num_drop = [k for k in file_compare\
                             if file_compare[k] != min(file_compare.values())]
            #print(file_num_drop)
            # The file(s) to drop.
            for i in file_num_drop:
                print(f'Removing duplicate file from list: "{dupe_files[i]}"')
                files_list.remove(dupe_files[i])
                #print(len(files_list))

        df_new = pd.DataFrame() # DataFrame to store values.
        for file in files_list[0:2]:
            print(file)
            df = pd.read_table(file, usecols = columns_we_want,
                       index_col = 0, skiprows = 1)
            df.index.name = None    # Remove index name.
            print(df.shape)
            df_temp = df[df.index == 'MIR1270'].copy()
            #print(df_temp)
            df_temp['case_id'] = Path(file).parents[0].name
            df_temp['tumor_stage'] = Path(file).parents[1].name
            #print(df_temp)

            ### Next step: Merge df_temp to df_new (now it's empty)
            df_new.merge(df_temp)
            print(df_new)




if __name__ == '__main__':
    main()