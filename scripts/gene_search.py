'''
To do list:
1. search within star_gene_counts.tsv for LINC00094, MIR1270,
calculate mean value for each stage folder.
'''

import os
from pathlib import Path
import pandas as pd
import numpy as np
from glob import glob
import matplotlib.pyplot as plt
import seaborn as sns
from file_sorter import create_folder, replace_special_chars


def main():
    ### Input parameters. ###
    project_name = 'TCGA_SKCM'
    output_folder = Path('C:/Repositories/Melanoma_TCGA/analysis/')
    project_folder = output_folder / project_name
    files_re_list = ['*.rna_seq.augmented_star_gene_counts*.tsv']
    columns_to_check = ['gene_id', 'unstranded', 'stranded_first', 'stranded_second']
    columns_we_want = ['gene_name', 'tpm_unstranded', 'fpkm_unstranded', 'fpkm_uq_unstranded']
    target_gene = ['MIR1270']


    ### Create project analysis folder if not already exist. ###
    project_name = replace_special_chars(project_name)
    analysis_name = f'{project_name}_analysis'
    analysis_folder = create_folder(analysis_name, output_folder,
                                   verbose = True)


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
        files_removed = []  # A list to store removed duplicate file path.
    #    for folder in target_folder:
    #        #print(folder)
    #        target_path = project_folder / folder / file_re
    #        dupe_files = glob(str(target_path))
    #        file_num = 0    # Give number to each file.
    #        file_compare = {}   # Dictionary for file compare.
    #        for file in dupe_files:
    #            df = pd.read_table(file, usecols = columns_to_check,
    #                               index_col = 0, skiprows = 1, nrows = 4)
    #            #print(df)
    #            #print(df.iloc[0, 0])
    #            # Add file_num and value to dictionary for comparison.
    #            file_compare[file_num] = df.iloc[0, 0]
    #            file_num += 1
    #        #print(file_compare)
    #        # Find the "key"(s) other than the smallest value.
    #        file_num_drop = [k for k in file_compare\
    #                         if file_compare[k] != min(file_compare.values())]
    #        #print(file_num_drop)
    #        # The file(s) to drop.
    #        for i in file_num_drop:
    #            # Add file path to files_removed.
    #            files_removed.append(dupe_files[i])
    #            print(f'Removing duplicate file from list: "{dupe_files[i]}"')
    #            files_list.remove(dupe_files[i])
    #            #print(len(files_list))
        for gene in target_gene:
            # Create gene folder in analysis_folder if not already exist.
            gene_name = replace_special_chars(gene)
            gene_folder = create_folder(gene_name, analysis_folder,
                                           verbose = True)
        #    print(f'Performing gene search for: {gene}...')
        #    df_gene = gene_search(files_list, gene = gene,
        #                          usecols = columns_we_want,
        #                          index_col = 0, skiprows = 1)
        #    print(df_gene.shape)
        #    ### Save to files. ###
        #    # Removed files.
        #    print(f'Writing file: "{gene_folder}/files_removed.txt"...')
        #    with open(f'{gene_folder}/files_removed.txt', 'w') as f:
        #        f.write(f'>Duplicate files removed for "{file_re}"\n')
        #        for line in files_removed:
        #            f.write(f'{line}\n')
        #    # Gene dataframe.
        #    print(f'Writing file: "{gene_folder}/{gene}.txt"...')
        #    df_gene.to_csv(f'{gene_folder}/{gene}.txt', sep = '\t', index = False)

            # Load file to skip above operation.
            file_path = f'{gene_folder}/{gene}.txt'
            df_gene = pd.read_csv(file_path, sep='\t', header=0)
            #print(df_gene.shape)

            tumor_stage = df_gene.groupby('tumor_stage')
            # Count for each stage.
            #print(tumor_stage.size())
            print(tumor_stage.median())
            print(tumor_stage.mean())
            ### Save to file. ###
            tumor_stage.median().to_csv(f'{gene_folder}/stage_median.txt', sep = '\t')
            tumor_stage.mean().to_csv(f'{gene_folder}/stage_mean.txt', sep = '\t')

            # Remove 0
            df_gene.iloc[:, 0:3] = df_gene.iloc[:, 0:3].replace(0, np.nan)
            #print(df_gene.head())

    #        # Generate histogram for column.
    #        # Set plot style
    #        plt.style.use('ggplot')
    #        for column in columns_we_want[1:4]:
    #            print(f'Histogram for column: "{column}"')
    #            histogram(df_gene, column)
            
            tumor_stage = df_gene.groupby('tumor_stage')
            # Count for each stage.
            #print(tumor_stage.size())
            print(tumor_stage.median())
            print(tumor_stage.mean())
            ### Save to file. ###
            tumor_stage.median().to_csv(f'{gene_folder}/stage_median_rm0.txt', sep = '\t')
            tumor_stage.mean().to_csv(f'{gene_folder}/stage_mean_rm0.txt', sep = '\t')           

        
def gene_search(files_list, gene, usecols, 
                index_col = 0, skiprows = 0):
    df_new = pd.DataFrame() # DataFrame to store values.
    for file in files_list:
        #print((f'Working on file: "{file}"'))
        df = pd.read_table(file, usecols = usecols,
                   index_col = index_col, skiprows = skiprows)
        df.index.name = None    # Remove index name.
        #print(df.shape)
        df_temp = df[df.index == gene].copy()
        #print(df_temp)
        df_temp['case_id'] = Path(file).parents[0].name
        df_temp['tumor_stage'] = Path(file).parents[1].name
        #print(df_temp)
        # Merge df_temp to df_new.
        df_new = pd.concat([df_new, df_temp], axis=0, join='outer')
    return df_new


def histogram(df, x):
    fig, ax1 = plt.subplots()
    sns.histplot(df, x=x, ax=ax1)
    plt.tight_layout() ### Rescale the fig size to fit the data
    plt.show()


if __name__ == '__main__':
    main()