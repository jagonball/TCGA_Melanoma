'''
Search files "star_gene_counts.tsv" for "LINC00094" "BRD3OS" or "MIR1270",
calculate median, mean value for each stage.
'''
from pathlib import Path
import pandas as pd
import numpy as np
from glob import glob
import matplotlib.pyplot as plt
import seaborn as sns
from file_sorter import create_folder, replace_special_chars

__author__ = "Johnathan Lin <jagonball@gmail.com>"
__email__ = "jagonball@gmail.com"

def main():
    ### Input parameters. ###
    project_name = 'TCGA_SKCM'
    output_folder = Path('C:/Repositories/Melanoma_TCGA/analysis/')
    project_folder = output_folder / project_name
    # List of target file's name with regular expression.
    files_re_list = ['*.rna_seq.augmented_star_gene_counts*.tsv']
    # Columns for comparing duplicate files.
    columns_to_check = ['gene_id', 'unstranded',
                        'stranded_first', 'stranded_second']
    columns_we_want = ['gene_name', 'tpm_unstranded',
                       'fpkm_unstranded', 'fpkm_uq_unstranded']
    target_gene = ['MIR1270', 'BRD3OS', 'BRD3']
    # Main data with all cases' information.
    main_data = Path('C:/Repositories/Melanoma_TCGA/data/clinical_patient_skcm.txt')
    # The columns we want.
    columns_main = ['bcr_patient_barcode', 'ajcc_tumor_pathologic_pt', 
                    'ajcc_nodes_pathologic_pn', 'ajcc_metastasis_pathologic_pm']

    ### Create project analysis folder if not already exist. ###
    project_name = replace_special_chars(project_name)
    analysis_name = f'{project_name}_analysis'
    analysis_folder = create_folder(analysis_name, output_folder,
                                   verbose = True)


    # Read main_data into DataFrame.
    df_main = pd.read_table(main_data, sep = '\t', header = 0,
                            index_col = 'bcr_patient_barcode',
                            usecols = columns_main, skiprows = [1,2])
    #print(df_main.shape)
    # Rename columns.
    columns_rename = {'ajcc_tumor_pathologic_pt': 'stage_T',
                      'ajcc_nodes_pathologic_pn': 'stage_N',
                      'ajcc_metastasis_pathologic_pm': 'stage_M'}
    df_main = df_main.rename(columns = columns_rename)


    # Get target files list. 
    for file_re in files_re_list:
        print(f'Searching "{file_re}"...')
        # Path for glob search.
        search_path = project_folder / '*' / '*' / file_re
        # List of matching files in project folder.
        files_list = glob(str(search_path))
        print(f'"{len(files_list)}" matching files found for "{search_path}"')
        # Turn all "path" into "stage/case_id".
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
        for folder in target_folder:
            #print(folder)
            # Path for glob search.
            target_path = project_folder / folder / file_re
            dupe_files = glob(str(target_path))
            file_num = 0    # Give a number to each file.
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
            # Remove file(s) to drop from files_list.
            for f in file_num_drop:
                # Add file path to files_removed.
                files_removed.append(dupe_files[f])
                print(f'Removing duplicate file from list: "{dupe_files[f]}"')
                files_list.remove(dupe_files[f])
                #print(len(files_list))
        

        # Find target_gene within files in files_list.
        for gene in target_gene:
            print(f'Working on: "{gene}"')
            # Create gene folder in analysis_folder if not already exist.
            gene_name = replace_special_chars(gene)
            gene_folder = create_folder(gene_name, analysis_folder,
                                           verbose = True)
            print(f'Performing gene search for: "{gene}"...')
            df_gene = gene_search(files_list, gene = gene,
                                  usecols = columns_we_want,
                                  index_col = 0, skiprows = 1)
            print(df_gene.shape)

            # Annotate T, M, N stage information.
            df_gene = df_gene.set_index('case_id')
            print(df_gene.shape)
            #print(df_gene.head(-5))
            df_gene = pd.concat([df_gene, df_main], axis=1, join='inner')
            print(df_gene.shape)

            ### Save to files. ###
            # Removed files.
            print(f'Writing file: "{gene_folder}/files_removed.txt"...')
            with open(f'{gene_folder}/files_removed.txt', 'w') as f:
                f.write(f'>Duplicate files removed for "{file_re}"\n')
                for line in files_removed:
                    f.write(f'{line}\n')
            # Gene dataframe.
            print(f'Writing file: "{gene_folder}/{gene}.txt"...')
            df_gene.to_csv(f'{gene_folder}/{gene}.txt', sep = '\t')


            ## Load file to skip the above gene search. ##
            file_path = f'{gene_folder}/{gene}.txt'
            df_gene = pd.read_csv(file_path, sep='\t', header=0)
            print(df_gene.shape)
            #print(df_gene.iloc[0,0])
            #print(type(df_gene.iloc[0,0]))


            ### Analysis. ###
            tumor_stage = df_gene.groupby('tumor_stage')
            stage_T = df_gene.groupby('stage_T')
            stage_N = df_gene.groupby('stage_N')
            stage_M = df_gene.groupby('stage_M')
            # Count for each stage.
            #print(tumor_stage.size())
            #print(tumor_stage.median())
            #print(tumor_stage.mean())

            
            ### Save to file. ###
            # For "tumor_stage".
            tumor_stage.size().to_csv(gene_folder /
                                      'stage_count.txt', 
                                      sep = '\t', header = ['cases'])
            tumor_stage.median().to_csv(gene_folder /
                                        'stage_median.txt',
                                        sep = '\t')
            tumor_stage.mean().to_csv(gene_folder /
                                      'stage_mean.txt',
                                      sep = '\t')
            # For "stage_T".
            stage_T.size().to_csv(gene_folder /
                                  'stage_T_count.txt', 
                                  sep = '\t', header = ['cases'])
            stage_T.median().to_csv(gene_folder /
                                    'stage_T_median.txt',
                                    sep = '\t')
            stage_T.mean().to_csv(gene_folder /
                                  'stage_T_mean.txt',
                                  sep = '\t')
            # For "stage_N".
            stage_N.size().to_csv(gene_folder /
                                  'stage_N_count.txt', 
                                  sep = '\t', header = ['cases'])
            stage_N.median().to_csv(gene_folder /
                                    'stage_N_median.txt',
                                    sep = '\t')
            stage_N.mean().to_csv(gene_folder /
                                  'stage_N_mean.txt',
                                  sep = '\t')
            # For "stage_M".
            stage_M.size().to_csv(gene_folder /
                                  'stage_M_count.txt', 
                                  sep = '\t', header = ['cases'])
            stage_M.median().to_csv(gene_folder /
                                    'stage_M_median.txt',
                                    sep = '\t')
            stage_M.mean().to_csv(gene_folder /
                                  'stage_M_mean.txt',
                                  sep = '\t')


            ### Remove "0" if any row contain "0". ###
            # Check if any value is "0" in first 3 columns.
            # Returns a Series of bool.
            df_gene_bool = (df_gene.iloc[:, 0:3] == 0).any()
            if df_gene_bool.any():
                print(f'Removing "0" for {gene}')
                # Remove "0".
                df_gene_rm0 = df_gene.copy()
                df_gene_rm0.iloc[:, 0:3] = df_gene_rm0.iloc[:, 0:3].\
                                           replace(0, np.nan)
                #print(df_gene_rm0.head())
                # Drop rows with missing value(nan).
                tumor_stage_rm0 = df_gene_rm0.dropna().groupby('tumor_stage')
                stage_T_rm0 = df_gene_rm0.dropna().groupby('stage_T')
                stage_N_rm0 = df_gene_rm0.dropna().groupby('stage_N')
                stage_M_rm0 = df_gene_rm0.dropna().groupby('stage_M')
                #print(tumor_stage_rm0.median())
                #print(tumor_stage_rm0.mean())


                ### Save to file. ###
                # For "tumor_stage".
                tumor_stage_rm0.size().to_csv(gene_folder /
                                              'stage_count_rm0.txt', 
                                              sep = '\t', header = ['cases'])
                tumor_stage_rm0.median().to_csv(gene_folder /
                                                'stage_median_rm0.txt',
                                                sep = '\t')
                tumor_stage_rm0.mean().to_csv(gene_folder /
                                              'stage_mean_rm0.txt',
                                              sep = '\t')
                # For "stage_T".
                stage_T_rm0.size().to_csv(gene_folder /
                                          'stage_T_count_rm0.txt', 
                                          sep = '\t', header = ['cases'])
                stage_T_rm0.median().to_csv(gene_folder /
                                            'stage_T_median_rm0.txt',
                                            sep = '\t')
                stage_T_rm0.mean().to_csv(gene_folder /
                                          'stage_T_mean_rm0.txt',
                                          sep = '\t')
                # For "stage_N".
                stage_N_rm0.size().to_csv(gene_folder /
                                          'stage_N_count_rm0.txt', 
                                          sep = '\t', header = ['cases'])
                stage_N_rm0.median().to_csv(gene_folder /
                                            'stage_N_median_rm0.txt',
                                            sep = '\t')
                stage_N_rm0.mean().to_csv(gene_folder /
                                          'stage_N_mean_rm0.txt',
                                          sep = '\t')
                # For "stage_M".
                stage_M_rm0.size().to_csv(gene_folder /
                                          'stage_M_count_rm0.txt', 
                                          sep = '\t', header = ['cases'])
                stage_M_rm0.median().to_csv(gene_folder /
                                            'stage_M_median_rm0.txt',
                                            sep = '\t')
                stage_M_rm0.mean().to_csv(gene_folder /
                                          'stage_M_mean_rm0.txt',
                                          sep = '\t')
            else:
                print(f'No "0" found for {gene}')

            # Generate histogram for each column.
            # Set plot style
            plt.style.use('ggplot')
            for column in columns_we_want[1:4]:
                print(f'Histogram for "df_gene" column: "{column}"')
                histogram(df_gene, column, gene_folder, '')
            for column in columns_we_want[1:4]:
                print(f'Histogram for "df_gene_rm0" column: "{column}"')
                histogram(df_gene_rm0, column, gene_folder, '_rm0')
    

            ### Merge the same stage. ###
            # Dictionary for stage rename.
            stage_merge = {'Stage_IA': 'Stage_I',
                           'Stage_IB': 'Stage_I',
                           'Stage_IIA': 'Stage_II',
                           'Stage_IIB': 'Stage_II',
                           'Stage_IIC': 'Stage_II',
                           'Stage_IIIA': 'Stage_III',
                           'Stage_IIIB': 'Stage_III',
                           'Stage_IIIC': 'Stage_III'}

            df_gene_merge = df_gene.replace(stage_merge).copy()

            # All stages: '[Not Available]', 'T0', 'T1~4', 'Tis', 'TX'.
            stage_merge_T = {'T1a': 'T1',
                             'T1b': 'T1',
                             'T2a': 'T2',
                             'T2b': 'T2',
                             'T3a': 'T3',
                             'T3b': 'T3',
                             'T4a': 'T4',
                             'T4b': 'T4'}

            # All stages: '[Not Available]', 'N0', 'N1~3', 'NX'.
            stage_merge_N = {'N1a': 'N1',
                             'N1b': 'N1',
                             'N2a': 'N2',
                             'N2b': 'N2',
                             'N2c': 'N2'}

            # All stages: '[Not Available]', 'M0', 'M1'.
            stage_merge_M = {'M1a': 'M1',
                             'M1b': 'M1',
                             'M1c': 'M1'}

            df_gene_merge['stage_T'] = df_gene['stage_T'].\
                                       replace(stage_merge_T).copy()
            df_gene_merge['stage_N'] = df_gene['stage_N'].\
                                       replace(stage_merge_N).copy()
            df_gene_merge['stage_M'] = df_gene['stage_M'].\
                                       replace(stage_merge_M).copy()
            df_gene_merge = df_gene_merge.set_index('Unnamed: 0')
            df_gene_merge.index.name = None    # Remove index name.
            #print(df_gene_merge.head(-5))


            ### Analysis. ###
            tumor_stage_merge = df_gene_merge.groupby('tumor_stage')
            stage_T_merge = df_gene_merge.groupby('stage_T')            
            stage_N_merge = df_gene_merge.groupby('stage_N')
            stage_M_merge = df_gene_merge.groupby('stage_M')
            # Count for each stage.
            #print(tumor_stage_merge.size())
            #print(tumor_stage_merge.median())
            #print(tumor_stage_merge.mean())
            # Dictionary for stage rename.


            ### Save to file. ###
            df_gene_merge.to_csv(f'{gene_folder}/merge_{gene}.txt', sep = '\t')
            tumor_stage_merge.size().to_csv(gene_folder /
                                            'merge_stage_count.txt',
                                            sep = '\t', header = ['cases'])
            tumor_stage_merge.median().to_csv(gene_folder /
                                              'merge_stage_median.txt',
                                              sep = '\t')
            tumor_stage_merge.mean().to_csv(gene_folder /
                                            'merge_stage_mean.txt',
                                            sep = '\t')
            # For "stage_T".
            stage_T_merge.size().to_csv(gene_folder /
                                        'merge_stage_T_count.txt', 
                                        sep = '\t', header = ['cases'])
            stage_T_merge.median().to_csv(gene_folder /
                                          'merge_stage_T_median.txt',
                                          sep = '\t')
            stage_T_merge.mean().to_csv(gene_folder /
                                        'merge_stage_T_mean.txt',
                                        sep = '\t')
            # For "stage_N".
            stage_N_merge.size().to_csv(gene_folder /
                                        'merge_stage_N_count.txt', 
                                        sep = '\t', header = ['cases'])
            stage_N_merge.median().to_csv(gene_folder /
                                          'merge_stage_N_median.txt',
                                          sep = '\t')
            stage_N_merge.mean().to_csv(gene_folder /
                                        'merge_stage_N_mean.txt',
                                        sep = '\t')
            # For "stage_M".
            stage_M_merge.size().to_csv(gene_folder /
                                        'merge_stage_M_count.txt', 
                                        sep = '\t', header = ['cases'])
            stage_M_merge.median().to_csv(gene_folder /
                                          'merge_stage_M_median.txt',
                                          sep = '\t')
            stage_M_merge.mean().to_csv(gene_folder /
                                        'merge_stage_M_mean.txt',
                                        sep = '\t')


            ### Remove "0" if any row contain "0". ###
            # Check if any value is "0" in first 3 columns.
            # Returns a Series of bool.
            df_gene_bool = (df_gene_merge.iloc[:, 0:3] == 0).any()
            if df_gene_bool.any():
                print(f'Removing "0" for {gene}')
                # Remove "0".
                df_gene_merge_rm0 = df_gene_merge.copy()
                df_gene_merge_rm0.iloc[:, 0:3] = df_gene_merge_rm0.iloc[:, 0:3].\
                                                 replace(0, np.nan)
                #print(df_gene_merge_rm0.head())
                # Drop rows with missing value(nan).
                tumor_stage_merge_rm0 = df_gene_merge_rm0.dropna().groupby('tumor_stage')
                stage_T_merge_rm0 = df_gene_merge_rm0.dropna().groupby('stage_T')
                stage_N_merge_rm0 = df_gene_merge_rm0.dropna().groupby('stage_N')
                stage_M_merge_rm0 = df_gene_merge_rm0.dropna().groupby('stage_M')
                #print(tumor_stage_merge_rm0.median())
                #print(tumor_stage_merge_rm0.mean())


                ### Save to file. ###
                # For "tumor_stage".
                tumor_stage_merge_rm0.size().to_csv(gene_folder /
                                                    'merge_stage_count_rm0.txt',
                                                    sep = '\t',
                                                    header = ['cases'])
                tumor_stage_merge_rm0.median().to_csv(gene_folder /
                                                      'merge_stage_median_rm0.txt',
                                                      sep = '\t')
                tumor_stage_merge_rm0.mean().to_csv(gene_folder /
                                                    'merge_stage_mean_rm0.txt',
                                                    sep = '\t')
                # For "stage_T".
                stage_T_merge_rm0.size().to_csv(gene_folder /
                                                'merge_stage_T_count_rm0.txt', 
                                                sep = '\t', header = ['cases'])
                stage_T_merge_rm0.median().to_csv(gene_folder /
                                                  'merge_stage_T_median_rm0.txt',
                                                  sep = '\t')
                stage_T_merge_rm0.mean().to_csv(gene_folder /
                                                'merge_stage_T_mean_rm0.txt',
                                                sep = '\t')
                # For "stage_N".
                stage_N_merge_rm0.size().to_csv(gene_folder /
                                                'merge_stage_N_count_rm0.txt', 
                                                sep = '\t', header = ['cases'])
                stage_N_merge_rm0.median().to_csv(gene_folder /
                                                  'merge_stage_N_median_rm0.txt',
                                                  sep = '\t')
                stage_N_merge_rm0.mean().to_csv(gene_folder /
                                                'merge_stage_N_mean_rm0.txt',
                                                sep = '\t')
                # For "stage_M".
                stage_M_merge_rm0.size().to_csv(gene_folder /
                                                'merge_stage_M_count_rm0.txt', 
                                                sep = '\t', header = ['cases'])
                stage_M_merge_rm0.median().to_csv(gene_folder /
                                                  'merge_stage_M_median_rm0.txt',
                                                  sep = '\t')
                stage_M_merge_rm0.mean().to_csv(gene_folder /
                                                'merge_stage_M_mean_rm0.txt',
                                                sep = '\t')

        
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


def histogram(df, x, save_path, name):
    fig, ax1 = plt.subplots()
    sns.histplot(df, x=x, ax=ax1)
    plt.tight_layout() ### Rescale the fig size to fit the data
    #plt.show()
    plt.savefig(f'{save_path}/hist_{x}{name}.png')


if __name__ == '__main__':
    main()