'''
To sort files downloaded from TCGA
customized for project TCGA-SKCM
at https://portal.gdc.cancer.gov/projects/TCGA-SKCM
'''
import os
import subprocess
import re
from pathlib import Path
import glob
import pandas as pd
import numpy as np
import shutil
import json

def main():
    ### input parameters ###
    project_name = 'TCGA_SKCM'
    output_folder = Path('C:/Repositories/Melanoma_TCGA/analysis/')
    main_data = Path('C:/Repositories/Melanoma_TCGA/data/clinical_patient_skcm.txt')
    json_files = Path('C:/Repositories/Melanoma_TCGA/data/files.2023-01-07.json')
    json_cases = Path('C:/Repositories/Melanoma_TCGA/data/cases.2023-01-07.json')
    download_folder = Path('C:/Repositories/Melanoma_TCGA/data/0106_test/')
    download_compressed = 'gdc_download_20230106_023217.521073.tar.gz'


    ### Create project folder if not already exist
    project_name = replace_special_chars(project_name)
    project_folder = create_folder(project_name, output_folder)


    ### Unzip downloaded file ###
    # Create 'extracted' folder if not already exist.
    extracted_folder = create_folder('extracted', download_folder)
    #if not os.path.exists(f'{download_folder}extracted/'):
    #    print(f'Folder: {download_folder}extracted/ not found, creating...')
    #    os.mkdir(f'{download_folder}extracted/')
    # Unzip compressed file to extracted_folder
    print(f'## Extracting the compressed file to {extracted_folder}...')
    compressed_file = download_folder / download_compressed
    subprocess.run(['tar', '-xf', compressed_file, '-C', extracted_folder])
    # Read the 'MANIFEST.txt' file for file information
    file_info = pd.read_table(extracted_folder/'MANIFEST.txt',
                              low_memory=False)
    #print(file_info)
    # Filter for state 'validated', to ignore the 'annotations.txt'
    file_info = file_info[file_info['state'] == 'validated']
    #print(file_info)
    ## Move the files into 'temp_folder'
    # Create 'temp_folder' folder if not already exist.
    temp_folder = create_folder('temp_folder', download_folder)
    #if not os.path.exists(f'{download_folder}temp_folder/'):
    #    print(f'Folder: {download_folder}temp_folder/ not found, creating...')
    #    os.mkdir(f'{download_folder}temp_folder/')
    print(f'## Moving the extracted files to {temp_folder}...')
    for file in file_info['filename']:
        #print(file)
        file_name = file.split('/')[-1]
        #print(file_name)
        shutil.move(extracted_folder/file,
                    temp_folder/file_name)


    ### Read main_data into DataFrame
    main_df = pd.read_table(main_data, sep = '\t', header = [0, 1, 2])
    #print(main_df['ajcc_pathologic_tumor_stage'])
    ## Select tumor stage column, squeeze into Series,
    ## then reorder the unique stage values.
    unique_stages = np.sort(main_df['ajcc_pathologic_tumor_stage']\
                            .squeeze().unique())
    print(f'Unique stages: {unique_stages}')
    ## Create a folder for each tumor stage,
    ## then for each patient create folder with bcr_patient_barcode
    for stage in unique_stages:
        #print(stage)
        # Check if the stage folder exist, if not, create one.
        stage = replace_special_chars(stage)
        #print(stage)
        stage_folder = create_folder(stage, project_folder)
        # Subset main DataFrame with specific stage
        #stage_df = main_df[main_df['ajcc_pathologic_tumor_stage'] == stage]

    # Read json.
    with open(json_files) as file:
        files_contents = file.read()
    #print(files_contents)
    with open(json_cases) as file:
        cases_contents = file.read()

    # Load json contents into a list of dictionaries.
    parsed_json_files = json.loads(files_contents)
    #print(parsed_json_files[0:2])#['case_id'])
    parsed_json_cases = json.loads(cases_contents)
    #print(parsed_json_cases)

    # Files that do not require rename.
    files_re_list_0 = ['*.PDF',
                       '*_RPPA_data.tsv']
    # Files to rename at first '.'.
    files_re_list_1 = ['*.wxs.aliquot_ensemble_masked.maf.gz',
                       '*.rna_seq.augmented_star_gene_counts.tsv',
                       '*.mirbase21.isoforms.quantification.txt',
                       '*.mirbase21.mirnas.quantification.txt']
    # Files to rename at second '.'.
    files_re_list_2 = ['*.gene_level_copy_number.v36.tsv']


    # Create a list to store target files' path.
    target_files_list_0 = search_target_files(files_re_list_0, temp_folder)
    print(f'Numbers of files found that do not require rename: '
          f'{len(target_files_list_0)}')
    target_files_list_1 = search_target_files(files_re_list_1, temp_folder)
    print(f'Numbers of files found to rename at first ".": '
          f'{len(target_files_list_1)}')
    #print(target_files_list_1)
    target_files_list_2 = search_target_files(files_re_list_2, temp_folder)
    print(f'Numbers of files found to rename at second ".": '
          f'{len(target_files_list_2)}')
    rename_target_files(target_files_list_1, temp_folder, json_files,
                        files_json = parsed_json_files,
                        cases_json = parsed_json_cases,
                        delimiter = '.', id_pos = 0, name_pos = 1)
    rename_target_files(target_files_list_2, temp_folder, json_files,
                        files_json = parsed_json_files,
                        cases_json = parsed_json_cases,
                        delimiter = '.', id_pos = 1, name_pos = 2)

    
    




def create_folder(folder_name, output_folder):
    """Create folder if not already exist.

    :param folder_name: The folder name.
    :type folder_name: str
    :param output_folder: The path to check for the folder.
    :type output_folder: str
    :return: The full path to inside the folder.
    :rtype: str
    """
    final_folder = output_folder / folder_name
    print(f'Checking if folder "{folder_name}" exists...')
    if folder_name not in os.listdir(output_folder):
        print(f'## Folder "{folder_name}" not found, creating...')
        os.mkdir(final_folder)    
    return final_folder


def replace_special_chars(str):
    mod_str = re.sub('[^a-zA-Z0-9 \n\.]', '_', str)
    mod_str = mod_str.replace(" ", "_")
    return mod_str


def search_target_files(file_list, folder_path):
        """Search files matching "file_list" in "folder_path".

        :param file_list: a list of files. (accept regular expression)
        :type file_list: list
        :param folder_path: The folder path to search for.
        :type folder_path: Path or str
        :return: A list of target files.
        :rtype: list
        """
        target_files_list = []
        for name in file_list:
            target_files = str(folder_path/name)
            target_files_list += glob.glob(target_files)
        return target_files_list


def rename_target_files(files_list, folder_path, json_files,
                        files_json, cases_json,
                        delimiter = '.', id_pos = 0, name_pos = 1):
        for i in files_list:
            target_file_name = Path(i).name
            # Files to keep the name after first ".".
            file_id = target_file_name.split(delimiter)[id_pos]
            #print(file_id)
            name_keep = target_file_name.split(delimiter)[name_pos:]
            name_keep = delimiter.join(name_keep)
            #print(name_keep)
            # Iterate through files_json for matching file_name.
            target_case_id = None
            target_submitter_id = None
            for i in files_json:
                if i['file_name'] == target_file_name:
                    #print(i)
                    if len(i['cases']) != 1:
                        print(f'Warning: list "cases" for file '
                              f'"{target_file_name}" in "{json_files}" '
                              f'are empty or more than one items')
                    target_case_id = i['cases'][0]['case_id']
                    #print(target_case_id)
                    for i in cases_json:
                        if i['case_id'] == target_case_id:
                            #print(i)
                            target_submitter_id = i['submitter_id']
                            #print(target_submitter_id)
            # Rename target file if there's match.
            if target_submitter_id:
                os.rename(folder_path/target_file_name,
                          folder_path/f'{target_submitter_id}.{name_keep}')


'''
To do list:
1. Create folders based on tumor stage information
from file './data/clinical_patient_skcm.txt'
    [Not Available]
    I/II NOS
    Stage 0
    Stage I
    Stage IA
    Stage IB
    Stage II
    Stage IIA
    Stage IIB
    Stage IIC
    Stage III
    Stage IIIA
    Stage IIIB
    Stage IIIC
    Stage IV
2. 7 files for each case:
First 4 rename with 'submitter_id' + originl name after first '.':
    # Data Type: Masked Somatic Mutation
    '*.wxs.aliquot_ensemble_masked.maf.gz' (unzip then delete .gz)
    # Data Type: Gene Expression Quantification
    '*.rna_seq.augmented_star_gene_counts.tsv'
    # Data Type: Isoform Expression Quantification
    '*.mirbase21.isoforms.quantification.txt'
    # Data Type: miRNA Expression Quantification
    '*.mirbase21.mirnas.quantification.txt'
Keep these 2:
    # Data Type: Pathology Report
    '*.80F131FA-7E24-4210-800F-3E6F442CAB6F.PDF'
    # Data Type: Protein Expression Quantification
    '*TCGA-D3-A3ML-06A-21-A241-20_RPPA_data.tsv'
Rename with 'submitter_id' + originl name after second '.':
    # Data Type: Gene Level Copy Number
    '*.gene_level_copy_number.v36.tsv'
3. Move and rename files in 
f'{download_folder}temp_folder/' to the above created folders
based on json_files 'case_id', 'file_name',
rename file with: json_cases 'submitter_id' and json_files 'data_category'
4. search within star_gene_counts.tsv for LINC00094, MIR1270,
calculate mean value for each stage folder.
'''




if __name__ == '__main__':
    main()