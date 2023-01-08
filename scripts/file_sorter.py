'''
To sort files downloaded from TCGA, based on tumor stage and case barcode.
customized for project TCGA-SKCM
at https://portal.gdc.cancer.gov/projects/TCGA-SKCM

Download 7 types of files:
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
    '*.PDF'
# Data Type: Protein Expression Quantification
    '*_RPPA_data.tsv'
Rename with 'submitter_id' + originl name after second '.':
# Data Type: Gene Level Copy Number
    '*.gene_level_copy_number.v36.tsv'

To do list:
1. search within star_gene_counts.tsv for LINC00094, MIR1270,
calculate mean value for each stage folder.
'''
import os
import sys
import subprocess
import re
from pathlib import Path
from glob import glob
import pandas as pd
import numpy as np
import shutil
import json

def main():
    ### Input parameters. ###
    project_name = 'TCGA_SKCM'
    output_folder = Path('C:/Repositories/Melanoma_TCGA/analysis/')
    json_files = Path('C:/Repositories/Melanoma_TCGA/data/files.2023-01-07.json')
    json_cases = Path('C:/Repositories/Melanoma_TCGA/data/cases.2023-01-07.json')
    # Download data is ".tar.gz" or with manifest and folders.
    download_is_tar_gz = False
    ## !!! Warning: download_folder must only contain
    ## either downloaded ".tar.gz" file
    ## or manifest file plus folders downloaded with the manifest.
    ## Other files will be deleted by this code.
    download_folder = Path('C:/Repositories/Melanoma_TCGA/data/0107_all/')
    download_compressed = 'gdc_download_20230106_023217.521073.tar.gz'
    manifest_file = 'gdc_manifest_20230107_155741.txt'
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
    # Main data with all cases' information.
    main_data = Path('C:/Repositories/Melanoma_TCGA/data/clinical_patient_skcm.txt')
    stage_colname = 'ajcc_pathologic_tumor_stage'
    barcode_colname = 'bcr_patient_barcode'


    ### Create project folder if not already exist. ###
    project_name = replace_special_chars(project_name)
    project_folder = create_folder(project_name, output_folder,
                                   verbose = True)


    ### Managing downloaded file. ###
    # Create "extracted" folder if not already exist.
    extracted_folder = create_folder('extracted', download_folder,
                                      verbose = True)
    # Unzip downloaded ".tar.gz"
    if download_is_tar_gz:
        file_check(download_compressed, 'compressed file', download_folder)
        # Unzip compressed file to extracted_folder.
        print(f'## Extracting "{download_compressed}" into "{extracted_folder}"...')
        compressed_file = download_folder / download_compressed
        subprocess.run(['tar', '-xf', compressed_file, '-C', extracted_folder])
        # Set the "manifest_file" name.
        manifest_file = 'MANIFEST.txt'
    # Move the manifest file and folders to "extracted" folder.
    else:
        download_list = os.listdir(download_folder)
        file_check(manifest_file, 'manifest file', download_folder)
        # Remove "extracted" from download_list if exists.
        if 'extracted' in download_list:
            download_list.remove('extracted')
        #print(len(download_list))
        print(f'## Moving the downloaded files to "{extracted_folder}"...')
        move_files_in_list(download_list, download_folder, extracted_folder)


    # Read the "manifest_file" for file information.
    file_info = pd.read_table(extracted_folder / manifest_file,
                              low_memory=False)
    #print(file_info)
    # Filter for state 'validated', to ignore the 'annotations.txt'.
    file_info = file_info[file_info['state'] == 'validated']
    #print(file_info)
    # Create 'temp_folder' folder if not already exist.
    temp_folder = create_folder('temp_folder', download_folder,
                                 verbose=True)
    print(f'## Moving the extracted files to "{temp_folder}"...')
    move_files_in_list(file_info['filename'], extracted_folder, temp_folder)
    # Remove extracted_folder.
    print(f'## Move completed, deleting folder "{extracted_folder}"...')
    shutil.rmtree(extracted_folder) #, ignore_errors=True)


    ### Create a list to store target files' path. ###
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


    ### Rename files with case id. ###
    # Read json files.
    parsed_json_files = read_json(json_files)
    #print(parsed_json_files[0:2])
    parsed_json_cases = read_json(json_cases)
    # Rename files.
    rename_target_files(target_files_list_1, temp_folder, json_files,
                        files_json = parsed_json_files,
                        cases_json = parsed_json_cases,
                        delimiter = '.', id_pos = 0, name_pos = 1)
    rename_target_files(target_files_list_2, temp_folder, json_files,
                        files_json = parsed_json_files,
                        cases_json = parsed_json_cases,
                        delimiter = '.', id_pos = 1, name_pos = 2)


    ### Create folders based on tumor stage and case barcode, ###
    ### then move the renamed files from "temp_folder" into ###
    ### the corresponding folder. ###
    # Read main_data into DataFrame.
    main_df = pd.read_table(main_data, sep = '\t',
                            header = 0, skiprows = [1,2])
    #print(main_df.head())
    # Select tumor stage column, then reorder the unique stage values.
    unique_stages = np.sort(main_df[stage_colname].unique())
    print(f'Unique stages: {unique_stages}')
    ## Create a folder for each tumor stage,
    ## then for each patient create folder with bcr_patient_barcode.
    for stage in unique_stages:
        print(f'Working on: {stage}')
        # Check if the stage folder exist, if not, create one.
        stage_mod = replace_special_chars(stage)
        #print(stage_mod)
        stage_folder = create_folder(stage_mod, project_folder)
        # Subset main DataFrame with specific stage.
        stage_df = main_df[main_df[stage_colname] == stage]
        # Iterate "bcr_patient_barcode" to create and move files.
        for barcode in stage_df[barcode_colname]:
            #print(barcode)
            # Create case folder with barcode within "stage_folder".
            case_folder = create_folder(barcode, stage_folder)
            ## Move files from temp_folder to case_folder.
            case_list = [f'{barcode}*']
            case_files = search_target_files(case_list, temp_folder)
            #print(f'Found "{len(case_files)}" files for case "{barcode}"')
            move_files_in_list(case_files, temp_folder, case_folder)
            
    ## Remove temp_folder.
    #print(f'## Move completed, deleting folder {temp_folder}...')
    #shutil.rmtree(temp_folder) #, ignore_errors=True)


def create_folder(folder_name, path_to_folder, verbose = False):
    """Create folder if not already exist.

    :param folder_name: The folder name.
    :type folder_name: str
    :param path_to_folder: The path to check for the folder.
    :type path_to_folder: str
    :return: The full path to inside the folder.
    :rtype: str
    """
    final_folder = path_to_folder / folder_name
    if verbose:
        print(f'Checking if folder "{folder_name}" '
              f'exists in "{path_to_folder}"...')
    if folder_name not in os.listdir(path_to_folder):
        if verbose:
            print(f'## Folder "{folder_name}" not found, creating...')
        os.mkdir(final_folder)    
    return final_folder


def replace_special_chars(str):
    """Replace special characters with "_".

    :param str: Input string.
    :type str: str
    :return: Output string.
    :rtype: str
    """
    mod_str = re.sub('[^a-zA-Z0-9 \n\.]', '_', str)
    mod_str = mod_str.replace(" ", "_")
    return mod_str


def move_files_in_list(files_list, from_folder, to_folder):
    """Move files in the list or Series, from "from_folder" to "to_folder".
    
    :param files_list: A list or Series of file names.
    :type files_list: list or Series
    :param from_folder: The original folder path.
    :type from_folder: Path
    :param to_folder: The destination folder path.
    :type to_folder: Path
    """
    for file in files_list:
        #print(file)
        file_name = Path(file).name
        #print(f'file_name: {file_name}')
        shutil.move(from_folder / file,
                    to_folder / file_name)


def file_check(file, file_type, folder):
    """Check if file is in folder.

    :param file: The file name.
    :type file: str
    :param file_type: The type of file.
    :type file_type: str
    :param folder: The folder to look for.
    :type folder: Path or str
    """
    if file not in os.listdir(folder):
        print(f'Error: Cannot find the {file_type} '
              f'"{file}", please check again.')
        sys.exit()


def read_json(json_file):
    """Read json into a list of dictionaries.

    :param json: Path to json file.
    :type json: Path or str
    :return: Parsed json contents.
    :rtype: list
    """
    with open(json_file) as file:
        contents = file.read()
    parsed_json = json.loads(contents)
    return parsed_json


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
        target_files = str(folder_path / name)
        target_files_list += glob(target_files)
    return target_files_list


def rename_target_files(files_list, folder_path, json_files,
                        files_json, cases_json,
                        delimiter = '.', id_pos = 0, name_pos = 1):
    """Rename files in "files_list" within "folder_path", 
    based on information of two json files.

    :param files_list: A list of file names, accept regular expression.
    :type files_list: list
    :param folder_path: Folder path of the files.
    :type folder_path: Path or str
    :param json_files: File path of json files.
    :type json_files: Path or str
    :param files_json: Parsed files json.
    :type files_json: list
    :param cases_json: Parsed cases json.
    :type cases_json: list
    :param delimiter: Delimiter of file name, defaults to '.'
    :type delimiter: str, optional
    :param id_pos: Position of "ID" of the specified delimiter, defaults to 0
    :type id_pos: int, optional
    :param name_pos: Position of "Name to keep" of the specified delimiter, defaults to 1
    :type name_pos: int, optional
    """
    for i in files_list:
        target_file_name = Path(i).name
        # Files to keep the name after first ".".
        file_id = target_file_name.split(delimiter)[id_pos]
        #print(file_id)
        name_keep = target_file_name.split(delimiter)[name_pos:]
        name_keep = delimiter.join(name_keep)
        #print(name_keep)
        # Iterate files_json for matching file_name.
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
            old_name = folder_path / target_file_name
            new_name = folder_path/f'{target_submitter_id}.{name_keep}'
            if not new_name.exists():
                os.rename(old_name, new_name)


if __name__ == '__main__':
    main()