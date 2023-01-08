'''
To check the "duplicate" files.
Copy the "duplicate" files into folder "dupe_files_check".
'''

import os
import shutil
from pathlib import Path
from glob import glob
from file_sorter import create_folder, read_json, search_target_files


def main():
    ### Input parameters. ###
    project_name = 'TCGA_SKCM'
    output_folder = Path('C:/Repositories/Melanoma_TCGA/analysis/')
    download_folder = Path('C:/Repositories/Melanoma_TCGA/data/0107_all/')
    temp_folder = download_folder / 'temp_folder'
    json_files = Path('C:/Repositories/Melanoma_TCGA/data/files.2023-01-07.json')
    json_cases = Path('C:/Repositories/Melanoma_TCGA/data/cases.2023-01-07.json')
    # Files to rename at first '.'.
    files_re_list_1 = ['*.wxs.aliquot_ensemble_masked.maf',
                       '*.rna_seq.augmented_star_gene_counts.tsv',
                       '*.mirbase21.isoforms.quantification.txt',
                       '*.mirbase21.mirnas.quantification.txt']
    # Files to rename at second '.'.
    files_re_list_2 = ['*.gene_level_copy_number.v36.tsv']

    check_folder = create_folder('dupe_files_check', download_folder,
                                 verbose = True)

    # Read json files.
    parsed_json_files = read_json(json_files)
    #print(f'Number of files in files.json: {len(parsed_json_files)}')
    #print(parsed_json_files[0:2])
    parsed_json_cases = read_json(json_cases)
    #print(f'Number of cases cases.json: {len(parsed_json_cases)}')

    # Check temp_folder for remaining files.
    files_list = os.listdir(temp_folder)
    print(f'Number of files found: {len(files_list)}')
    target_files_list_1 = search_target_files(files_re_list_1, temp_folder)
    print(f'Numbers of files found to rename at first ".": '
          f'{len(target_files_list_1)}')
    #print(target_files_list_1)
    target_files_list_2 = search_target_files(files_re_list_2, temp_folder)
    print(f'Numbers of files found to rename at second ".": '
          f'{len(target_files_list_2)}')


    copy_check_files(target_files_list_1,
                     files_json = parsed_json_files,
                     cases_json = parsed_json_cases,
                     output_folder = output_folder,
                     project_name = project_name,
                     check_folder = check_folder,
                     delimiter = '.', name_pos = 1)
    copy_check_files(target_files_list_2,
                     files_json = parsed_json_files,
                     cases_json = parsed_json_cases,
                     output_folder = output_folder,
                     project_name = project_name,
                     check_folder = check_folder,
                     delimiter = '.', name_pos = 2)

    
def copy_check_files(files_list, files_json, cases_json,
                     output_folder, project_name, check_folder,
                     delimiter = '.', name_pos = 1):
    # Iterate files_list, find "duplicate" file in output_folder.
    for file in files_list:
        target_file_name = Path(file).name
        name_keep = target_file_name.split(delimiter)[name_pos:]
        name_keep = delimiter.join(name_keep)
        # Add ".gz" for the ".maf" files.
        if target_file_name.endswith('.maf'):
            target_file_name = target_file_name + '.gz'
        # Iterate files_json for matching file_name.
        target_case_id = next((dict['cases'][0]['case_id']\
                               for dict in files_json\
                               if dict['file_name'] == target_file_name),
                              None)
        target_submitter_id = next((dict['submitter_id']\
                                    for dict in cases_json\
                                    if dict['case_id'] == target_case_id),
                                   None)
        target_folder = str(output_folder / project_name / '*' \
                            / target_submitter_id / '*')
        # List of files in "target_submitter_id" folder.
        case_files = glob(target_folder)
        # Copy the "duplicate" file to "check_folder".
        file_to_check = [p for p in case_files if f'{target_submitter_id}.{name_keep}' in p]
        print(file_to_check)
        file_to_check_name = Path(file_to_check[0]).name
        #print(file_to_check)
        shutil.copyfile(file_to_check[0], check_folder / file_to_check_name)
        # Copy and rename the "file" to "check_folder"
        new_name = check_folder / f'{target_submitter_id}.2.{name_keep}'
        shutil.copyfile(file, new_name)


if __name__ == '__main__':
    main()