'''
To sort files downloaded from TCGA
customized for project TCGA-SKCM
at https://portal.gdc.cancer.gov/projects/TCGA-SKCM
'''
import os
import subprocess
import pandas as pd
import shutil
import json

def main():
    ### input parameters ###
    main_data = './data/nationwidechildrens.org_clinical_patient_skcm.txt'
    json_files = './data/files.2022-12-27.json'
    json_cases = './data/cases.2022-12-27.json'
    download_folder = './data/downloads/'
    download_compressed = 'gdc_download_20221227_162931.686561.tar.gz'
    

    ### Unzip downloaded file ###
    # Check if 'extracted' folder exist, if not, create one.
    if not os.path.exists(f'{download_folder}extracted/'):
        print(f'Folder: {download_folder}extracted/ not found, creating...')
        os.mkdir(f'{download_folder}extracted/')
#    # Unzip to 'download_folder/extracted'
#    subprocess.run(f'tar -xf {download_folder}{download_compressed}\
#                         -C {download_folder}extracted/')
    # Read the 'MANIFEST.txt' file for file information
    file_info = pd.read_table(f'{download_folder}extracted/MANIFEST.txt',
                              low_memory=False)
    #print(file_info)
    # Filter for state 'validated', to ignore the 'annotations.txt'
    file_info = file_info[file_info['state'] == 'validated']
    #print(file_info)
    ## Move the files into 'temp_folder'
    # Check if 'temp_folder' folder exist, if not, create one.
    if not os.path.exists(f'{download_folder}temp_folder/'):
        print(f'Folder: {download_folder}temp_folder/ not found, creating...')
        os.mkdir(f'{download_folder}temp_folder/')
    for file in file_info['filename']:
        print(file)
        file_name = file.split('/')[1]
        print(file_name)
        shutil.move(f'{download_folder}extracted/{file}',
                    f'{download_folder}temp_folder/{file_name}')

'''
To do list:
1. Create folders based on tumor stage information
from file './data/nationwidechildrens.org_clinical_patient_skcm.txt'
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
2. Move and rename files in 
f'{download_folder}temp_folder/' to the above created folders
based on json_files 'case_id', 'file_name',
rename file with: json_cases 'submitter_id' and json_files 'data_category'
'''

#    # Read json
#    with open(json_files) as file:
#        files_contents = file.read()
#    #print(files_contents)
#    with open(json_cases) as file:
#        cases_contents = file.read()
#    
#
#    # Load json contents into a list of dictionaries
#    parsed_json_files = json.loads(files_contents)
#    print(parsed_json_files[0]['case_id'])
#    parsed_json_cases = json.loads(cases_contents)


if __name__ == '__main__':
    main()