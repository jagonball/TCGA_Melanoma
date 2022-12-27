'''
To sort files downloaded from TCGA
customized for project TCGA-SKCM
at https://portal.gdc.cancer.gov/projects/TCGA-SKCM
'''
import os
import subprocess
import shutil
import json

def main():
    ### input parameters ###
    main_data = './data/nationwidechildrens.org_clinical_patient_skcm.txt'
    json_files = './data/files.2022-12-27.json'
    json_cases = './data/cases.2022-12-27.json'
    download_folder = './data/downloads/'
    download_compressed = 'gdc_download_20221227_162931.686561.tar.gz'
    

    # Unzip downloaded file
    subprocess.run(tar -xf download_compressed)
    

    # Read json
    with open(json_files) as file:
        files_contents = file.read()
    #print(files_contents)
    with open(json_cases) as file:
        cases_contents = file.read()
    

    # Load json contents into a list of dictionaries
    parsed_json_files = json.loads(files_contents)
    print(parsed_json_files[0]['case_id'])
    parsed_json_cases = json.loads(cases_contents)


if __name__ == '__main__':
    main()