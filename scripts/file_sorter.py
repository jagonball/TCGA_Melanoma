'''
To sort files downloaded from TCGA
'''
import json

def main():
    # Read json
    json_file = './data/cases.2022-12-27.json'
    with open(json_file) as file:
        file_contents = file.read()
    #print(file_contents)

    # Load into a list of dictionaries
    parsed_json = json.loads(file_contents)
    print(parsed_json[0]['case_id'])

if __name__ == '__main__':
    main()