#!/usr/bin/env python3

import json
import requests

def parse_json(input_file):
    ''' Wrapper to colect the data from parse_node. Currently, results are organized as a list for each isolate. First field contains isolate name, second field mutation data and third field all the metadata'''
    
    data = []

    def parse_node(input_file, array):    
        '''
        Search recursively across a json dictionary and return the keys corresponding to name, branch_attrs (mutations info), node_attrs (metadata, e.g. authors).
        ''' 
        for record in input_file['children']:
            isolate = record['name']                   
            if 'branch_attrs' in record.keys():
                mutations = record['branch_attrs']['mutations'] 
            else:
                mutations = []
            if 'node_attrs' in record.keys():
                metadata = record['node_attrs'] 
            else:
                metadata = []
            array.append([isolate, mutations, metadata])
            if 'children' in record.keys():
                parse_node(record, array)
        return array

    results = parse_node(input_file, data)
    return results

def get_nextstrain_data():
    ''' Fetch data from nextstrain raw data site.'''

    url = 'https://data.nextstrain.org/ncov.json'
    r = requests.get(url)
    nextstrain_json = r.json()
    results = parse_json(nextstrain_json['tree'])
    return results

if __name__ == '__main__':
    results = get_nextstrain_data()
    print(results)
    
