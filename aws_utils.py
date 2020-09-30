#! /usr/bin/env python


from pathlib import Path
import os
import smart open

def get_bucket_path(bucket_name, object_key):

    """format path and credentials to read from S3"""
    home = str(Path.home())
    keys = {}
    for line in open(home + '/.aws/credentials').readlines():
        if line.startswith('['):
            pass
        else:
            key, val = line.rstrip().split(' = ')
            keys[key] = val
    path = 's3://{}:{}@{}/{}'.format(keys['aws_access_key_id'],keys['aws_secret_access_key'],
                                    bucket_name, object_key,)

    #path = get_bucket_path('kg-data-raw', 'pe_compounds.csv')
    #cp = pd.read_csv(smart_open.open(path))

    return (path)
