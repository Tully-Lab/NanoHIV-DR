#!/usr/bin/env python3

'''
subtype_query.py: Short script to perform a custom sierrapy web query, parse
the returned json and extract subtype information + DRM list
'''

__author__ = "Samantha Campbell"
__email__= "samantha.campbell@glasgow.ac.uk"

import argparse
import subprocess
import os.path as osp
import datetime


parser = argparse.ArgumentParser()
parser.add_argument('--fasta', required=True, help='input fasta file')
parser.add_argument('--json', required=False, help='query results json filename')
parser.add_argument('--output', required=False, help='output filename')
args = parser.parse_args()

inputFasta = ""
output_json = ""

query_file = "custom_query.gql"


# Get user defined file names
if args.fasta:
    inputFasta = args.fasta
else:
    print ("Fasta file required!")

file_base = osp.basename(inputFasta)
file_id = osp.splitext(file_base)[0]

now = datetime.datetime.now()
time_now = (str(now.day) + "-" +str(now.month) + "-" + str(now.year))

base = time_now + "_" + file_id

if args.json:
    output_json = args.json
else:
    output_json = base + "_subtypes.json"


# check if query file exists - if not, generate it
if osp.exists(query_file):
    print ("Using custom query file:\t"+query_file+"\n")
else:
    print ("Generating custom query file...\n")
    customQuery= "inputSequence { header, \n}, \nsubtypeText, \nvalidationResults { \nlevel, \nmessage \n}, \
                \nalignedGeneSequences { \nfirstAA, lastAA, \ngene{ \nname, length }, \nDRM:mutations{text}\n},\ndrugResistance \
                {\ngene {\nname\n},\ndrugScores {\ndrugClass { name }, \ndrug {\nname,\ndisplayAbbr,\nfullName\n},\n \
                score,\npartialScores { \nmutations { \ntext,\n primaryType, \ncomments {type, text} \n}, \nscore \n}, \ntext \n}, \
                \nmutationsByTypes{ \nmutationType, \nmutations { \ntext} \n} \n}"

    
    with open (query_file, "w+") as out:
        out.write(customQuery)
    print ("Custom query file:\t" + query_file+"\n")    


# call sierrapy and run on input data
print("Performing sierrapy query...\n")

sierrapy_command = [
    "sierrapy",
    "fasta",
    inputFasta,
    "-q", query_file,
    "-o", output_json,
    "--no-sharding"
]

subprocess.run(sierrapy_command, check=True)

print(f"SierraPy result written to: {output_json}")

