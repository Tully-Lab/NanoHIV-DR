#!/usr/bin/env python3
'''
Short Python script to parse json returned from sierrapy query and store sample metadata in tab delimited text file.
Fields include sampleID, month, actualMonth and year plus a column for each drug and corresponding resistance mutation (if any).
Modified so it can handle UVRI style labels
'''

import json
import argparse
import datetime
import pandas as pd


parser = argparse.ArgumentParser()
parser.add_argument('--json', required=True, help='input json file containing query results')
parser.add_argument('--output', required=False, help='name of tab-delimited text file containing sample metadata')
args = parser.parse_args()

# get user input and set file names
json_in = ""
output_file=""

if args.json:
    json_in = args.json
else:
    print ("json file cannot be read by parser")

now = datetime.datetime.now()

if now.month < 10:
    time_now = (str(now.day) + "-0" +str(now.month) + "-" + str(now.year))
elif now.day <10:
    time_now = ("0" +str(now.day) + "-" +str(now.month) + "-" + str(now.year))
elif now.month < 10 and now.day <10: 
    time_now = ("0" +str(now.day) + "-0" +str(now.month) + "-" + str(now.year))
else:
    time_now = (str(now.day) + "-" +str(now.month) + "-" + str(now.year))

if args.output:
    output_file = args.output
else:
    output_file = time_now + "_DRM-overview.txt"


# parse sierrapy json output for relevant info

start_month = 11
start_year = 2004

append_list = []
with open(json_in) as json_file:
    data = json.load(json_file)
    for i in data:
        sample = i["inputSequence"]["header"]
        
        name_list = sample.split("_")
        if len(name_list)==4:
          day = name_list[0]
          month = name_list[1]
          year = name_list[2]
          ID = name_list[3]
          sample_month = ((int(year) - start_year) *12) + (int(month) - start_month) + 1
          append_dict = {'SampleID': ID, 'Year': year, 'Real_Month': month, 'Month': sample_month}
          header_list = ['SampleID', 'Year', 'Real_Month', 'Month', 'ATV', 'DRV', 'FPV', 'IDV' , 'LPV', 'NFV', 'SQV', 'TPV', 'ABC', 'AZT', 'D4T', 'DDI', 'FTC', 'LMV', 'TDF', 'DOR', 'EFV', 'ETR', 'NVP', 'RPV']

        else:
          append_dict = {'SampleID': sample}
          header_list = ['SampleID', 'ATV', 'DRV', 'FPV', 'IDV' , 'LPV', 'NFV', 'SQV', 'TPV', 'ABC', 'AZT', 'D4T', 'DDI', 'FTC', 'LMV', 'TDF', 'DOR', 'EFV', 'ETR', 'NVP', 'RPV']

        for j in i["drugResistance"]:
            for k in j["drugScores"]:
                drug_name = k["drug"]["name"]
                drug_score = k["score"]
                if drug_score != 0.0:
                    mutation_list=[]
                    for m in k["partialScores"]:
                        for key,value in m.items():
                            if key =='mutations':
                                for x, y in value[0].items():
                                    if x == 'text':
                                        mutation_list.append(y)
                    append_dict[drug_name] = mutation_list
                else:
                    append_dict[drug_name] = 0
        append_list.append(append_dict)
        
df = pd.DataFrame(append_list, columns = header_list)
df.set_index('SampleID', inplace=True)
df.to_csv(output_file, sep='\t')



