import pandas as pd
import ast
import csv
import itertools
from calculations.python.paths import *


def distances_tocsv(file=data_path(EXTRA_COMPARISONS_RESULTS), output_file=data_path(COMPARISONS_CSV),
                    compare_to_ref=False):

    with open(file, 'r') as file, open(output_file, 'w') as csvfile:
        next(file)
        counter = 0
        for line in file:
            counter+=1

            if type(eval(line)) == type(()):
                temp_dict = (eval(line))[0]

                if compare_to_ref:
                    temp_dict_ref = eval(line)[1]

            else:
                temp_dict = eval(line)

                compare_to_ref = False

            temp_dict2 = {gene['name'] + '_' + k: gene[k] for gene in temp_dict['genes' if 'genes' in temp_dict else 'tests'] for k in gene if k != 'name'}
            print(temp_dict2)

            del temp_dict['genes' if 'genes' in temp_dict else 'tests']
            for k, v in temp_dict.items():
                temp_dict2[k] = v
            writer = csv.DictWriter(csvfile, fieldnames=temp_dict2.keys())
            if counter == 1:
                writer.writeheader()
            writer.writerow(temp_dict2)
    df = pd.read_csv(output_file)
    for i in df:
        if 'mutations' in i:
            for x in df[i]:
                print(x)
            print(df)

# distances_tocsv()
#distances_tocsv(data_path(CALCULATED_GENOMES_DATA), data_path(DISTANCES_CSV), True)
distances_tocsv('../../data/tests_results.txt', '../../data/pcr.csv', True)
#df = pd.read_csv('../../data/distances.csv')


def merge_csv(file1, file2, output_file):
    df1 = pd.read_csv(file1)  # that one should be classic csv with N's
    df2 = pd.read_csv(file2)  # 2 jest mniejsza

    df2.rename(columns={'compared': 'header'}, inplace=True)
    df2['date'] = df2.apply(lambda row: row.header.split('|')[2], axis=1)

    common = df1.merge(df2, on='header', how='inner')
    common.to_csv(output_file)


# merge_csv(data_path(DISTANCES_CSV), data_path(COMPARISONS_CSV), data_path(MERGED_CSV))

if __name__ == "__main__":
    pass
    # merge_csv('genome_data-test.csv', 'Test-new-format.csv')
