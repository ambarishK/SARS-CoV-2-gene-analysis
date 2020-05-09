import pandas as pd
import ast
import csv
import itertools


def distances_tocsv(file='extra_comparisons_results.txt', output_file='comparisons.csv'):
    with open(file, 'r') as file:
        data_list = []

        next(file)

        counter = 0
        for line in file:
            if type(eval(line)) == int:
                pass

            counter += 1

            if type(eval(line)) == type(()):
                temp_dict = (eval(line))[0]

            else:
                temp_dict = eval(line)
            print(temp_dict)

            temp_dict2 = {gene['name'] + '_' + k: gene[k] for gene in temp_dict['genes'] for k in gene if k != 'name'}
            del temp_dict['genes']
            for k, v in temp_dict.items():
                temp_dict2[k] = v

            data_list.append(temp_dict2)
        csv_columns = data_list[0].keys()
        with open(output_file, 'w') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=csv_columns)
            writer.writeheader()
            for data in data_list:
                writer.writerow(data)
    # df = pd.read_csv(output_file)


distances_tocsv()
#distances_tocsv('distances3.txt', 'distances.csv')


def merge_csv(file1, file2, output_file):
    df1 = pd.read_csv(file1)  # that one should be classic csv with N's
    df2 = pd.read_csv(file2)  # 2 jest mniejsza

    df2.rename(columns={'compared': 'header'}, inplace=True)
    df2['date'] = df2.apply(lambda row: row.header.split('|')[2], axis = 1 )

    common = df1.merge(df2, on='header', how='inner')
    common.to_csv(output_file)



merge_csv('distances.csv', 'comparisons.csv', 'genome_neigh.csv')



if __name__ == "__main__":
    pass
    # merge_csv('genome_data-test.csv', 'Test-new-format.csv')

