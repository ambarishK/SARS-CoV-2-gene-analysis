import pandas as pd
import matplotlib.pyplot as plt
import ast
import csv
import itertools

with open('distances.txt', 'r') as file:
    list = []

    for header, genome in itertools.zip_longest(*[file] * 2):

        a = header.rstrip().split('|')[-1]
        List = ast.literal_eval(genome)
        dict = {'Accession ID' : header.rstrip().split('|')[1]}
        dict['Date'] = a
        for gene in List:
            dict[gene[0]] = gene[3]
        list.append(dict)
    print(list)
    csv_columns = ['Accession ID', 'Date', 'E_gene', 'M_gene', 'S_gene', 'N_gene', 'orf1ab', 'ORF3a', 'ORF6', 'ORF7', 'ORF8', 'ORF10','S_gene_translation', 'ORF7_translation', 'ORF3a_translation', 'orf1ab_translation', 'ORF6_translation', 'M_gene_translation', 'ORF8_translation', 'N_gene_translation', 'ORF10_translation', 'E_gene_translation', 'whole_genome']
    with open('csv_to_merge.csv', 'w') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=csv_columns)
        writer.writeheader()
        for data in list:
           writer.writerow(data)

##### PREVIOUS SCRIPT TO MAKE CSV BASED O DATASET FROM GISAID SITE
# df1 = pd.read_csv('csv_to_merge.csv').set_index('Accession ID')
# df2 = pd.read_csv('quality.csv').set_index('Accession ID')
# df3 = df1.join(df2, how = 'inner')
# df3 = df3.loc[df3['Sequence Quality'] != 'Low']
# df3 = df3.loc[df3['Nuc.Completeness'] == 'Complete']
# df3.to_csv('Genome_data.csv')
    # with open('test.csv', 'w') as f:
    #     for key in dict.keys():
    #         f.write("%s,%s\n" % (key, dict[key]))