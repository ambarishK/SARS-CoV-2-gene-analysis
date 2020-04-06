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
        dict['date'] = a
        for gene in List:
            dict[gene[0]] = gene[3]
            dict[gene[0] + '_start'] = gene[1]
            dict[gene[0] + '_end'] = gene[2]
            dict[gene[0] + '_len'] = gene[2]-gene[1]
        list.append(dict)
    temp_cols = []
    for i in ['E_gene', 'M_gene', 'S_gene', 'N_gene', 'orf1ab', 'ORF3a', 'ORF6', 'ORF7', 'ORF8', 'ORF10','S_gene_translation', 'ORF7_translation', 'ORF3a_translation', 'orf1ab_translation', 'ORF6_translation', 'M_gene_translation', 'ORF8_translation', 'N_gene_translation', 'ORF10_translation', 'E_gene_translation', 'whole_genome']:
        temp_cols.append(i+'_start')
        temp_cols.append(i+'_end')
        temp_cols.append(i+'_len')
    csv_columns = ['Accession ID', 'date', 'E_gene', 'M_gene', 'S_gene', 'N_gene', 'orf1ab', 'ORF3a', 'ORF6', 'ORF7', 'ORF8', 'ORF10','S_gene_translation', 'ORF7_translation', 'ORF3a_translation', 'orf1ab_translation', 'ORF6_translation', 'M_gene_translation', 'ORF8_translation', 'N_gene_translation', 'ORF10_translation', 'E_gene_translation', 'whole_genome']
    csv_columns+=temp_cols

    with open('Genome_Data-31_03_2020.csv', 'w') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=csv_columns)
        writer.writeheader()
        for data in list:
           writer.writerow(data)

df = pd.read_csv('Genome_Data-31_03_2020.csv')


Ref_len = {'E_gene': 227, 'M_gene': 668, 'S_gene': 3821, 'N_gene': 1259, 'orf1ab': 21289, 'ORF3a': 827, 'ORF6': 185, 'ORF7': 365, 'ORF8': 365, 'ORF10': 116}


for gene in ['E_gene', 'M_gene', 'S_gene', 'N_gene', 'orf1ab', 'ORF3a', 'ORF6', 'ORF7', 'ORF8', 'ORF10']:

    df = df[df[gene+'_len'] == Ref_len[gene]+1]

    x= gene+'_translation'
    df.pop(x+'_start')
    df.pop(x+'_end')
    df.pop(x+'_len')
print(df['whole_genome_end'].mean())

df.pop('whole_genome' + '_start')
df.pop('whole_genome' + '_end')
df.to_csv('Genome_data-X.csv', index=False)





##### PREVIOUS SCRIPT TO MAKE CSV BASED O DATASET FROM GISAID SITE
# df1 = pd.read_csv('csv_to_merge.csv').set_index('Accession ID') ### csv_to_merge.csv is from bigd.big.ac.cn ncov
# df2 = pd.read_csv('quality.csv').set_index('Accession ID')
# df3 = df1.join(df2, how = 'inner')
# df3 = df3.loc[df3['Sequence Quality'] != 'Low']
# df3 = df3.loc[df3['Nuc.Completeness'] == 'Complete']
# df3.to_csv('Genome_data.csv')
    # with open('test.csv', 'w') as f:
    #     for key in dict.keys():
    #         f.write("%s,%s\n" % (key, dict[key]))
