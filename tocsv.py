import pandas as pd
import ast
import csv
import itertools

def distances_tocsv(file='distances.txt'):
    output_file = 'Genome_Data-X.csv'
    with open(file, 'r') as file:
        list = []
        for header, genome in itertools.zip_longest(*[file] * 2):
            a = header.rstrip().split('|')[-1]
            List = ast.literal_eval(genome)
            dict = {'Accession ID' : header.rstrip().split('|')[1]}
            dict['date'] = a
            dict['name'] = header.replace(' ', '').rstrip().split('|')[0][:-1].replace('>', '')
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
        csv_columns = ['name','date','Accession ID', 'E_gene', 'M_gene', 'S_gene', 'N_gene', 'orf1ab', 'ORF3a', 'ORF6', 'ORF7', 'ORF8', 'ORF10','S_gene_translation', 'ORF7_translation', 'ORF3a_translation', 'orf1ab_translation', 'ORF6_translation', 'M_gene_translation', 'ORF8_translation', 'N_gene_translation', 'ORF10_translation', 'E_gene_translation', 'whole_genome']
        csv_columns+=temp_cols
        with open(output_file, 'w') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=csv_columns)
            writer.writeheader()
            for data in list:
               writer.writerow(data)
    df = pd.read_csv(output_file)
    Ref_len = {'E_gene': 227, 'M_gene': 668, 'S_gene': 3821, 'N_gene': 1259, 'orf1ab': 21289, 'ORF3a': 827, 'ORF6': 185, 'ORF7': 365, 'ORF8': 365, 'ORF10': 116}

    for gene in ['E_gene', 'M_gene', 'S_gene', 'N_gene', 'orf1ab', 'ORF3a', 'ORF6', 'ORF7', 'ORF8', 'ORF10']:
        df = df[df[gene+'_len'] == Ref_len[gene]+1]
        x= gene+'_translation'

        #print(df.sort_values(by=gene, ascending=False))
        df.pop(x + '_start')
        df.pop(x + '_end')
        df.pop(x + '_len')

    df['genomic_len'] = df['ORF10_end'] - df['orf1ab_start']
    #print(df.sort_values(by='S_gene', ascending=False))
    df= df[df['genomic_len'] == 29409]
    #print(df['genomic_len'])
    df.pop('whole_genome' + '_start')
    df.pop('whole_genome' + '_end')
    df.set_index('name')
    #df['whole_genome'] = df.apply(lambda row: row.E_gene + row.M_gene + row.N_gene + row.S_gene + row.orf1ab + row.ORF3a+ row.ORF6 + row.ORF7 + row.ORF8 +row.ORF10, axis = 1)
    df.to_csv(output_file, index=False)

if __name__ == "__main__":
    distances_tocsv()