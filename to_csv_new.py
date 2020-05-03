import pandas as pd
import ast
import csv
import itertools


def distances_tocsv(file='extra_comparisons_results.txt', output_file = 'Test-new-format.csv'):

    with open(file, 'r') as file:
        data_list = []

        next(file)

        counter =0
        for line in file:
            if type(eval(line)) == int:
                pass


            counter+=1

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
            if counter > 10:
                break
        csv_columns =data_list[0].keys()
        with open(output_file, 'w') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=csv_columns)
            writer.writeheader()
            for data in data_list:
               writer.writerow(data)
    df = pd.read_csv(output_file)
distances_tocsv()
distances_tocsv('distances3.txt','genome_data-test.csv')

def merge_csv(file1,file2):
    df1= pd.read_csv(file1) # that one should be classic csv with N's
    df2 = pd.read_csv(file2)
    print(df1)
    df1.set_index('header')
    df2.set_index('compared')
    common = pd.merge(df1, df2, how='left')

    # for gene in ['E_gene', 'M_gene', 'S_gene', 'N_gene', 'orf1ab', 'ORF3a', 'ORF6', 'ORF7', 'ORF8', 'ORF10']:
    #     df = df[df[gene+'_len'] == Ref_len[gene]+1]
    #     x= gene+'_translation'
    #     y= gene+'_N'
    #     #print(df.sort_values(by=gene, ascending=False))
    #     df.pop(x + '_start')
    #     df.pop(x + '_end')
    #     df.pop(x + '_len')
    #     df.pop(y + '_start')
    #     df.pop(y + '_end')
    #     df.pop(y + '_len')
    # for x in ['E_gene_mutations', 'M_gene_mutations', 'S_gene_mutations', 'N_gene_mutations', 'orf1ab_mutations', 'ORF3a_mutations', 'ORF6_mutations', 'ORF7_mutations', 'ORF8_mutations', 'ORF10_mutations', 'E_gene_translation_changes', 'M_gene_translation_changes', 'S_gene_translation_changes', 'N_gene_translation_changes', 'orf1ab_translation_changes', 'ORF3a_translation_changes', 'ORF6_translation_changes', 'ORF7_translation_changes', 'ORF8_translation_changes', 'ORF10_translation_changes']:
    #     df.pop(x + '_start')
    #     df.pop(x + '_end')
    #     df.pop(x + '_len')
    # df['genomic_len'] = df['ORF10_end'] - df['orf1ab_start']
    # #print(df.sort_values(by='S_gene', ascending=False))
    # #df= df[df['genomic_len'] == 29409]
    # df.pop('whole_genome' + '_start')
    # df.pop('whole_genome' + '_end')
    # df.set_index('name')
    # df['sum'] = df.apply(lambda
    #                          row: row.E_gene + row.M_gene + row.N_gene + row.S_gene + row.orf1ab + row.ORF3a + row.ORF6 + row.ORF7 + row.ORF8 + row.ORF10,
    #                      axis=1)
    # df['sum_N'] = df.apply(lambda
    #                            row: row.E_gene_N + row.M_gene_N + row.N_gene_N + row.S_gene_N + row.orf1ab_N + row.ORF3a_N + row.ORF6_N + row.ORF7_N + row.ORF8_N + row.ORF10_N,
    #                        axis=1)
    #
    # #df['whole_genome'] = df.apply(lambda row: row.E_gene + row.M_gene + row.N_gene + row.S_gene + row.orf1ab + row.ORF3a+ row.ORF6 + row.ORF7 + row.ORF8 +row.ORF10, axis = 1)
    # df.to_csv(output_file, index=False)


if __name__ == "__main__":
    merge_csv('genome_data-test.csv','Test-new-format.csv')
