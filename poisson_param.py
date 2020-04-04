import pandas as pd

df = pd.read_csv('Genome_data-X.csv')
#df = df.sort_values(by='Date', ascending=False)
a = str(df['Date'].head(1)).split('-')
print(a)
df['filter'] = df.apply(lambda row: len(str(row.Date).split('-')), axis=1)
df = df[df['filter'] == 3]
df['Year'] = df.apply(lambda row: str(row.Date).replace('\n','-').split('-')[0].split(' ')[-1], axis=1)
df['Month']  = df.apply(lambda row: int(str(row.Date).split('-')[1]), axis=1)
df['Day']  = df.apply(lambda row: int(str(row.Date).split('-')[2]), axis=1)
df['sum'] = df.apply(lambda row: row.E_gene + row.M_gene + row.N_gene + row.S_gene + row.orf1ab + row.ORF3a+ row.ORF6 + row.ORF7 + row.ORF8 +row.ORF10, axis = 1)
FirstVal = df[df['Month'] == 2]
FirstVal = FirstVal[FirstVal['Day'] <= 7 ]
print(len(FirstVal))


SecondVal = df[df['Month'] == 2]
SecondVal = SecondVal[SecondVal['Day'] <= 28]
SecondVal = SecondVal[SecondVal['Day'] >21]


Genes = ['E_gene', 'M_gene', 'S_gene', 'N_gene', 'orf1ab', 'ORF3a', 'ORF6', 'ORF7', 'ORF8', 'ORF10','whole_genome','sum']

vals = {}


for gene in Genes:
    print(gene)
    print('1st')
    print(FirstVal[gene].mean())
    print('2nd')
    print(SecondVal[gene].mean())
    print('____')
    vals[gene] = SecondVal[gene].mean()/(FirstVal[gene].mean() + 0.00000000001)

print('....')
print(vals['sum']) ## it should be around 1.35 as it is because 4 of February is around 45 day of Epidemy and 24 is around 66