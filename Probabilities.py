

import pandas as pd

df = pd.read_csv('Genome_Data.csv')

Genes = ['E_gene', 'M_gene', 'S_gene', 'N_gene', 'orf1ab', 'ORF3a', 'ORF6', 'ORF7', 'ORF8', 'ORF10']
Genes_len= {'E_gene' : 26245 - 26472,
'M_gene'  : 26523 - 27191,
'S_gene'  :21563 - 25384,
'N_gene' :  28274 - 29533,
'orf1ab' : 266 - 21555,
'ORF3a' : 25393 - 26220,
'ORF6' : 27202 - 27387,
'ORF7' : 27394 - 27759,
'ORF8' : 27894 - 28259,
'ORF10' : 29558 - 29674}


Genome = 29902
Mutations = df['whole_genome'].mean()

Data_values = {}
for gene in Genes:
    Data_values[gene] = df[gene].mean()

Expected_values = {}

for gene in Genes:
    Expected_values[gene] = -Genes_len[gene]/Genome*Mutations


Comparsion = {}

for gene in Genes:
    Comparsion[gene] =Data_values[gene]- Expected_values[gene]
print(Comparsion)