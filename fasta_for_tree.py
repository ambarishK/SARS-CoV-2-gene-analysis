import pandas as pd


df = pd.read_csv('Genome_data-X.csv')
dfi = df['Accession ID']

accession_list = dfi.to_list()
#for f in accession_list:
   #print(f)

accession_list_set = set(accession_list)

df = df.set_index('Accession ID')

output = []
with open('Cleaned_up_genes.fasta', 'r') as infile, open('filtered.fasta', 'w') as outfile:
    bad = True
    sequence = ''
    aid = ''
    temp_dict = {}
    bad = True
    accession_id = ''
    for line in infile:
        if line.startswith('>'):

            accession_id = line.split('|')[1]

            if accession_id in accession_list_set:
                bad = False
                key =line.replace(' ','')
                temp_dict[key] = ''

            elif not accession_id in accession_list_set:
                bad=True
        else:
            if not bad:
                aid = accession_id
                temp_dict[key] += line[int(df.loc[aid]['orf1ab_start']):int(df.loc[aid]['orf1ab_end'])]
                bad = True
    for key in temp_dict:
        output.append(key.rstrip())
        output.append(temp_dict[key].rstrip())
    #print(Dict)
print(len(temp_dict))


with open("for_tree.fasta", "w") as f:
    for line in output:
        f.write("%s\n" % line)

with open("for_tree.fasta", "r") as f:
    counter = 0
    for line in f:
        counter+=1



