import pandas as pd

#previous_genomes = pd.read_csv('sth.csv')
genomes_published = pd.read_csv('quality.csv')
output = []
with open('new_genes.fasta', 'r') as fasta_file:
    Temp_dict = {}
    for line in fasta_file:
        if line.startswith('>'):
            temp_header = line
            Temp_dict[temp_header] = ''

        if not line.startswith('>'):
            Temp_dict[temp_header] += line.upper()

    for header in Temp_dict:
        genome = Temp_dict[header].rstrip()
        head = header.rstrip().split('|')
        acession_id = header.rstrip().split('|')[1]
        temp_df = genomes_published.loc[genomes_published['Accession ID'] == acession_id]
        temp_df = temp_df.loc[temp_df['Nuc.Completeness'] == 'Complete']
        temp_df = temp_df.loc[temp_df['Sequence Quality'] != 'Low']
        if(len(temp_df) != 0 ):

            output.append(header)
            output.append(genome)

with open("Cleaned_up_genes.fasta", "w") as f:
    for line in output:
        f.write("%s\n" % line)