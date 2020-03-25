#import pandas as pd
a=0 ## checking len of fasta
#previous_genomes = pd.read_csv('sth.csv')
#genomes_published = pd.read_csv('quality.csv')
output = []
with open('new_genes.fasta', 'r') as fasta_file:
    Temp_dict = {}
    for line in fasta_file:
        if line.startswith('>'):
            temp_header = line
            Temp_dict[temp_header] = ''

        if not line.startswith('>'):
            Temp_dict[temp_header] += line.upper()
    bad_seq = []
    for header in Temp_dict:
            #Quality Check
        if len(Temp_dict[header].split('N')) >3 or len(Temp_dict[header].split('-')) >2 or len(Temp_dict[header].rstrip())<29000:
            bad_seq.append(header)
            #Completeness Check
        elif len(Temp_dict[header].split('B'))>1 or len(Temp_dict[header].split('H'))>1 or len(Temp_dict[header].split('D'))>1 or len(Temp_dict[header].split('V'))>1 or len(Temp_dict[header].split('W'))>1 or len(Temp_dict[header].split('S'))>1 or len(Temp_dict[header].split('K'))>1 or len(Temp_dict[header].split('W'))>1 or len(Temp_dict[header].split('R'))>1 or len(Temp_dict[header].split('Y'))>1:
            bad_seq.append(header)
    for bad in bad_seq:
        Temp_dict.pop(bad)


    for header in Temp_dict:
        genome = Temp_dict[header].rstrip()
        head = header.rstrip().split('|')
        # acession_id = header.rstrip().split('|')[1]
        # temp_df = genomes_published.loc[genomes_published['Accession ID'] == acession_id]
        # temp_df = temp_df.loc[temp_df['Nuc.Completeness'] == 'Complete']
        # temp_df = temp_df.loc[temp_df['Sequence Quality'] != 'Low']
        # if(len(temp_df) != 0 ):
        a+=1
        output.append(header)
        output.append(genome)
print(a)
with open("Cleaned_up_genes.fasta", "w") as f:
    for line in output:
        f.write("%s\n" % line)