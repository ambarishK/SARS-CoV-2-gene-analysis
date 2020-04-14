def filter_fasta_file(file):
    output = []
    a=0
    with open(file, 'r') as fasta_file:
        Temp_dict = {}
        for line in fasta_file:
            if line.startswith('>'):
                temp_header = line
                Temp_dict[temp_header] = ''
            else:
                Temp_dict[temp_header] += line.upper().rstrip()
        bad_seq = []
        b = len(Temp_dict) #number of all genomes passed
        for header in Temp_dict:
                #Quality Check
            if  len(Temp_dict[header].rstrip())<29000 or len(Temp_dict[header].split('-')) >1 or len(Temp_dict[header].split('N')) >15:
                bad_seq.append(header)
                #Completeness Check
            elif len(Temp_dict[header].split('B'))>1 or len(Temp_dict[header].split('H'))>1 or len(Temp_dict[header].split('D'))>1 or len(Temp_dict[header].split('V'))>1 or len(Temp_dict[header].split('W'))>1 or len(Temp_dict[header].split('S'))>1 or len(Temp_dict[header].split('K'))>1 or len(Temp_dict[header].split('W'))>1 or len(Temp_dict[header].split('R'))>1 or len(Temp_dict[header].split('Y'))>1:
                bad_seq.append(header)
        for bad in bad_seq:
            Temp_dict.pop(bad)

        for header in Temp_dict:
            genome = Temp_dict[header].rstrip()
            # head = header.rstrip().split('|')
            # acession_id = header.rstrip().split('|')[1]
            # temp_df = genomes_published.loc[genomes_published['Accession ID'] == acession_id]
            # temp_df = temp_df.loc[temp_df['Nuc.Completeness'] == 'Complete']
            # temp_df = temp_df.loc[temp_df['Sequence Quality'] != 'Low']
            # if(len(temp_df) != 0 ):
            a+=1
            output.append(header.replace(' ','').rstrip())
            output.append(genome.rstrip())
    print('The total number of valid genomes is {}, {}% of all genomes passed'.format(a, round(a/b*100)))
    with open("Cleaned_up_genes.fasta", "w") as f:
        for line in output:
            f.write("%s\n" % line)

if __name__ == "__main__":
    filter_fasta_file('gisaid_cov2020_sequences.fasta')
