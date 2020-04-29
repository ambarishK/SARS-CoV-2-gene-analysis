# import importlib
# importlib.reload(module)
import filter_fasta
import findgenemutations
import tocsv

def process_fasta(file):
    filter_fasta.filter_fasta_file(file)
    findgenemutations.count_all()
    tocsv.distances_tocsv()

    print('Fasta file parsed, results are in Genome_Data-X.csv')

if __name__ == "__main__":
    process_fasta('gisaid_cov2020_sequences.fasta')