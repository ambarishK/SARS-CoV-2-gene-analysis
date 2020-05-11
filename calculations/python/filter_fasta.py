import argparse


def filter_fasta_file(input_filename: str, output_filename: str) -> None:
    MIN_GENOME_LEN = 29000

    genomes = {}
    with open(input_filename, 'r') as fasta_file:
        for line in fasta_file:
            line = line.rstrip()
            if len(line) == 0:
                continue
            if line.startswith('>'):
                header = line
                genomes[header] = ''
            else:
                genomes[header] += line.upper()
    bad_seq = []
    all_genomes = len(genomes)
    for header, genome in genomes.items():
        if len(genome) < MIN_GENOME_LEN or any((g in genome) for g in ['-', 'B', 'H', 'D', 'V', 'W', 'S', 'K', 'R', 'Y']):
            bad_seq.append(header)
    for bad in bad_seq:
        del genomes[bad]
    valid_genomes = 0
    with open(output_filename, "w") as f:
        for header, genome in genomes.items():
            valid_genomes += 1
            print(header.replace(' ', '').rstrip(), file=f)
            print(genome, file=f)
    print('The total number of valid genomes is {}, {}% of all genomes passed'.format(valid_genomes, round(valid_genomes / all_genomes * 100)))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", type=str, help="input file", default="gisaid_hcov-19_2020_04_29_22.fasta")
    parser.add_argument("-o", type=str, help="output file", default="Cleaned_up_genes.fasta")
    args = parser.parse_args()
    filter_fasta_file(args.i, args.o)
