import argparse
from calculations.python.paths import *
from typing import Callable, Optional

MIN_GENOME_LEN = 29000

def perform_on_genomes(action: Callable[[str, str], None], input_filename: str) -> None:
    with open(input_filename, 'r') as file:
        header = None
        contents = ''
        for line in file:
            if line.startswith(">"):
                if header is not None:
                    action(header, contents)
                contents = ''
                header = line.rstrip()
                continue
            else:
                contents += line.rstrip()
        action(header, contents)

def filter_genome(genome: str) -> Optional[str]:
    if len(genome) < MIN_GENOME_LEN or any((g in genome) for g in ['-', 'B', 'H', 'D', 'V', 'W', 'S', 'K', 'R', 'Y']):
        return None
    return genome.upper()

def filter_fasta_file(input_filename: str, output_filename: str) -> None:
    valid_genomes = 0
    all_genomes = 0
    with open(output_filename, 'w') as f:
        def filter_and_save_genome(header: str, genome: str):
            nonlocal all_genomes, valid_genomes
            all_genomes += 1
            filtered = filter_genome(genome)
            if filtered is not None:
                valid_genomes += 1
                print(header.replace(' ', '').rstrip(), file=f)
                print(filtered, file=f)
        perform_on_genomes(filter_and_save_genome, input_filename)
    print('The total number of valid genomes is {}, {}% of all genomes passed'.format(valid_genomes, round(valid_genomes / all_genomes * 100)))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", type=str, help="input file", default=data_path(UNFILTERED_GENOMES))
    parser.add_argument("-o", type=str, help="output file", default=data_path(FILTERED_GENOMES))
    args = parser.parse_args()
    filter_fasta_file(args.i, args.o)
