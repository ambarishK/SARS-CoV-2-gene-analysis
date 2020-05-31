from calculations.python.paths import *
from calculations.python.util import *
import Levenshtein

with open(data_path(REFERENCE_GENOME), 'r') as file:
    reference = ''.join(file.read().split('\n')[1:])

with open(data_path(REFERENCE_GENES), 'r') as file:
    genes = [g.split() for g in file.read().split('\n')]
    genes = [(g[0], translate(transcribe(reference[int(g[1]):int(g[2])]))) for g in genes]

with open(data_path("netmhcpan_bare_peptides.txt"), 'r') as input, open(data_path("netmhcpan_peptides.txt"), 'w') as output:
    for entry in input:
        entry = entry.strip()
        m = min(((gene_name, gene_data, begin) for gene_name, gene_data in genes for begin in range(len(gene_data) - len(entry))), key=lambda g: Levenshtein.distance(g[1][g[2]:g[2]+len(entry)], entry))
        print(f"{entry}\t{m[0]}\t{m[2]}", file=output)
