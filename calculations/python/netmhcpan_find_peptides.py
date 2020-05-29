from calculations.python.paths import *
from calculations.python.util import *
import Levenshtein

with open(data_path(REFERENCE_GENOME), 'r') as file:
    reference = translate(transcribe(''.join(file.read().split('\n')[1:])))

with open(data_path("netmhcpan_bare_peptides.txt"), 'r') as input, open(data_path("netmhcpan_peptides.txt"), 'w') as output:
    for entry in input:
        m = min(range(len(reference) - len(entry)), key=lambda x: Levenshtein.distance(reference[x:x + len(entry)], entry))
        print(f"{entry.strip()}\t{m}", file=output)
