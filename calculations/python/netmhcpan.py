from calculations.python.paths import *
from calculations.python.util import *
from calculations.python.find_tests import load_genomes as load_genomes_raw
import ast, tempfile, subprocess, sys, Levenshtein

MAX_ENTRIES = 40
NETMHCPAN_LOCATION = data_path('netMHCpan')

if len(sys.argv) < 2:
    print('USAGE:\nnetmhcpan HLAs...', file=sys.stderr)
    exit(1)

def get_mhc():
    return (sys.argv[i] for i in range(1, len(sys.argv)))

class GenomeError(Exception):
    pass

def find_peptide_in_genome(gene_protein, peptide_data, peptide_location, search_radius=10, min_peptide_len=8):
    m = min((gene_protein[peptide_location + x:peptide_location + x + len(peptide_data)] for x in range(-search_radius, search_radius + 1)), key=lambda x: Levenshtein.distance(x, peptide_data))
    if len(m) >= min_peptide_len:
        return m
    raise GenomeError

with open(data_path('netmhcpan_peptides.txt'), 'r') as file:
    peptides = [(p[0], p[1], int(p[2])) for p in (p.split('\t') for p in file.read().split('\n') if len(p) > 0)]

def get_genome_part_to_analyze(genome, contents):
    result = {}
    for gene in genome[0]["genes"]:
        gene_protein = translate(transcribe(contents[gene["begin"]:gene["end"]]))
        for peptide_data, peptide_gene, peptide_location in peptides:
            if peptide_gene == gene["name"] and gene["invalid_nucleotides"] == 0:
                try:
                    result[peptide_data] = find_peptide_in_genome(gene_protein, peptide_data, peptide_location)
                except GenomeError:
                    pass
    return result

def load_genomes():
    for v, e in enumerate(load_genomes_raw()):
        if v > MAX_ENTRIES:
            break
        yield e

class NetMHCpanResult:
    fields = ["pos", "core", "of", "gp", "gl","ip", "il", "icore", "score_el", "bindlevel"]

    def __init__(self, pos, core, of, gp, gl, ip, il, icore, score_el, bindlevel):
        self.pos = pos
        self.core = core
        self.of = of
        self.gp = gp
        self.gl = gl
        self.ip = ip
        self.il = il
        self.icore = icore
        self.score_el = score_el
        self.bindlevel = bindlevel

    def __repr__(self):
        return str(self)

    def __str__(self):
        return str(self.to_dict())

    def to_dict(self):
        return {name: getattr(self, name) for name in NetMHCpanResult.fields}

def load_netmhcpan_results(file) -> [NetMHCpanResult]:
    results = []
    f = (line for line in file if len(line.rstrip().lstrip()) > 0 and not line.startswith('#'))
    file = iter(f)
    next(file) #skip Distance to training data
    next(file) #skip dashes
    next(file) #skip labels
    next(file) #skip dashes
    for line in file:
        if line.startswith('-'):
            break
        data = line.split()
        data = [e for e in data if len(e) > 0]
        if len(data) != 13:
            print(line)
            raise ValueError("Wrong format?")
        results.append(NetMHCpanResult(int(data[0]), data[3], int(data[4]), int(data[5]), int(data[6]), int(data[7]), int(data[8]), data[9], float(data[11]), float(data[12])))
    return results

with open(data_path(CALCULATED_GENOMES_DATA), 'r') as genomes, tempfile.NamedTemporaryFile('w') as temp_input, tempfile.NamedTemporaryFile('r') as temp_output:
    next(genomes)
    good = 0
    skipped = 0
    parts = {}
    for counter, e in enumerate(zip(genomes, load_genomes())):
        line, raw_genome = e
        header, contents = raw_genome
        genome = ast.literal_eval(line)
        part = get_genome_part_to_analyze(genome, contents)
        good += len(part)
        skipped += len(peptides) - len(part)
        for reference_peptide, peptide_in_genome in part.items():
            if peptide_in_genome not in parts:
                parts[peptide_in_genome] = [(counter, header, reference_peptide)]
            else:
                parts[peptide_in_genome].append((counter, header, reference_peptide))
    if len(parts) == 0:
        raise ValueError("Empty input or all skipped")
    for part in parts:
        print(part, file=temp_input)
    temp_input.flush()
    print(','.join(["header", "mhc", "reference_peptide", "matched_peptide"] + NetMHCpanResult.fields))
    for mhc in get_mhc():
        print(f"Starting netMHCpan {mhc}", file=sys.stderr)
        subprocess.check_output(f'{NETMHCPAN_LOCATION} -a {mhc} -p -f {temp_input.name} > {temp_output.name}', shell=True)
        print(f"netMHCpan finished {mhc}", file=sys.stderr)
        temp_output.seek(0)
        for result, parts_entry in zip(load_netmhcpan_results(temp_output), parts.items()):
            result_csv = ','.join(str(getattr(result, field)) for field in NetMHCpanResult.fields)
            for _, header, reference_peptide in parts_entry[1]:
                print(header, end=',')
                print(mhc, end=',')
                print(reference_peptide, end=',')
                print(parts_entry[0], end=',')
                print(result_csv)
    print(f"{skipped} entries skipped, {good} good", file=sys.stderr)
