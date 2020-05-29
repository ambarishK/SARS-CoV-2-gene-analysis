from calculations.python.paths import *
from calculations.python.find_tests import load_genomes as load_genomes_raw
import ast, tempfile, subprocess, sys

MAX_ENTRIES = 4
NETMHCPAN_LOCATION = data_path('netMHCpan')
MHC_ALLELES = ['HLA-A02:01', 'HLA-A01:01']

def get_genome_part_to_analyze(genome, contents):
    for gene in genome[0]["genes"]:
        if gene["name"] == "E_gene":
            return contents[gene["begin"]:gene["begin"] + 8]
    raise ValueError

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
        try:
            part = get_genome_part_to_analyze(genome, contents)
            good += 1
        except ValueError:
            skipped += 1
            continue
        if len(part) < 8:
            raise ValueError("Lengths too short (<8)")
        if part not in parts:
            parts[part] = [(counter, header)]
        else:
            parts[part].append((counter, header))
    if len(parts) == 0:
        raise ValueError("Empty input or all skipped")
    for part in parts:
        print(part, file=temp_input)
    temp_input.flush()
    print(','.join(["header", "mhc"] + NetMHCpanResult.fields))
    for mhc in MHC_ALLELES:
        print(f"Starting netMHCpan {mhc}", file=sys.stderr)
        subprocess.check_output(f'{NETMHCPAN_LOCATION} -a {mhc} -p -f {temp_input.name} > {temp_output.name}', shell=True)
        print(f"netMHCpan finished {mhc}", file=sys.stderr)
        temp_output.seek(0)
        for result, parts_entry in zip(load_netmhcpan_results(temp_output), parts.items()):
            result_csv = ','.join(str(getattr(result, field)) for field in NetMHCpanResult.fields)
            for _, header in parts_entry[1]:
                print(header, end=',')
                print(mhc, end=',')
                print(result_csv)
    print(f"{skipped} genomes skipped, {good} good", file=sys.stderr)