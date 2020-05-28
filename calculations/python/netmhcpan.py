from calculations.python.paths import *
from calculations.python.find_tests import load_genomes as load_genomes_raw
import ast, tempfile, subprocess, sys

MAX_ENTRIES = 100
NETMHCPAN_LOCATION = data_path('netMHCpan')
MHC_ALLELE = 'HLA-A02:01'
PEPTIDE = False

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
    next(file)
    for _ in file:
        next(file) #skip dashes
        next(file) #skip labels
        line = next(file) #skip dashes
        data = line.split()
        data = [e for e in data if len(e) > 0]
        if len(data) != 13:
            print(line)
            raise ValueError("Wrong format?")
        results.append(NetMHCpanResult(int(data[0]), data[3], int(data[4]), int(data[5]), int(data[6]), int(data[7]), int(data[8]), data[9], float(data[11]), float(data[12])))
        line = next(file)
        if not line.startswith('-'):
            raise ValueError("Many entries?")
        next(file)
        next(file)
    return results

'''
---------------------------------------------------------------------------------------------------------------------------
 Pos         MHC        Peptide      Core Of Gp Gl Ip Il        Icore        Identity  Score_EL %Rank_EL BindLevel
---------------------------------------------------------------------------------------------------------------------------
   1 HLA-A*02:01     SSQCVNLTTR SQCVNLTTR  0  1  1  0  0   SSQCVNLTTR           asdf0 0.0001840   34.071
   1 HLA-A*02:01      SSQCVNLTT SSQCVNLTT  0  0  0  0  0    SSQCVNLTT           asdf0 0.0001460   37.000
   2 HLA-A*02:01      SQCVNLTTR SQCVNLTTR  0  0  0  0  0    SQCVNLTTR           asdf0 0.0016100   14.832
   1 HLA-A*02:01       SSQCVNLT SSQ-CVNLT  0  0  0  3  1     SSQCVNLT           asdf0 0.0000100   75.000
   2 HLA-A*02:01       SQCVNLTT SQC-VNLTT  0  0  0  3  1     SQCVNLTT           asdf0 0.0000640   48.250
   3 HLA-A*02:01       QCVNLTTR QCVN-LTTR  0  0  0  4  1     QCVNLTTR           asdf0 0.0000220   63.889
---------------------------------------------------------------------------------------------------------------------------

Protein asdf0. Allele HLA-A*02:01. Number of high binders 0. Number of weak binders 0. Number of peptides 6

-----------------------------------------------------------------------------------
'''

with open(data_path(CALCULATED_GENOMES_DATA), 'r') as genomes, tempfile.NamedTemporaryFile('w') as temp_input, tempfile.NamedTemporaryFile('r') as temp_output:
    next(genomes)
    substrings_length = None
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
        if substrings_length is not None:
            if len(part) != substrings_length:
                raise ValueError("Uneven lengths")
        else:
            substrings_length = len(part)
        if part not in parts:
            parts[part] = [(counter, header)]
        else:
            parts[part].append((counter, header))
    if substrings_length is None:
        raise ValueError("Empty input or all skipped")
    if substrings_length < 8:
        raise ValueError("Lengths too short (<8)")
    for part in parts:
        print(">a", file=temp_input)
        print(part, file=temp_input)
    temp_input.flush()
    print(f"Starting netMHCpan", file=sys.stderr)
    subprocess.check_output(f'{NETMHCPAN_LOCATION} -a {MHC_ALLELE} {"-p " if PEPTIDE else ""}-l {substrings_length} -f {temp_input.name} > {temp_output.name}', shell=True)
    print(f"netMHCpan finished", file=sys.stderr)
    temp_output.seek(0)
    print(','.join(["header", "mhc"] + NetMHCpanResult.fields))
    for result, parts_entry in zip(load_netmhcpan_results(temp_output), parts.items()):
        result_csv = ','.join(str(getattr(result, field)) for field in NetMHCpanResult.fields)
        for _, header in parts_entry[1]:
            print(header, end=',')
            print(MHC_ALLELE, end=',')
            print(result_csv)
    print(f"{skipped} genomes skipped, {good} good", file=sys.stderr)