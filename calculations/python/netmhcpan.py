from calculations.python.paths import *
import tempfile, subprocess, sys

NETMHCPAN_LOCATION = data_path('netMHCpan')

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

with tempfile.NamedTemporaryFile('r') as temp_output:
    print(','.join(["header", "mhc"] + NetMHCpanResult.fields))
    for i in range(1, len(sys.argv)):
        mhc = sys.argv[i]
        print(f"Starting netMHCpan {mhc}", file=sys.stderr)
        subprocess.check_output(f'{NETMHCPAN_LOCATION} -a {mhc} -p -f {data_path("netmhcpan_input.txt")} > {temp_output.name}', shell=True)
        print(f"netMHCpan finished {mhc}", file=sys.stderr)
        temp_output.seek(0)
        with open(data_path("netmhcpan_input.txt"), 'r') as input:
            for result, header in zip(load_netmhcpan_results(temp_output), input):
                result_csv = ','.join(str(getattr(result, field)) for field in NetMHCpanResult.fields)
                print(header.rstrip(), end=',')
                print(mhc, end=',')
                print(result_csv)