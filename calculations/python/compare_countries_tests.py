import operator
import pandas as pd
from calculations.python.paths import *
from ast import literal_eval
import itertools

TOP_RESULTS = 20

pcr_data = pd.read_csv(data_path('pcr.csv'))
pcr_data['country'] = pcr_data.apply(lambda row: row.header.split('/')[1], axis=1)

def get_primers() -> {str: str}:
    with open(data_path('primers_public.fas'), 'r') as ref_primers_file:
        names = list()
        primers = list()
        for line in ref_primers_file:
            if line.startswith('>'):
                names.append(line.rstrip()[1:])
            else:
                primers.append(line.rstrip())
    if len(names) != len(primers):
        raise ValueError("Invalid primers file.")
    return dict(zip(names, primers))

primers = get_primers()

def create_hist():
    hist_dict = {}
    for country, rows in pcr_data.groupby(['country']):
        for primer_name in primers:
            values = {}
            for mutations in rows[primer_name + "_mutations"].values:
                for m in literal_eval(mutations):
                    if m["position"] in values:
                        values[m["position"]] += 1
                    else:
                        values[m["position"]] = 0
            hist_dict[primer_name + "/" + country] = [values[i] if i in values else 0 for i in range(len(primers[primer_name]))]
    return hist_dict

hist_dict = create_hist()

def create_differences():
    differences = {}
    last_percent = 0
    counter = 0
    total = len(hist_dict) * (len(hist_dict) - 1) // 2
    for first, second in itertools.combinations(hist_dict, 2):
        counter += 1
        if counter * 100 // total != last_percent:
            last_percent = counter * 100 // total
            print(f"{last_percent}%", flush=True)
        hf = hist_dict[first]
        hs = hist_dict[second]
        differences[first + "." + second] = sum((abs(a - b) for a, b in zip(hf, hs))) + sum(hf[len(hs):]) + sum(hs[len(hf):])
    return differences

differences = create_differences()

top_differences = dict(sorted(differences.items(), key=operator.itemgetter(1), reverse=True)[:TOP_RESULTS])

print(top_differences)