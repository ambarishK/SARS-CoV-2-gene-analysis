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
    hist_dict2 = {}
    for country, rows in pcr_data.groupby(['country']):
        for primer_name in primers:
            values = {}
            for mutations in rows[primer_name + "_mutations"].values:
                for m in literal_eval(mutations):
                    if m["position"] in values:
                        values[m["position"]] += 1
                    else:
                        values[m["position"]] = 0
            hist_dict[primer_name + "/" + country] = [values[i] / len(rows) if i in values else 0 for i in range(len(primers[primer_name])-5,len(primers[primer_name]))]
            hist_dict2[primer_name + "/" + country] = [values[i] if i in values else 0 for i in range(len(primers[primer_name])-5,len(primers[primer_name]))]
    return hist_dict, hist_dict2

hist_dict, for_df = create_hist()
hist_df = pd.DataFrame(for_df)

import matplotlib.pyplot as plt


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
#
# print(hist_df['RdRp_HongKong_R/England'])#['RdRP_HongKong_R/England'])
# temp_data = hist_df['RdRp_HongKong_R/England'].to_list()
# plt.hist(temp_data)
# plt.show()
plt.style.use('ggplot')
def hist_comparison(data, country1, country2):
    plt.subplot(1, 2, 1)
    temp_data = data[country1].to_list()
    axis = [i for i in range(len(hist_df[country1].to_list()))]
    plt.bar(axis, temp_data)
    plt.title(country1)
    plt.subplot(1, 2, 2)
    temp_data = data[country2].to_list()
    axis = [i for i in range(len(hist_df[country1].to_list()))]
    plt.bar(axis,temp_data)
    plt.title(country2)
    plt.tight_layout()
    plt.show()
# a=0
# for x in top_differences:
#     if 'RdRp_HongKong_F' in x.split('.')[0] and 'RdRp_HongKong_F' in x.split('.')[1]:
#         hist_comparison(hist_df,x.split('.')[0], x.split('.')[1])
#         a+=1
#         if a==10:
#             break