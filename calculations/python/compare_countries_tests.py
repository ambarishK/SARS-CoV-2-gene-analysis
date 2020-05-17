import operator
import pandas as pd
from calculations.python.paths import *
from ast import literal_eval
import itertools
import matplotlib.pyplot as plt

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

def create_hist() -> ({str: [int]}, {str: [int]}):
    hist_dict_normalized = {}
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
            hist_dict_normalized[primer_name + "/" + country] = [values[i] / len(rows) if i in values else 0 for i in range(len(primers[primer_name]))]
            hist_dict[primer_name + "/" + country] = [values[i] if i in values else 0 for i in range(len(primers[primer_name]))]
    return hist_dict_normalized, hist_dict

def hist_to_suffixes(hist_dict: {str: [int]}, suffix_len: int=5) -> {str: [int]}:
    return {k: v[-suffix_len:] for k, v in hist_dict.items()}

hist_dict, for_df = create_hist()
hist_dict = hist_to_suffixes(hist_dict)
hist_df = pd.DataFrame(hist_to_suffixes(for_df))

def create_differences(hist_dict: {str: [int]}) -> {str: int}:
    differences = {}
    for test in primers:
        for country1, country2 in itertools.combinations(pcr_data['country'].unique(), 2):
            differences[f"{test}/{country1}.{test}/{country2}"] = sum((abs(a - b) for a, b in zip(hist_dict[f"{test}/{country1}"], hist_dict[f"{test}/{country2}"])))
    return differences

differences = create_differences(hist_dict)
top_differences = dict(sorted(differences.items(), key=operator.itemgetter(1), reverse=True)[:TOP_RESULTS])

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
