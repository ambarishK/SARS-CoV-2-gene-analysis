from calculations.python.paths import *
from ast import literal_eval
from scipy import stats
import plotly.graph_objects as go
import plotly.express as px
import sys

GENES = ['E_gene', 'M_gene', 'S_gene', 'N_gene', 'orf1ab', 'ORF3a', 'ORF6', 'ORF7', 'ORF8', 'ORF10']
'''
//not tested
def calc_hist(filename: str) -> {str: [int]}:
    histogram_dict = {gene: {} for gene in GENES}
    with open(filename, 'r') as file:
        next(file)
        for line in file:
            entry = literal_eval(line)
            if len(entry[1]["genes"]) != len(GENES):
                continue
            gene_dict = {gene["name"]: gene for gene in entry[0]["genes"]}
            for gene in entry[1]["genes"]:
                if gene_dict[gene["name"]]["invalid_nucleotides"] != 0:
                    continue
                for mutation in gene["mutations"]:
                    if mutation["position"] in histogram_dict[gene["name"]]:
                        histogram_dict[gene["name"]][mutation["position"]] += 1
                    else:
                        histogram_dict[gene["name"]][mutation["position"]] = 1
    return {gene: [data[i] if i in data else 0 for i in range(max(data))] for gene, data in histogram_dict.items()}
'''

def del_zscore(list, zscore=20):
    z = abs(stats.zscore(list))
    listold = list
    list = []
    i = 0
    for item in z:
        if item > zscore:
            pass
        else:
            list.append(listold[i])
        i += 1
    return list

def icolorbar2(mut_dict):
    fig = go.Figure()
    i = 0
    for gene in mut_dict:
        data1 = [del_zscore(mut_dict[gene], 10)]
        gname = gene.replace('_', ' ')
        fig.add_trace(
            px.imshow(data1, labels=dict(x="Length of " + gname, color="No. of mutations")).data[0])
        i+=1

    def get_buttons():
        result=[]
        i=0
        for gene in mut_dict:
            ltestt = [False]*len(GENES)
            ltestt[i] = True
            a = dict(label=gene.replace('_', ' '), method="update", args=[{"visible": ltestt}, {"title": gene.replace('_', ' '), "annotations": []}])
            result.append(a)
            i+=1
        return result

    fig.update_layout(
        updatemenus=[
            dict(
                buttons=get_buttons(),
            )
        ])
    fig.update_traces(hovertemplate='Position: %{x} <br>No. of mutations: %{z} <extra></extra>')
    fig.update_yaxes(showticklabels=False)
    fig.write_html(sys.stdout, full_html=False, include_plotlyjs='cdn')
