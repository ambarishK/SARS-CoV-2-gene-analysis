import ast
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import statistics
from calculations.python.paths import *

# pd.set_option("display.max_rows", None, "display.max_columns", None)

# df = pd.read_csv('../../data/genome_data.csv')
# print(df)
# exit()
# print(os.getcwd())
# file= data_path(CALCULATED_GENOMES_DATA)
# file = '../../data/genome_data.csv'
DF = pd.read_csv(data_path('distances.csv'))
# DF = pd.read_csv(file, escapechar='\\')
DF['date'] = pd.to_datetime(DF['date'])  # converting data column to data type

DF['country'] = DF['header'].apply(lambda x: str(x).split('/')[-3])
# '''
# for x in DF:
#     if x != 'country' and x != 'header':
#         DF.drop(x, axis=1, inplace=True)
#
# df = DF.groupby(by='country').count()
# df.sort_values(by='header', ascending=False, inplace=True)
# print(df.head(30))
# '''

# print(len(DF))
# DF.sort_values(by='date', ascending=False, inplace=True)
# DF = DF.head((len(DF) - 6))  # deleting outlier data
# print(len(DF))
# string = 'S_gene_translation_changes'


# DF['orf1ab_mutations'] = DF.apply(lambda
#                                row: (ast.literal_eval(row.orf1ab_mutations)),
#                            axis=1)
# hist_data = []
# for cols in DF:
#     if '_mutations' in  cols or '_changes' in cols:
#         for sth in DF[cols]:
#             sth = ast.literal_eval(sth)

def change_format(list1):  # takes list and returns list counting how many there are 0s, 1s, 2s,... n's
    max_val = max(list1)
    value_counter = []
    for i in range(0, max_val + 1):
        value_counter.append(list1.count(i))
    return value_counter


# print(change_format([1,5,9]))

def save_dict(dict, filename):
    with open(filename, 'w') as f:
        for item in dict:
            f.write(str(item) + '\n')
            f.write(str(dict[item]) + '\n')
    print('Dict saved as ' + filename)


def read_dict(file, mode='normal'):
    dic = {}
    with open(file, 'r') as f:
        while True:
            line1 = f.readline().rstrip()
            line2 = f.readline().rstrip()
            if mode == 'list':  # if values in dict are lists
                if not line2:
                    break
                else:
                    line2 = ast.literal_eval(line2)
            if not line1: break
            if not line2: break
            dic[line1] = line2
    return dic


# plt.style.use('ggplot')
# fig, ax = plt.subplots(nrows=20,ncols=1, sharex=True, sharey=True)
genes = [column.split('_protein_mutations')[0] for column in DF.columns if column.endswith('_protein_mutations')]

def get_mutations(gene: str):
    return (mutation
            for mutations, invalid_nucleotides in zip(DF[gene + '_mutations'], DF[gene + '_invalid_nucleotides'])
            if invalid_nucleotides == 0
            for mutation in ast.literal_eval(mutations))

def make_histogram_dict(l):
    values = {}
    for el in l:
        if el in values:
            values[el] += 1
        else:
            values[el] = 1
    return values

def histogram_dict_to_list(d: {int: int}, length=None):
    return [d[i] if i in d else 0 for i in range(length)]

def calc_hist():
    hist_dict = {gene: histogram_dict_to_list(make_histogram_dict((g['position'] for g in get_mutations(gene))), max(abs(DF[gene + "_end"] - DF[gene + "_begin"])) + 1) for gene in genes}
    return hist_dict

# save_dict(hist_dict, 'dict_mut.txt')

mut = []
trans = []
# for x in genes:
#     mut.append(x + '_mutations')
#     trans.append(x + '_translation_changes')


# mut = ['E_gene_mutations', 'M_gene_mutations', 'S_gene_mutations', 'N_gene_mutations', 'orf1ab_mutations', 'ORF3a_mutations',
#        'ORF6_mutations', 'ORF7_mutations', 'ORF8_mutations', 'ORF10_mutations']

# trans = ['E_gene_translation_changes', 'M_gene_translation_changes', 'S_gene_translation_changes', 'N_gene_translation_changes',
#          'orf1ab_translation_changes', 'ORF3a_translation_changes', 'ORF6_translation_changes', 'ORF7_translation_changes', 'ORF8_translation_changes',
#          'ORF10_translation_changes']


def filter_country(country='Poland', data=DF):
    data['country'] = data['name'].apply(lambda x: str(x).split('/')[-3])
    dataframe2 = data[data['country'] == country]
    print('Number of genomes in {} is {}, {}% of whole'.format(country, len(dataframe2),
                                                               round(len(dataframe2) / len(data), 3)))
    return dataframe2


def get_lists_colorbar():
    lists = []
    for genename in genes:
        hist_data = []
        m=0
        DF1 = DF[DF[genename + '_invalid_nucleotides'] == 0]
        for location in DF1[genename + '_mutations']:
            m+=1
            try:
                location = ast.literal_eval(location)
            except: #when there is naan or other value- for debugging straight away
                print(type(location))
                print(location)
                print(genename)
                print(m)
                exit()
            if (len(location) != 0):
                for list in location:
                    print(list)
                    exit()
                    hist_data.append(list[0])
        new_list = change_format(hist_data)
        lists.append(new_list)
    return lists


# def llwrite(list, filename='mutation_density.pg'):
#     with open(filename, 'w') as f:
#         for _list in list:
#             f.write('#' + '\n')
#             for inty in _list:
#                 f.write(str(inty) + '\n')
#     print('File saved as ' + filename)


# llwrite(get_lists_colorbar(), 'mutation_density22test.pg')


# def llread(filename='mutation_density4.pg'):
#     with open(filename, 'r') as f:
#         result = []
#         i = 0
#         for line in f:
#             if line == '#\n':
#                 if i == 0:
#                     i += 1
#                 else:
#                     result.append(a)
#                 a = []
#             else:
#                 a.append(int(line[:-1]))
#         result.append(a)
#         return result


def del_zscore(list, zscore=20):
    z = abs(stats.zscore(list))
    # print('Previous max zscore was {}'.format(max(z)))
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


def colorbar():
    i = 0
    data = llread()
    fig, axs = plt.subplots(10, sharey=True, sharex=False)
    for list in data:
        list = del_zscore(list, 10)
        print(mut[i])
        label = mut[i][:-10]
        arr = np.array(list)
        arr = np.stack((arr, arr), axis=0)  # thicc1
        axs[i].imshow(arr, aspect='auto', cmap='plasma')
        axs[i].set_ylabel(label, rotation=0, labelpad=25, y=0.3)
        axs[i].axes.get_xaxis().set_ticks([])
        axs[i].axes.get_yaxis().set_ticks([])
        i += 1
    plt.suptitle('Mutation density in genes of SARS-CoV-2 virus', y=0.95)
    plt.savefig('./plots/density/mutation_density2.png', dpi=600)
    plt.show()


def icolorbar(gene):
    import plotly.express as px
    i = 0
    for gene1 in genes:
        if gene1 == gene:
            break
        i += 1
    data = llread()
    data = [data[i]]
    gname = gene.replace('_', ' ')
    fig = px.imshow(data,
                    labels=dict(x="Length of " + gname, color="No. of mutations"),
                    y=[gname]
                    )
    # fig.update_xaxes(side="top")
    fig.update_traces(hovertemplate='Place in ' + gname + ': %{x} <br>No. of mutations: %{z} <extra></extra>')  #
    # fig.update_traces(hovertemplate=None, selector={'name': 'Europe'})
    name = 'div_density_chart_' + gene + '.html'
    fig.write_html('./plots/html/' + name, full_html=False, include_plotlyjs='cdn')
    print('File saved in directory plots/html/ as ' + name)
    fig.show()


# icolorbar('S_gene')
# for i in genes:
#     icolorbar(i)


def icolorbar2(mut_dict):
    import plotly.graph_objects as go
    import plotly.express as px
    # data = llread(file)
    fig = go.Figure()
    i=0
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
            ltest = [False]*len(genes)
            ltestt = ltest
            ltestt[i] = True
            a = dict(label=gene.replace('_', ' '), method="update", args=[{"visible": ltestt}, {"title": gene.replace('_', ' '), "annotations": []}])
            result.append(a)
            i+=1

        return result

    fig.update_layout(
        updatemenus=[
            dict(
                active=4,
                buttons=get_buttons(),
            )
        ])

    fig.update_traces(hovertemplate='Place in gene: %{x} <br>No. of mutations: %{z} <extra></extra>')
    name = 'div_density_chart_allgenes.html'
    fig.update_yaxes(showticklabels=False)
    fig.write_html('' + name, full_html=False, include_plotlyjs='cdn')
    print('File saved in directory plots/html/ as ' + name)

    fig.show()


# icolorbar2(read_dict('dict_mut.txt', 'list'))


# a = read_dict('dict_mut.txt', 'list')
#
# for x in a:
#     print(x)
#     print(max(a[x]))

def colorbar1(gene):
    i = 0
    for gene1 in mut:
        if gene1.startswith(gene):
            break
        i += 1
    data = llread()
    fig, ax = plt.subplots()
    arr = np.array(del_zscore(data[i], 20))
    print('Mean number of mutations is {}'.format(round(statistics.mean(map(float, arr))), 2))
    arr = np.stack((arr, arr), axis=0)  # thicc1
    im = ax.imshow(arr, extent=(0, 1, 0, 0.1), cmap='plasma')
    ax.axes.get_yaxis().set_ticks([])
    ax.set_xlabel('Length of ' + gene + ' (' + str(len(data[i])) + ')', labelpad=10)
    fig.colorbar(im, ax=ax)
    plt.yticks([])
    plt.xticks([])
    plt.suptitle('Mutation density in ' + gene + ' of SARS-CoV-2 virus', y=0.95)
    plt.savefig('./plots/density/mutation_density_' + gene + '.png', dpi=600)
    print('File saved in plots/density as mutation_density_{}.png'.format(gene))
    # plt.show()
    # print(label)

# colorbar1('M_gene')


# for i in genes:
#     colorbar1(i)
