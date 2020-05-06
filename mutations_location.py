import ast
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import statistics

pd.set_option("display.max_rows", None, "display.max_columns", None)
x = '[ "A","B","C" , " D"]'
x = ast.literal_eval(x)
DF = pd.read_csv('Genome_Data-X.csv')
DF['date'] = pd.to_datetime(DF['date']) #converting data column to data type
#DF = DF.loc[DF['Nuc.Completeness'] == 'Complete'] #only using rows which have complete genome data

DF.sort_values(by ='date',ascending = False, inplace= True)
DF = DF.head((len(DF)-6)) #deleting outlier data
# print(len(DF))
DF = DF[DF['N_gene']==0]
# print(len(DF))
string = 'S_gene_translation_changes'
# DF['orf1ab_mutations'] = DF.apply(lambda
#                                row: (ast.literal_eval(row.orf1ab_mutations)),
#                            axis=1)
# hist_data = []
# for cols in DF:
#     if '_mutations' in  cols or '_changes' in cols:
#         for sth in DF[cols]:
#             sth = ast.literal_eval(sth)


def change_format(list1):
    max_val = max(list1)
    value_counter = []
    for i in range(max_val + 10):
        value_counter.append(list1.count(i))
    return value_counter


# plt.style.use('ggplot')
# fig, ax = plt.subplots(nrows=20,ncols=1, sharex=True, sharey=True)
mut = ['E_gene_mutations', 'M_gene_mutations', 'S_gene_mutations', 'N_gene_mutations', 'orf1ab_mutations', 'ORF3a_mutations',
       'ORF6_mutations', 'ORF7_mutations', 'ORF8_mutations', 'ORF10_mutations']
trans = ['E_gene_translation_changes', 'M_gene_translation_changes', 'S_gene_translation_changes', 'N_gene_translation_changes',
         'orf1ab_translation_changes', 'ORF3a_translation_changes', 'ORF6_translation_changes', 'ORF7_translation_changes', 'ORF8_translation_changes',
         'ORF10_translation_changes']


def filter_country(country='Poland', dataframe=DF):
    dataframe['country'] = dataframe['name'].apply(lambda x: str(x).split('/')[-3])
    dataframe2 = dataframe[dataframe['country'] == country]
    print('Number of genomes in {} is {}, {}% of whole'.format(country, len(dataframe2),
                                                               round(len(dataframe2) / len(dataframe), 3)))
    return dataframe2


def get_lists_colorbar():
    lists = []
    for cols in mut:
        hist_data=[]
        for location in DF[cols]:
            location = ast.literal_eval(location)
            if len(location) != 0:
                for list in location:
                    hist_data.append(list[0])
        new_list = change_format(hist_data)
        lists.append(new_list)
    return lists

def llwrite(list,filename='mutation_density.pg'):
    with open(filename, 'w') as f:
        for _list in list:
            f.write('#' + '\n')
            for inty in _list:
                f.write(str(inty) + '\n')

# llwrite(get_lists_colorbar(),'mutation_density2.pg')

def llread(filename='mutation_density2.pg'):
    with open(filename, 'r') as f:
        result=[]
        i=0
        for line in f:
            if line == '#\n':
                if i==0:
                    i+=1
                else:
                    result.append(a)
                a=[]
            else:
                a.append(int(line[:-1]))
        result.append(a)
        return result


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
    # a = []
    for list in data:
        list = del_zscore(list, 10)
        print(mut[i])
        label = mut[i][:-10]
        arr = np.array(list)
        arr = np.stack((arr, arr), axis=0) #thicc1
        # arr = np.concatenate((arr, arr), axis=0)
        # arr = np.concatenate((arr, arr), axis=0)
        # arr = np.concatenate((arr, arr), axis=0) #thicc4
        # ax.subplot(20,1,t)
        im = axs[i].imshow(arr,aspect='auto',cmap='plasma')
        axs[i].set_ylabel(label, rotation=0, labelpad=25, y=0.3)
        # axs[i].xaxis.labelpad=100
        # axs[i].set_yticklabels([])
        # plt.ylabel(cols[:-10])
        # axs[i].axis('off')
        axs[i].axes.get_xaxis().set_ticks([])
        axs[i].axes.get_yaxis().set_ticks([])
        i+=1
    # fig.colorbar(im, ax=axs)
    # plt.yticks([])
    # plt.xticks([])
    plt.suptitle('Mutation density in genes of SARS-CoV-2 virus',y=0.95)
    plt.savefig('./plots/density/mutation_density2.png',dpi=600)
    plt.show()


def icolorbar():
    import plotly.express as px
    data = [[1, 25, 30, 50, 1], [20, 1, 60, 80, 30], [30, 60, 1, 5, 20]]
    fig = px.imshow(data,
                    labels=dict(x="Day of Week", y="Time of Day", color="Productivity"),
                    x=['Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday'],
                    y=['Morning', 'Afternoon', 'Evening']
                    )
    fig.update_xaxes(side="top")
    fig.show()


# icolorbar()

# colorbar()


def colorbar1(gene):
    i=0
    for gene1 in mut:
        if gene1.startswith(gene):
            break
        i+=1
    data=llread()
    fig, ax = plt.subplots()
    label = mut[i][:-10]
    print(data[i])
    arr = np.array(del_zscore(data[i],20))
    # print(len(arr))
    # print((arr))
    print('Mean number of mutations is {}'.format(round(statistics.mean(map(float, arr))),2))
    # indexNames = dfObj[(dfObj['Age'] >= 30) & (dfObj['Country'] == 'India')].index
    arr = np.stack((arr,arr),axis=0) #thicc1
    im = ax.imshow(arr,extent=(0,1,0,0.1),cmap='plasma')
    ax.axes.get_yaxis().set_ticks([])
    ax.set_xlabel('Length of ' + label + ' (' + str(len(data[i])) + ')', labelpad=10)
    fig.colorbar(im, ax=ax)
    plt.yticks([])
    plt.xticks([])
    plt.suptitle('Mutation density in ' + label +  ' of SARS-CoV-2 virus', y=0.95)
    plt.savefig('./plots/density/mutation_density_' + label + '.png', dpi=600)
    plt.show()
    # print(label)

# colorbar1('N_gene')


for i in mut:
    colorbar1(i)