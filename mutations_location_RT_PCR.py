import ast
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import numpy as np

pd.set_option("display.max_rows", None, "display.max_columns", None)
x = '[ "A","B","C" , " D"]'
x = ast.literal_eval(x)
DF = pd.read_csv('Genome_Data-X2.csv')
DF['date'] = pd.to_datetime(DF['date']) #converting data column to data type
#DF = DF.loc[DF['Nuc.Completeness'] == 'Complete'] #only using rows which have complete genome data

DF.sort_values(by ='date',ascending = False, inplace= True)
DF = DF.head((len(DF)-6)) #deleting outlier data
primers = ['N_China_F', 'N_China_R', 'N_China_P', 'N_HongKong_F', 'N_HongKong_R', 'N_HongKong_P', 'N_Japan_F', 'N_Japan_R', 'N_Japan_P', 'N_Thailand_F', 'N_Thailand_R', 'N_USA1_F', 'N_USA1_R', 'N_USA2_F', 'N_USA2_R', 'N_USA2_P', 'N_USA3_F', 'N_USA3_R', 'N_USA3_P', 'RdRp_China_F', 'RdRp_China_R', 'RdRp_China_P', 'RdRp_Germany_F', 'RdRp_Germany_R', 'RdRp_Germany_P', 'RdRp_HongKong_F', 'RdRp_HongKong_R', 'RdRp_HongKong_P']

for prime in primers:
    DF= DF[DF[prime+'_N']==0]

def filter_country(country = 'Poland', dataframe=DF):
    dataframe['country'] = dataframe['name'].apply(lambda x: str(x).split('/')[-3])
    dataframe2 = dataframe[dataframe['country'] == country]
    print('Number of genomes in {} is {}, {}% of whole'.format(country, len(dataframe2), round(len(dataframe2)/len(dataframe),3)))
    return dataframe2

filter_country()
# DF['orf1ab_mutations'] = DF.apply(lambda
#                                row: (ast.literal_eval(row.orf1ab_mutations)),
#                            axis=1)
# hist_data = []
# for cols in DF:
#     if '_mutations' in  cols or '_changes' in cols:
#         for sth in DF[cols]:
#             sth = ast.literal_eval(sth)

def change_format(list):
    if(len(list) != 0):
        max_val = max(list)
        value_counter = []
        for i in range(max_val + 10):
            value_counter.append(list.count(i))
        return value_counter
    else: return [-1]
# plt.style.use('ggplot')
# fig, ax = plt.subplots(nrows=20,ncols=1, sharex=True, sharey=True)
mut = [i +'_mutations' for i in primers]
trans = ['E_gene_translation_changes', 'M_gene_translation_changes', 'S_gene_translation_changes', 'N_gene_translation_changes',
         'orf1ab_translation_changes', 'ORF3a_translation_changes', 'ORF6_translation_changes', 'ORF7_translation_changes', 'ORF8_translation_changes',
         'ORF10_translation_changes']

# def get_lists_colorbar():
#     lists = []
#     for cols in mut:
#         hist_data=[]
#         for location in DF[cols]:
#             location = ast.literal_eval(location)
#             if len(location) != 0:
#                 #IF name == USA2 etc, the same for other countries
#                 if location[0][1]=='G' and location[0][2] == 'A':
#                     location.remove(location[0])
#                 # if location[1][1] == 'A' and location[1][2] == 'G':
#                 #     location.remove(location[1])
#                 for list in location:
#                     if list[1] =='Y' and (list[2]=='C' or list[2]=='T'):
#                         print(list[2])
#                     elif list[1]=='R' and(list[2]=='A' or list[2]=='G'):
#                         print(list[2])
#
#                     else:
#                         hist_data.append(list[0])
#         new_list = change_format(hist_data)
#         lists.append(new_list)
#     return lists


def get_lists_colorbar():
    lists = []
    for cols in mut:
        hist_data = []
        for location in DF[cols]:
            location = ast.literal_eval(location)

            if len(location) != 0:
                # IF name == USA2 etc, the same for other countries

                for list in location:

                    if 'N_HongKong_P' in cols:
                        if list[0] == 9 and location[0][2] == 'C' and location[0][1] == 'T':
                            location.remove(list)

                        print(location)

                    if list[1] == 'Y' and (list[2] == 'C' or list[2] == 'T'):
                        print(list[2])
                    elif list[1] == 'R' and (list[2] == 'A' or list[2] == 'G'):
                        print(list[2])

                    else:
                        hist_data.append(list[0])

        new_list = change_format(hist_data)
        lists.append(new_list)
    return lists



# get_lists_colorbar()

def llwrite(list,filename='mutation_density.pg'):
    with open(filename, 'w') as f:
        for _list in list:
            f.write('#' + '\n')
            for int in _list:
                f.write(str(int) + '\n')

# llwrite(get_lists_colorbar())

def llread(filename='mutation_density.pg'):
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
        return(result)


def colorbar():
    i = 0
    data = llread()
    fig, axs = plt.subplots(10,sharey=True,sharex=False)
    # a = []
    for list in data:
        label = mut[i][:-10]
        arr = np.array(list)
        arr = np.stack((arr,arr),axis=0) #thicc1
        # arr = np.concatenate((arr, arr), axis=0)
        # arr = np.concatenate((arr, arr), axis=0)
        # arr = np.concatenate((arr, arr), axis=0) #thicc4
        # ax.subplot(20,1,t)
        im = axs[i].imshow(arr,aspect='auto',cmap='plasma')
        axs[i].set_ylabel(label, rotation=0, labelpad=25)
        # axs[i].xaxis.labelpad=100
        # axs[i].set_yticklabels([])
        # plt.ylabel(cols[:-10])
        # axs[i].axis('off')
        axs[i].axes.get_xaxis().set_ticks([])
        axs[i].axes.get_yaxis().set_ticks([])
        i+=1
    fig.colorbar(im, ax=axs)
    # plt.yticks([])
    # plt.xticks([])
    plt.suptitle('Mutation density in genes of SARS-CoV-2 virus',y=0.95)
    plt.savefig('./plots/density/mutation_density.png')
    plt.show()

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
    arr = np.array(data[i])
    # print(len(arr))
    print(max(arr))
    arr = np.stack((arr,arr),axis=0) #thicc1
    # arr = np.concatenate((arr, arr), axis=0)
    # arr = np.concatenate((arr, arr), axis=0)
    # arr = np.concatenate((arr, arr), axis=0) #thicc4
    # ax.subplot(20,1,t)
    im = ax.imshow(arr,extent=(0,1,0,0.1),cmap='plasma')
    # ax.set_ylabel(label, rotation=0, labelpad=25)
    # axs[i].xaxis.labelpad=100
    # axs[i].set_yticklabels([])
    # plt.ylabel(cols[:-10])
    # axs[i].axis('off')
    # ax.axes.get_xaxis().set_ticks([])
    ax.axes.get_yaxis().set_ticks([])
    ax.set_xlabel('Length of ' + label, labelpad=10)
    fig.colorbar(im, ax=ax)
    # plt.yticks([])
    plt.xticks([])
    plt.suptitle('Mutation density in ' + label +  ' of SARS-CoV-2 virus', y=0.95)
    plt.savefig('./plots/density/mutation_density_' + label + '.png', dpi=300)
    # plt.show()
    # print(label)

# for prime in primers:
#     colorbar1(prime)

# for i in mut:
#     colorbar1(i)