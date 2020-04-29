import ast
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
pd.set_option("display.max_rows", None, "display.max_columns", None)
x = '[ "A","B","C" , " D"]'
x = ast.literal_eval(x)
DF = pd.read_csv('Genome_Data-X.csv')
DF['date'] = pd.to_datetime(DF['date']) #converting data column to data type
#DF = DF.loc[DF['Nuc.Completeness'] == 'Complete'] #only using rows which have complete genome data

DF.sort_values(by ='date',ascending = False, inplace= True)
DF = DF.head((len(DF)-6)) #deleting outlier data
string = 'S_gene_translation_changes'
# DF['orf1ab_mutations'] = DF.apply(lambda
#                                row: (ast.literal_eval(row.orf1ab_mutations)),
#                            axis=1)
# hist_data = []
# for cols in DF:
#     if '_mutations' in  cols or '_changes' in cols:
#         for sth in DF[cols]:
#             sth = ast.literal_eval(sth)
def change_format(list):
    max_val = max(list)
    value_counter = []
    for i in range(max_val + 10):
        value_counter.append(list.count(i))
    return value_counter
# plt.style.use('ggplot')
# fig, ax = plt.subplots(nrows=20,ncols=1, sharex=True, sharey=True)
mut = ['E_gene_mutations', 'M_gene_mutations', 'S_gene_mutations', 'N_gene_mutations', 'orf1ab_mutations', 'ORF3a_mutations',
       'ORF6_mutations', 'ORF7_mutations', 'ORF8_mutations', 'ORF10_mutations']
trans = ['E_gene_translation_changes', 'M_gene_translation_changes', 'S_gene_translation_changes', 'N_gene_translation_changes',
         'orf1ab_translation_changes', 'ORF3a_translation_changes', 'ORF6_translation_changes', 'ORF7_translation_changes', 'ORF8_translation_changes',
         'ORF10_translation_changes']
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
            for int in _list:
                f.write(str(int) + '\n')

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
    plt.savefig('./plots/density/mutation_density_' + label + '.png', dpi=600)
    plt.show()
    # print(label)

colorbar1('N_gene')

# for i in mut:
#     colorbar1(i)



# for cols in DF:
#     hist_data = []
#     if '_mutations' in cols or '_changes' in cols:
#         for location in DF[cols]:
#             location = ast.literal_eval(location)
#             if len(location) !=0:
#                 for list in location:
#                      hist_data.append(list[0])
#         new_list = change_format(hist_data)

        # arr1 = np.array(new_list)
        # arr = np.stack((arr1,arr1),axis=0)
        # arr = np.concatenate((arr, arr), axis=0)
        # arr = np.concatenate((arr, arr), axis=0)
        # arr = np.concatenate((arr, arr), axis=0)
        # plt.subplot(20,1,i)
        # ax.imshow(arr)
        # ax.colorbar()
        # plt.ylabel(cols[:-10])
        # i+=1
        # print(i)
        # plt.show()
        # print(arr)
        # break
        # plt.savefig('plots/distributions/' + cols + '.png')
# print(columny)
# plt.colorbar()
# plt.show()