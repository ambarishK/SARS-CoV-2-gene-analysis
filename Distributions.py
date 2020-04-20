import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# DF = pd.read_csv('Genome_Data-25_03_2020.csv')
# DF['Date'] = pd.to_datetime(DF['Date']) #converting data column to data type
# # DF = DF.loc[DF['Nuc.Completeness'] == 'Complete'] #only using rows which have complete genome data
#
# DF.sort_values(by ='Date',ascending = False, inplace= True)
# DF = DF.head((len(DF)-6)) #deleting outlier data

def histograms(): #function plotting distribution of mutation and comparing it to geometric distribution
    for gene in DF:
        if(DF.dtypes[gene] == np.int64 and 'translation' not in gene):
            a = DF[gene].tolist()
            biny = 100
            outliers_perc = 0.6
            maxi = np.percentile(a,100-outliers_perc)
            a[:] = (x for x in a if (x<maxi))
            try: #some distribution are quite different and has mean that doesnt work there, so we check if it makes sense (is 0<p<1 in geom)
                b= np.random.geometric(1/(np.mean(a)+1),size = 5000)-1
            except:
                break
            bina = np.arange(0,maxi,maxi/biny)
            plt.ylabel('% of samples')
            plt.suptitle('Distribution of mutations number in '+ gene)
            plt.hist(a, bins = bina, weights=np.ones(len(a)) / len(a)*100, alpha=0.4, label=gene, color='green', edgecolor='black', linewidth=1.2)
            plt.hist(b, bins=bina, weights=np.ones(len(b)) / len(b) * 100, alpha=0.4, label='geometric',color='red', edgecolor='black', linewidth=1.2)
            plt.legend(loc='upper right')
            plt.savefig('plots/histograms/Hist_of_' + gene + '_mutations_geo.png') #zapisywanie pliku
            plt.show()

def protein_lengths_n(fp):  # funkcja ktora na podstawie pliku fasta liczy rozklad jego N
    seq = []
    with open(fp) as file:
        for line in file:
            line = line.rstrip()
            if line.startswith('>'):  # liczba aminokwasow
                seq.append('')
            else:
                seq[-1] = '{}{}'.format(seq[-1], line)

    numbers = []
    for item in seq:
        numbers.append(len(item.split('N')) - 1)
    return numbers

N_genes = ['E_gene_N', 'M_gene_N', 'S_gene_N', 'N_gene_N', 'orf1ab_N', 'ORF3a_N', 'ORF6_N', 'ORF7_N', 'ORF8_N', 'ORF10_N']
# N_genes = ['ORF8_N']

def n_histogram(genes=N_genes):
    df = pd.read_csv('Genome_Data-X_N.csv')
    for gene in genes:
        data1 = df[gene].tolist()

        #modification of data
        minimum = 1
        outliers_perc = 0
        maximum = np.percentile(data1,100-outliers_perc)
        print(maximum)
        data = [x for x in data1 if (x>=minimum and x<=maximum)]

        biny = 100
        zeros = data1.count(minimum-1)
        zeros_p = round((zeros / len(data1) * 100), 2)

        bina = np.arange(start=0,stop=maximum+maximum/biny,step=(maximum/biny))
        plt.hist(data, bins= bina,weights=np.ones(len(data)) / len(data)*100)

        plt.suptitle('% Distribution of N in '+gene[:-2] +', '+ str(zeros_p) + '% is 0')
        plt.ylabel('% of genomes')
        plt.xlabel('Number of N | Total no. of genomes plotted is ' + str(len(data)) + ' out of ' + str(len(data1)))
        plt.savefig('Distribution of N in '+ gene[:-2])
        # plt.show()
        plt.clf()


n_histogram()

# a = protein_lengths_n('sequences_2.fasta') #function which shows distribution of the number of N
# b= []
# for item in a:
#   if (item <= 20):
#     pass
#   else:
#     b.append(item)
# print('Procent genomów z liczbą N równą 0 to:')
# print((1-len(b)/len(a))*100)
# a = b
# biny = 100
# # outliers_perc = 5
# maxi = np.percentile(a,85)
# bina = np.arange(0,maxi,maxi/biny)
# plt.hist(a, bins= bina,weights=np.ones(len(a)) / len(a)*100)
# plt.hist(a, bins= biny)