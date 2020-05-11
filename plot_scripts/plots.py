import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import numpy as np


DF = pd.read_csv('Genome_Data-X3.csv')
#converting data column to data type
#DF = DF.loc[DF['Nuc.Completeness'] == 'Complete'] #only using rows which have complete genome data

DF.sort_values(by ='date',ascending = False, inplace= True)
DF = DF.head((len(DF)-6)) #deleting outlier data
print(DF.tail(6)['date'])
#DF['date'] = pd.to_datetime(DF['date'])



def chart(gene,l_n, plot_color = 'blue'): #function which makes logarythmic/normal plot based on a column from dataframe


  if 'whole' in gene:
      temp_df = DF
  elif '_protein' in gene:
      temp_df = DF[DF[gene.split('_p')[0] + '_invalid_nucleotides'] == 0]
  else:temp_df = DF[DF[gene+'_invalid_nucleotides'] == 0]
  list_of_datetimes = temp_df['date'].tolist()
  dates = matplotlib.dates.date2num(list_of_datetimes)
  values = temp_df[gene+'_distance'].tolist()
  if l_n == 'l':
    values = np.log(values)


  print(values)
  x = plt.plot_date(dates, values, color=plot_color)

  plt.title( gene +' of SARS-CoV-2', loc='left',fontsize = 10)
  plt.ylabel('Mutations', fontsize = 10)
  plt.xlabel('Time')
  plt.xlim('2019-12-30', '2020-04-30')

  plt.xticks(fontsize = 10, rotation=35)
  plt.tight_layout()
  #plt.savefig(gene + '.png')



Genes = [ 'E_gene', 'M_gene', 'S_gene', 'N_gene','orf1ab', 'ORF6', 'ORF7', 'ORF8', 'ORF10','ORF3a']
a=1

for gene in Genes:

    plt.subplot(1,2,1)
    chart(gene, 'e')
    plt.subplot(1, 2, 2)
    chart(gene + '_protein', 'e', 'red')
    plt.savefig('plots/'+ gene + '.png')
    # plt.show()

