# -*- coding: utf-8 -*-
"""Hackathon_HackTheCrisis

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1ESbRMgDbY6NNdpENv_TFpuRjw5DmoVUW
"""

# Commented out IPython magic to ensure Python compatibility.
from google.colab import drive 

drive.mount('/content/drive')
# %cd /content/drive/My Drive/Colab Notebooks

import pandas as pd
DF = pd.read_csv('Genome_Data.csv')
DF['Date'] = pd.to_datetime(DF['Date']) #converting data column to data type
DF = DF.loc[DF['Nuc.Completeness'] == 'Complete'] #only using rows which have complete genome data

DF.sort_values(by ='Date',ascending = False, inplace= True)
DF = DF.head((len(DF)-6)) #deleting outlier data

import matplotlib.pyplot as plt
import matplotlib
import numpy as np

def chart(gene,l_n = 0): #function which makes logarythmic/normal plot based on a column from dataframe
  list_of_datetimes = DF['Date'].tolist()
  dates = matplotlib.dates.date2num(list_of_datetimes)
  values = DF[gene].tolist()
  if l_n == 'l':
    values = np.log(values)
  plt.plot_date(dates, values, color='red', alpha = 0.2)
  plt.title('Number of mutations in whole genome of SARS-CoV-2', fontsize = 14.5)
  plt.ylabel('Log Scale, no. mutations')
  plt.xlabel('Time')
  plt.xticks(rotation=35)
  plt.tight_layout()
  plt.savefig('SARS-CoV-2_Mutations_in_time.png')
chart('whole_genome','l')
