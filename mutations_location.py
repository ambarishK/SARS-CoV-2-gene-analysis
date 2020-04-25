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
str = 'S_gene_translation_changes'
# DF['orf1ab_mutations'] = DF.apply(lambda
#                                row: (ast.literal_eval(row.orf1ab_mutations)),
#                            axis=1)
# hist_data = []
# for cols in DF:
#     if '_mutations' in  cols or '_changes' in cols:
#         for sth in DF[cols]:
#             sth = ast.literal_eval(sth)
plt.style.use('ggplot')
for cols in DF:
    hist_data = []
    if '_mutations' in cols or '_changes' in cols:
        for location in DF[cols]:
            location = ast.literal_eval(location)
            if len(location) !=0:
                for list in location:
                     hist_data.append(list[0])
        plt.hist(hist_data, bins =max(hist_data)+10)
        plt.title(cols)
        plt.show()
        break
        plt.savefig('plots/distributions/' + cols + '.png')


