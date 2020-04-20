import pandas as pd
import matplotlib.pyplot as plt

def date_histogram(file):
    import matplotlib.pyplot as plt
    df = pd.read_csv(file)
    df['date'] = df['date'].astype('datetime64')
    for line in df:
        if (line != 'date'):
            df.drop(labels=line,inplace= True,axis=1)
    #specify below if you want it grouped by days/weeks/months
    # df.groupby([df['date'].dt.year, df['date'].dt.month, df['date'].dt.week]).count().plot(kind="bar")
    df = df.groupby([df['date'].dt.year, df['date'].dt.month, df['date'].dt.week]).count()
    print(df)
    # plt.show()

import datetime
date_histogram('Genome_Data-X.csv')
def return_ratio(month):
    df = pd.read_csv('Genome_Data-X2.csv')
    df['date'] = df['date'].astype('datetime64')
    # df = df.sort_values(by='date', ascending=False)
    # a = str(df['date'].head(1)).split('-')
    df['filter'] = df.apply(lambda row: len(str(row.date).split('-')), axis=1)
    df = df[df['filter'] == 3]
    df = df[df['sum_N'] == 0]
    df.drop(labels='filter', axis=1, inplace=True)

    df['sum'] = df.apply(lambda row: row.E_gene + row.M_gene + row.N_gene + row.S_gene + row.orf1ab + row.ORF3a + row.ORF6 + row.ORF7 + row.ORF8 + row.ORF10, axis=1)
    leng = (len(df[df['date'].dt.week == month]))
    FirstVal = df[df['date'].dt.week == month]
    # SecondVal = df[df['date'].dt.week == (month + 1)]

    # FirstVal = df[df['date'].dt.week == month]
    # SecondVal = df[df['date'].dt.week == (month + 1)]

    # FirstVal = FirstVal[FirstVal['date'].dt.day <= 14]
    # print(len(FirstVal))

    # SecondVal = SecondVal[SecondVal['date'].dt.day <= 28]
    # SecondVal = SecondVal[SecondVal['date'].dt.day > 14]
    # print(len(SecondVal))
    # print(SecondVal)
    # df['Year'] = df.apply(lambda row: str(row.date).replace('\n', '-').split('-')[0].split(' ')[-1], axis=1)
    # df['Month'] = df.apply(lambda row: int(str(row.date).split('-')[1]), axis=1)
    # df['Day'] = df.apply(lambda row: int(str(row.date).split('-')[2]), axis=1)
    # FirstVal = df[df['Month'] == 2]
    # FirstVal = FirstVal[FirstVal['Day'] <= 14]

    # SecondVal = df[df['Month'] == 3]
    # SecondVal = SecondVal[SecondVal['Day'] <= 28]
    # SecondVal = SecondVal[SecondVal['Day'] > 14]

    Genes = ['E_gene', 'M_gene', 'S_gene', 'N_gene', 'orf1ab', 'ORF3a', 'ORF6', 'ORF7', 'ORF8', 'ORF10', 'whole_genome',
             'sum']

    vals = {}
    for gene in Genes:
        # print(gene)
        # print('1st')
        # print(FirstVal[gene].mean())
        # print('2nd')
        # print(SecondVal[gene].mean())
        # print('____')
        # vals[gene] = SecondVal[gene].mean() / (FirstVal[gene].mean() + 0.00000000001)
        vals[gene] = FirstVal[gene].mean()
        # print(vals[gene])
        # print('XXXXXXXX')

    # print('....')
    # print(vals['sum'])
    return vals['sum'], leng

b=[]
length=[]
from sklearn.linear_model import LinearRegression
for i in range (1,15):
    c,d = return_ratio(i)
    b.append(c)
    length.append(d)
import numpy as np
# print(b)

print(length)
print(sum(length))

X = np.asarray(range(1,len(b)+1))

linear_regressor = LinearRegression()  #create object for the class

b = np.asarray(b)
X= X.reshape(-1,1)
linear_regressor.fit(X, b)  # perform linear regression
Y_pred = linear_regressor.predict(X)  # make predictions



plt.scatter(X, b)
plt.plot(X, Y_pred, color='red')
plt.show()

# print(length)