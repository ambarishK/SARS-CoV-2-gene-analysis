import plotly.express as px
import pandas as pd
import matplotlib
import numpy as np
import os

genes = ['E_gene', 'M_gene', 'S_gene', 'N_gene', 'orf1ab', 'ORF3a', 'ORF6', 'ORF7', 'ORF8', 'ORF10']

df = pd.read_csv('Genome_Data-X3.csv')
df['date'] = pd.to_datetime(df['date'])
df['week_number'] = df['date'].dt.week
df = df[df['week_number'] < 30]
df['counts'] = 1
df = df[df['sum_N'] == 0]

def bubble(df=df, gene = 'sum'):
    df = df.groupby(['week_number', gene]).size().reset_index(name='counts')
    df['counts_sqrt'] = np.power((df['counts']), 0.5)
    week_number_sums = df.groupby('week_number')['counts_sqrt'].sum()

    df['week_sum'] = df['week_number'].map(week_number_sums)

    df['val_share'] = np.floor(df['counts_sqrt'] / df['week_sum'] * 100)
    df['% of week data'] = df['val_share']
    df['Week'] = df['week_number']
    df['Number of mutations'] = df[gene]
    df['Number of analyzed genomes'] = df['counts']

    fig = px.scatter(df, x='Week', y='Number of mutations', size='% of week data',
                     hover_data=['Number of analyzed genomes'])
    if (all(v == 0 for v in df['Number of mutations'].tolist())):
    #     fig.update_yaxes(yaxis = dict(
    #     range=(0, 1, 2),
    #     constrain='domain'
    # ))
        fig.update_yaxes(rangemode="tozero")
    fig.update_layout(
        title='Mutation in ' + gene,
        hoverlabel=dict(
            bgcolor="white",
            font_size=16,
            font_family="Rockwell",
        ),
    )
    name = 'div_bubble_chart_' + gene + '.html'
    fig.write_html('./plots/html/' + name, full_html=False, include_plotlyjs='cdn')
    print('File saved in directory plots/html/ as ' + name)
    fig.show()

# for gene in genes:
#     bubble(df, gene)
print(os.path.dirname(os.path.abspath(__file__)))

if __name__ == "__main__":
    # bubble(df, 'ORF10')
    pass