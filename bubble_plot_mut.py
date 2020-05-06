import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import numpy as np

DF = pd.read_csv('Genome_Data-X.csv')
DF['date'] = pd.to_datetime(DF['date'])
DF['Week_Number'] = DF['date'].dt.week

# DF.sort_values(by='Week_Number', ascending=False, inplace=True)
# DF = DF.head((len(DF) - 6))  # deleting outlier data
DF['counter'] = 1
values = DF['orf1ab'].tolist()
dates = DF['Week_Number'].tolist()
time = []
mutations_list = []
counter = []

for time_step in list(set(dates)):
    temp_DF = DF[DF['date'] == time_step]
    print(time_step)
    for mutations_count in values:
        time.append(time_step)
        mutations_list.append(mutations_count)
        # counter.append(len(temp_DF[temp_DF['orf1ab'] == mutations_count]))

plt.scatter(time, mutations_list)
plt.show()