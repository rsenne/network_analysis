import pandas as pd
import os


def collapseCSVS(filepath):
    files = [filepath + '/' + f for f in os.listdir(filepath) if f.endswith('.csv')]
    df = pd.concat(map(pd.read_csv, files), axis=1)
    df = df.dropna(subset=['name'])
    return df


def getCounts(df):
    counts = df.groupby(['name', 'image'], as_index=False).size()
    counts = pd.pivot_table(counts, values='size', index='image', columns='name')
    df1 = counts.apply(lambda x: pd.Series(x.dropna().values))
    return df1


def split_and_sum(df):
    # [s.split(', ')[0] for s in df.T.index.values]
    d = df.T.groupby([s.split(', ')[0] for s in df.T.index.values]).sum().T
    return d


# df = collapseCSVS('/home/ryansenne/Data/Network/RData/N4_M2')
# counts = getCounts(df)
x = collapseCSVS('/home/ryansenne/Desktop/Network_Cells/Final')