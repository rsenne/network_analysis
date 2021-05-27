import pandas as pd
import os


def collapseCSVS(filepath):
    files = [filepath + '/' + f for f in os.listdir(filepath) if f.endswith('.csv')]
    df = pd.concat(map(pd.read_csv, files))
    df = df.dropna(subset=['name'])
    return df


def getCounts(df):
    counts = df.groupby(['name', 'image'], as_index=False).size()
    counts = pd.pivot_table(counts, values='size', index='image', columns='name')
    df1 = counts.apply(lambda x: pd.Series(x.dropna().values))
    return df1


df = collapseCSVS('/home/ryansenne/Data/Network/RData/N4_M2')
counts = getCounts(df)
