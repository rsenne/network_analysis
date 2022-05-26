import pandas as pd
import numpy as np
import pingouin as pg

def percentage_change(data_mean, norm_group = 'SAL'):
    percentage_change_dict = {}
    if norm_group not in data_mean.keys():
        print("The normalization group should be a key of the data_mean dictionary")
        exit()
    else:
        data_norm = data_mean[norm_group].values[0]
        cols_norm = list(data_mean[norm_group])
    for e in [exp for exp in data_mean.keys() if exp != norm_group]:
        data = data_mean[e].values[0]
        cols = list(data_mean[norm_group])
        if cols_norm != cols:
            print('The order of the names of the dataframe is different and needs to be rearanged')
            continue
        per_change = np.divide(np.subtract(data, data_norm), data_norm)*100
        dm = pd.DataFrame(per_change, index = cols, columns = ['Percentage of change'])
        dm = dm.sort_values(by = 'Percentage of change', ascending=False)
        #cols_1 = dm.iloc[0]
        #dm = dm[1:]
        #dm.columns = cols_1
        percentage_change_dict[e] = [per_change, cols, dm]
    return percentage_change_dict

def tstud(data, comp_group, threshold = 0.05):
    stats_df = []
    significant_areas = []
    for area, v in data.items():
        x = np.array(v)
        y = np.array(comp_group[area])
        n = np.concatenate([x,y])
        normality = pg.normality(n)
        param = normality['normal'][0]
        if param:
            df = pg.ttest(x,y)
            stats_df.append([area,df['p-val'][0]])
            if df['p-val'][0] <= threshold:
                significant_areas.append(area)
        else:
            df = pg.mwu(x, y)
            stats_df.append([area, df['p-val'][0]])
            if df['p-val'][0] <= threshold:
                significant_areas.append(area)
    stats_df = pd.DataFrame(stats_df, columns=['area', 'p-value'])
    stats_df = stats_df.set_index('area')
    return stats_df, significant_areas
