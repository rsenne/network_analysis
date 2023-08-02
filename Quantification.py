import pandas as pd
import numpy as np
import pingouin as pg

<<<<<<< Updated upstream
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
=======
def percentage_change(data_mean, norm_group = 'SAL', persex = False):
    percentage_change_dict = {}
    if not persex:
        if norm_group not in data_mean.keys():
            print("The normalization group should be a key of the data_mean dictionary")
        else:
            data_norm = data_mean[norm_group].values[0]
            cols_norm = list(data_mean[norm_group])
        for e in [exp for exp in data_mean.keys() if exp != norm_group]:
            data = data_mean[e].values[0]
            cols = list(data_mean[e])
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
    else:
        if norm_group not in data_mean.keys():
            print("The normalization group should be a key of the data_mean dictionary")
        for e in [exp for exp in data_mean.keys() if exp != norm_group]:
            for i,s in enumerate(list(data_mean[e].comparison)):
                data2 = data_mean[e].drop('comparison', axis = 1)
                data = data2.values[i]
                data_norm = data_mean[norm_group].drop('comparison', axis = 1).values[i]
                cols = list(data2)
                cols_norm = list(data_mean[norm_group].drop('comparison', axis = 1))
                if cols_norm != cols:
                    print('The order of the names of the dataframe is different and needs to be rearanged')
                    continue
                per_change = np.divide(np.subtract(data, data_norm), data_norm) * 100
                dm = pd.DataFrame(per_change, index=cols, columns=['Percentage of change'])
                dm = dm.sort_values(by='Percentage of change', ascending=False)
                # cols_1 = dm.iloc[0]
                # dm = dm[1:]
                # dm.columns = cols_1
                percentage_change_dict['_'.join([e,s])] = [per_change, cols, dm]

>>>>>>> Stashed changes
    return percentage_change_dict

def tstud(data, comp_group, threshold = 0.05):
    stats_df = []
    significant_areas = []
    for area, v in data.items():
<<<<<<< Updated upstream
=======
        if area in ['condition', 'name']:
            continue
>>>>>>> Stashed changes
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
<<<<<<< Updated upstream
=======


def twanova(combined_df, threshold = 0.05):
    stats_df = []
    posthocs = []
    significant_areas = []
    for area, v in combined_df.drop(['condition', 'name',
       'sex'], axis = 1).items():
        n = np.array(v)
        normality = pg.normality(n)
        param = normality['normal'][0]
        df = pg.anova(data = combined_df, dv = area, between = ['sex', 'condition'])
        df['area'] = area
        stats_df.append(df)
        if sum(threshold > np.array(list(df['p-unc']))) <= 0:
            significant_areas.append(area)
        ph = pg.pairwise_ttests(data = combined_df, dv = area, between = ['sex', 'condition'], parametric=param,
                                alpha=0.05, padjust='bonf', effsize='eta-square',
                                correction='auto',
                                nan_policy='listwise', return_desc=False, interaction=True,
                                within_first=False
                                )
        ph2 = pg.pairwise_ttests(data = combined_df, dv = area, between = ['condition', 'sex'], parametric=param,
                                alpha=0.05, padjust='bonf', effsize='eta-square',
                                correction='auto',
                                nan_policy='listwise', return_desc=False, interaction=True,
                                within_first=False
                                )
        ph['area'] = area
        ph2['area'] = area
        posthocs.append(ph)
        posthocs.append(ph2)
    stats_df = pd.concat(stats_df)
    posthocs = pd.concat(posthocs)
    return stats_df, significant_areas, posthocs


def percentage_change_sex(data_mean, norm_group = 'Female'):
    percentage_change_dict_sex = {}
    for e in data_mean.keys():
        for s in [sex for sex in list(data_mean[e].comparison) if sex != norm_group and sex != 'both']:
            data = data_mean[e].query('comparison==@s').drop('comparison', axis=1).values[0]
            data_norm = data_mean[e].query('comparison==@norm_group').drop('comparison', axis=1).values[0]
            cols = list(data_mean[e].drop('comparison', axis=1))
            per_change = np.divide(np.subtract(data, data_norm), data_norm)*100
            dm = pd.DataFrame(per_change, index = cols, columns = ['Percentage of change'])
            dm = dm.sort_values(by = 'Percentage of change', ascending=False)
            #cols_1 = dm.iloc[0]
            #dm = dm[1:]
            #dm.columns = cols_1
            percentage_change_dict_sex['_'.join([e,s])] = [per_change, cols, dm]
>>>>>>> Stashed changes
