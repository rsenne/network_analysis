#Before you run this code, you have to manually delete left and right "Medulla" and "Midbrain" from every .csv file from LCT
#There are areas that are needed that start with "Medullary" or "Midbrain" that will still get dropped if those are included

#Import libraries
import pandas as pd
import os

#Get the csv files and concatenate only the three columns you need in the analysis
#Yes Ryan, I know this is still a for-loop isntead of a list comp so bite me
def collect_csvs(filepath):
    files = [filepath + '/' + f for f in os.listdir(filepath) if f.endswith('.csv')]
    data = []
    for i in files:
        frame = (pd.read_csv(i))
        frame["filename"] = os.path.basename(i)
        data.append(frame)
    df = pd.concat(data, ignore_index= True)
    df = df.loc[:,["name","density (cells/mm^3)","filename"]]
    return df

#Get each df for Control and ChR2 group and sort them by name
df_Control = collect_csvs("/Users/kaitlyndorst/Desktop/Data_Analyses/Networks/Network Wrangling/Control/").sort_values(by = ["name","filename"])
df_ChR2 = collect_csvs("/Users/kaitlyndorst/Desktop/Data_Analyses/Networks/Network Wrangling/ChR2/").sort_values(by = ["name","filename"])

#Create a hard-code list of terms that can drop entries in the dataframe
discard = ["background","ayer","Cerebellar","Cerebellum","Cerebral","Thalamus","Hypothalamus","Pons","Brain stem",
           "Cochlear","plate","subplate","Epithalamus","Forel","A13","Hindbrain","Interbrain",
           "Interpeduncular nucleus,","septal complex","Lobule II","Lobule III","Lobules IV-V","mammilary nucleus",
           "Pretectal","Pallidum","Vestibular","retrorubral","Olfactory areas",
           "sland","Retrohippocampal", "Pons","Periventricular region","Somatomotor areas","Somatosensory areas","Striatum",
           "Mammilary body", "posterior complex","Visual areas","corpus callosum","cranial nerves","6","part","zone",
           "Ammon's","groups", "Basic","region","Crus","part","group","division","Auditory areas","-like","pathway",
           "ventricle","fasicle","root","matter","unassigned","direct tectospinal pathway","alveus","decussation",
            "general","system","formation","capsule","forceps","splenium","Cerebrum","Dorsal column nuclei","central nucleus",
           "dorsal nucleus","external nucleus","Isocortex","Perihypoglossal nuclei","Primary somatosensory area",
           "Thalamus,","nuclei","tectospinal pathway"]

#Create new dataframes after dropping the lines in the discard list
df_Control_deleted = df_Control[~df_Control.name.str.contains('|'.join(discard))]
df_ChR2_deleted = df_ChR2[~df_ChR2.name.str.contains('|'.join(discard))]

#Now split each Dataframe into left and right, strip the "left" and "right" strings in the columns
Split_cond = df_ChR2_deleted.name.str.contains("left")
Left_ChR2 = df_ChR2_deleted[Split_cond]
Right_ChR2 = df_ChR2_deleted[~Split_cond]

Split_cond_2 = df_Control_deleted.name.str.contains("left")
Left_Control = df_Control_deleted[Split_cond_2]
Right_Control = df_Control_deleted[~Split_cond_2]

#Strip the left and right labels in the names column
Left_ChR2["name"] = Left_ChR2["name"].map(lambda x: x.lstrip('left'))
Left_Control["name"] = Left_Control["name"].map(lambda x: x.lstrip('left'))

Right_ChR2["name"] = Right_ChR2["name"].map(lambda x: x.lstrip('right'))
Right_Control["name"] = Right_Control["name"].map(lambda x: x.lstrip('right'))

#After the strip and before the merge, drop all of the rows that have lowercase letters
#We do not need to include fiber tracts in the network analyses

Left_ChR2 = Left_ChR2[~Left_ChR2.name.str.islower()]
Right_ChR2 = Right_ChR2[~Right_ChR2.name.str.islower()]

Left_Control = Left_Control[~Left_Control.name.str.islower()]
Right_Control = Right_Control[~Right_Control.name.str.islower()]

#Merge each Left and Right

ChR2 = pd.merge(Left_ChR2,Right_ChR2,how = "inner",on = ["name","filename"],suffixes=("_Left","_Right"))
Control = pd.merge(Left_Control,Right_Control,how='inner',on = ["name","filename"],suffixes=("_Left","_Right"))

#Average the two density columns to get a bilateral reading from each animal
#Only take the three columns "filename","name", and the bilateral averages into the new dataframe

ChR2["Bilateral Density (cells/mm^3)"] = ChR2[["density (cells/mm^3)_Left","density (cells/mm^3)_Right"]].mean(axis = 1)
ChR2 = ChR2[["filename","name","Bilateral Density (cells/mm^3)"]]

Control["Bilateral Density (cells/mm^3)"] = Control[["density (cells/mm^3)_Left","density (cells/mm^3)_Right"]].mean(axis = 1)
Control = Control[["filename","name","Bilateral Density (cells/mm^3)"]]

#Use the pivot function to get the correct index that Ryan wants and just get the levels to have only the area names
ChR2 = ChR2.pivot(index = "filename", columns = "name", values =["Bilateral Density (cells/mm^3)"])
ChR2.columns = ChR2.columns.get_level_values(1)
ChR2.reset_index(drop = True, inplace = True)

Control = Control.pivot(index = "filename", columns = "name", values =["Bilateral Density (cells/mm^3)"])
Control.columns = Control.columns.get_level_values(1)
Control.reset_index(drop = True, inplace = True)

#Turn the final ChR2 and Control dfs into .csv files to send to Ryan
ChR2.to_csv("/Users/kaitlyndorst/Desktop/Data_Analyses/Networks/Network Wrangling/ChR2/ChR2_Large_Network.csv", index = False)
Control.to_csv("/Users/kaitlyndorst/Desktop/Data_Analyses/Networks/Network Wrangling/Control/Control_Large_Network.csv", index = False)
