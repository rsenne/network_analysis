#Before you run this code, you have to manually delete left and right "Medulla", "Pons" and "Midbrain" from every .csv file from LCT
#There are areas that are needed that start with "Medullary" or "Midbrain" that will still get dropped if those are included

#Import libraries
import pandas as pd
import os
import pickle as pkl

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
    df = df.loc[:,["name","acronym","density (cells/mm^3)","filename"]]
    return df

#Get each df for Control and ChR2 group and sort them by name
df_ChR2 = collect_csvs("/Users/kaitlyndorst/Desktop/ChR2_Large_Box/data").sort_values(by = ["name","filename"])
df_Control = collect_csvs("/Users/kaitlyndorst/Desktop/Control_Small_Box/data").sort_values(by = ["name","filename"])

#Create a hard-code list of terms that can drop entries in the dataframe
discard = ["background","layer","/Layer","5","6","Basic","Cerebrum","Cerebral","Isocortex","Cerebellar","Cerebellum","Thalamus",
           "Hypothalamus","Epithalamus","Interbrain","Pallidum","Striatum","Hindbrain","Brain stem","plate","pole","tectospinal pathway",
           "Lingula","II","III","IV-V","VI","VII","VIII","IX","X","Nucleus x","Nucleus y","Primary","Secondary","Supplemental",
           "Spinal","Accessory","dorsal part","ventral part","lateral part","medial part","anterior part","posterior part",
           "ventrolateral part","agranular part","ventromedial part","caudodorsal","rostroventral","capsular part",
           "parvicellular part","magnocellular part","Dorsal part","Ventral part","intermediate part","preoptic part","median part",
           "intermediate part","Interpeduncular nucleus,","Inferior colliculus,","division","like","lobule","visual","Abducens",
           "Ammon","Anterior area","olfactory","pretectal","posterior commisure","optic tract","Anterodorsal","Anteromedial",
           "Anteroventral nucleus","Arcuate","Supraoptic","Paraventricular hypothalamic","postrema","prostriata","Barrington",
           "ochlear","Crus","Culmen","Cuneate","Copula","Cortical amygdalar","Dorsal cochlear",
           "Dorsal column","vagus","peduncular area","auditory","Dorsal tegmental","Periventricular region","Ethmoid",
           "External","Facial","Forel","locculus","Gigantocellular","Gracile","Hippocamp","ypoglossal","lateral zone",
           "medial zone","olivary","salivatory","Infracerebellar","Interantero","leaflet","Intermedi","Intertrigeminal",
           "Lateral mammillary","Medial mammillary","Supramammillary","Tuberomammillary","Lateral posterior",
           "Lateral reticular","septal complex","vestibular","Latero","Linear","Magnocellular reticular","Medial accessory",
           "geniculate complex","Medullary","raphe nuclei","Midbrain,","Midline group","Motor nucleus","ambiguus","incertus",
           "Roller","lemniscus","solitary tract","trapezoid body","prepositus","Nucleus raphe","Paragigantocellular",
           "Parapyramidal","Parasolitary","Parataenial","Paratr","Parvicellular","Perireunesis","Peritrigeminal","Piriform",
           "Pontine","Posterior complex","Posterior intralaminar","limiting","Posterior triangular","Posterodorsal tegmental nucleus",
           "transition","Principal sensory","Retrohippocampal","Retroparafascicular","Subceruleus nucleus","Subcommissural",
           "Subgeniculate","Sublaterodorsal","Submedial","Suprage","Supratrigeminal","Taenia tecta","Tegmental reticular",
           "Ventral cochlear","Ventral group","Ventral posterolateral","Ventral posteromedial","Vestibular nuclei","Oculomotor",
           "Nucleus of the posterior","Postrhinal","Anterior group of the dorsal thalamus","Cajal","Precommissural",
           "Darkschewitsch","Supraoculomotor periaqueductal gray","Intralaminar nuclei","group of the dorsal","Geniculate group",
           ]


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
ChR2 = ChR2[["filename","name","acronym_Left","Bilateral Density (cells/mm^3)"]]

Control["Bilateral Density (cells/mm^3)"] = Control[["density (cells/mm^3)_Left","density (cells/mm^3)_Right"]].mean(axis = 1)
Control = Control[["filename","name","acronym_Left","Bilateral Density (cells/mm^3)"]]

ChR2["acronym_Left"] = ChR2["acronym_Left"].map(lambda x: x.replace('-L',""))
Control["acronym_Left"] = Control["acronym_Left"].map(lambda x: x.replace('-L',""))

#Use the pivot function to get the correct index that Ryan wants and just get the levels to have only the area names
ChR2 = ChR2.pivot(index = "filename", columns = "acronym_Left", values =["Bilateral Density (cells/mm^3)"])
ChR2.columns = ChR2.columns.get_level_values(1)
ChR2.reset_index(drop = True, inplace = True)

Control = Control.pivot(index = "filename", columns = "acronym_Left", values =["Bilateral Density (cells/mm^3)"])
Control.columns = Control.columns.get_level_values(1)
Control.reset_index(drop = True, inplace = True)

#Here are the ROIs by Allen Brain Group that we want to organize
ROIs = pd.read_csv("/Users/kaitlyndorst/Documents/GitHub/networkx/csv_files/ROIs.csv").loc[:,["Abbreviation","Allen Area"]]
ROIs_dict = dict(ROIs.values)

#Picking the dictionary
with open('Allen_Areas_dict.pickle','wb') as f:
    pkl.dump(ROIs_dict,f)

#Do not need any of this fluff, this is only for organizing the final .csv files based on Allen Anatomy
'''
ROIs = ROIs.sort_values("Allen Group Name").set_index("Abbreviation").T
cols = ROIs.columns.tolist()

#Need to reorganize the datframe by the index defined by our ROIs
ChR2 = ChR2[cols]
Control = Control[cols]'''

#Turn the final ChR2 and Control dfs into .csv files to send to Ryan
ChR2.to_csv("/Users/kaitlyndorst/Desktop/ChR2_Large_Box/ChR2_Large_Box.csv", index = False)
Control.to_csv("/Users/kaitlyndorst/Desktop/Control_Small_Box/Control_Small_Box.csv", index = False)