import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors

sns.set_theme(style="whitegrid", palette="cubehelix")

subspecies_with_colors = pd.read_csv("Data/subspecies_with_colors.csv")
mgcsts_samples = pd.read_csv("Data/mgCSTs.samples.df.csv")
mgCSTs_sort = pd.read_csv("Data/mgCSTs.sort.df.csv")
projects = pd.read_csv("Data/VIRGO2_projects.csv")


st.sidebar.subheader("Parameters")
deepsplit = st.sidebar.slider(label="deepSplit", min_value=0, max_value=4)
mincluster = st.sidebar.slider(label="minClusterSize", min_value=10, max_value=50)


data = mgcsts_samples[(mgcsts_samples['deepSplit'] == deepsplit) & (mgcsts_samples['minClusterSize'] == mincluster)]
data = data.reset_index(drop = True)
data = data.merge(projects, on = "sampleID", how = "left")


data2 = mgCSTs_sort[(mgCSTs_sort['deepSplit'] == deepsplit) & (mgCSTs_sort['minClusterSize'] == mincluster)]
data2 = data2.reset_index(drop = True)

count_sample = []
for element in data2['dtc'].values :
    count_sample.append(data.groupby(['dtc']).count()['sampleID'][element])
data2['count_sample'] = count_sample

color = []
for element in data2['domTaxa'].values :
    a = subspecies_with_colors[subspecies_with_colors['Subspecies'].apply(lambda x : x.replace(".", "_")) == element]['Color'].values
    if a.size > 0 :
        color.append(a[0])
    else :
        color.append("#8c8c8c")
data2['color'] = color


st.subheader("Clustering parameters")
col1, col2 = st.columns(2)

with col1 :
    g = sns.barplot(x = 'mgCST', y = 'count_sample', data = data2 , legend = True, hue = 'mgCST', palette=list(data2['color']))
    fig1 = g.figure
    plt.xlabel("mgCSTs")
    plt.ylabel("Number of samples")
    g.tick_params(axis='x', which='major', labelsize= 8, labelrotation=70)
    # g.legend(title = "Dominant Taxa", loc='upper center', bbox_to_anchor=(0.5, -0.2), fancybox=True, shadow=True, ncol = 2)
    new_legend = data2['color'].apply(lambda x : mcolors.to_rgba(x)).values

    new_patch = []
    for i,j in zip(new_legend, data2['domTaxa'].values) :
        new_patch.append(mpatches.Patch(color = i, label = j))

    g.legend(handles=new_patch, title = 'mgCSTs',loc='upper center', bbox_to_anchor=(0.5, -0.2), fancybox=True, shadow=True, ncol = 2)


    g.grid(False)
    st.pyplot(fig1)



# Project Colors
proj_cols = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#d1ba36', '#a65628', '#f781bf', '#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e']
df2 = data.groupby(["Project", "mgCST"]).size().reset_index().pivot(columns='Project', index = 'mgCST', values =0)

with col2 :
    h = df2.plot(kind='bar', stacked=True, color=proj_cols)
    fig2 = h.figure
    plt.xlabel("mgCSTs")
    plt.ylabel("Number of samples")
    h.legend(title = "Project",loc='upper right', fontsize="small")
    h.tick_params(axis='x', which='major', labelsize= 8, labelrotation=70)
    h.legend(title = "Project",loc='upper center', bbox_to_anchor=(0.5, -0.2), fancybox=True, shadow=True, ncol = 3)
    h.grid(False)
    st.pyplot(fig2)

st.subheader("Most abund species per mgCSTs")
with st.expander("See tables"):
    st.dataframe(data2[['mgCST','domTaxa','meanRelabund','color']])


# st.subheader("Colors signification")
# with st.expander("See colors"):
#     col1, col2 = st.columns(2)
#     a = int(len(data2['domTaxa'].unique())/2)

#     with col1 :
#         i = 0
#         for taxa in (data2['domTaxa'].unique()[:a+1]) :
#             color_taxa = data2['color'][data2['domTaxa'] == taxa].values[0]
#             mgCST = data2[data2['domTaxa'] == taxa]['mgCST'].values[0]
#             color = st.color_picker(taxa, color_taxa, key=i)
#             i+=1

#     with col2 :
#         i = -1
#         for taxa in data2['domTaxa'].unique()[-a:] :
#             color_taxa = data2['color'][data2['domTaxa'] == taxa].values[0]
#             color = st.color_picker(taxa, color_taxa, key=i*2)
#             i-=1
