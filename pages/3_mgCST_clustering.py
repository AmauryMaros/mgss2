import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
from dynamicTreeCut import cutreeHybrid
import seaborn as sns
import numpy as np

sns.set_theme(style="whitegrid", palette="cubehelix")

df = pd.read_csv("Data/mgCSTs.params.df.csv").rename(columns={'dtc' : 'mgCST'})
projects = pd.read_csv("Data/VIRGO2_projects.csv")
species_abund = pd.read_csv("Data/mgCSTs.df.csv").rename(columns={'dtc' : 'mgCST'})
colors_mgCSTs = pd.read_csv("Data/mgCST_color.csv")

st.sidebar.subheader("Parameters")
deepsplit = st.sidebar.slider(label="deepSplit", min_value=0, max_value=4)
mincluster = st.sidebar.slider(label="minClusterSize", min_value=10, max_value=50)

data = df[(df['deepSplit'] == deepsplit) & (df['minClusterSize'] == mincluster)]
data = data.reset_index(drop = True)
data = data.join(projects)
df2 = data.groupby(["Project", "mgCST"]).size().reset_index().pivot(columns='Project', index = 'mgCST', values =0)

most_abundant = species_abund[(species_abund['deepSplit'] == deepsplit) & (species_abund['minClusterSize'] == mincluster)]
most_abundant = most_abundant[['mgCST', 'domTaxa', 'meanRelabund']]


couleur = {k:v for k,v in zip(most_abundant['domTaxa'].unique(), colors_mgCSTs['color'][:len(most_abundant['domTaxa'].unique())])}
liste = []
for i in most_abundant['domTaxa'] :
    liste.append(couleur[i])

most_abundant['color'] = liste

st.subheader("Clustering parameters")

col1, col2 = st.columns(2)

with col1 :
    g = sns.countplot(x = 'mgCST', data = data, hue='mgCST', palette=most_abundant['color'].values, legend = False)
    fig1 = g.figure
    plt.xlabel("mgCSTs")
    plt.ylabel("Number of samples")
    g.tick_params(axis='x', which='major', labelsize= 8, labelrotation=70)
    # g.legend(title = 'mgCSTs',loc='upper center', bbox_to_anchor=(0.5, -0.2), fancybox=True, shadow=True, ncol = 5)
    g.grid(False)
    st.pyplot(fig1)

# Project Colors
proj_cols = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#d1ba36', '#a65628', '#f781bf', '#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e']

with col2 :
    h = df2.plot(kind='bar', stacked=True, color=proj_cols)
    fig2 = h.figure
    plt.xlabel("mgCSTs")
    plt.ylabel("Number of samples")
    h.legend(title = "Project",loc='upper right', fontsize="small")
    h.tick_params(axis='x', which='major', labelsize= 8, labelrotation=70)
    # h.legend(title = "Project",loc='upper center', bbox_to_anchor=(0.5, -0.2), fancybox=True, shadow=True, ncol = 3)
    h.grid(False)
    st.pyplot(fig2)

st.subheader("Most abund species per mgCSTs")
with st.expander("See tables"):
    col1,col2 = st.columns(2)
    with col1 :
        st.dataframe(most_abundant.set_index('mgCST'))
    with col2 :
        st.dataframe(most_abundant['domTaxa'].value_counts())

st.subheader("Colors signification")
col1, col2 = st.columns(2)
a = int(len(most_abundant['domTaxa'].unique())/2)

with col1 :
    i = 0
    for taxa in (most_abundant['domTaxa'].unique()[:a+1]) :
        color_taxa = most_abundant['color'][most_abundant['domTaxa'] == taxa].values[0]
        color = st.color_picker(taxa, color_taxa, key=i)
        i+=1

with col2 :
    i = -1
    for taxa in most_abundant['domTaxa'].unique()[-a:] :
        color_taxa = most_abundant['color'][most_abundant['domTaxa'] == taxa].values[0]
        color = st.color_picker(taxa, color_taxa, key=i*2)
        i-=1
