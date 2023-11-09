import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import pickle
from dynamicTreeCut import cutreeHybrid
import seaborn as sns

st.set_page_config(layout="centered")

with open("Data/samples_dist.pkl", "rb") as file:
    samples_dist = pickle.load(file)

with open("Data/samples_hc.pkl", "rb") as file:
    samples_hc = pickle.load(file)

mgss_coverage = pd.read_csv("Data/mgSs.coverage.stats.csv").rename(columns={"Unnamed: 0" : "sub_species", "as.factor(sample_cluster)" : "sample_cluster"})
species = mgss_coverage['sub_species'].apply(lambda x : x.split(".")[0]).unique()

st.sidebar.subheader("Parameters")
deepsplit = st.sidebar.slider(label="deepSplit", min_value=0, max_value=4)
mincluster = st.sidebar.slider(label="minClusterSize", min_value=10, max_value=50)

st.subheader("Clustering parameters")
option = st.selectbox("Species", species)

cut_tree = cutreeHybrid(samples_hc[option], samples_dist[option], minClusterSize=mincluster, deepSplit=deepsplit)

data = pd.DataFrame({"labels":cut_tree['labels']})

fig,ax = plt.subplots()
ax = sns.countplot(x = 'labels', data = data, hue = 'labels', legend=False)
plt.xlabel("Clusters")
plt.xticks(data['labels'].value_counts().index)
plt.ylabel("Number of samples")
ax.grid(False)
st.pyplot(fig)
