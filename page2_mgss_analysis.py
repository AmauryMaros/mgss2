import streamlit as st
import pandas as pd
from PIL import Image

st.set_page_config(layout="centered")

clust = pd.read_csv("Data/mgss.clustering.parameters.csv").rename(columns={"Unnamed: 0" : "species"})
mgss_coverage = pd.read_csv("Data/mgSs.coverage.stats.csv").rename(columns={"Unnamed: 0" : "sub_species", "as.factor(sample_cluster)" : "sample_cluster"})
species = mgss_coverage['sub_species'].apply(lambda x : x.split(".")[0]).unique()

st.title("Metagenomic Subspecies")

option = st.selectbox('Species', species)

tab1, tab2, tab3 = st.tabs(["Species coverage", "Subspecies stats", "Presence Absence"])


with tab1:
        
        col1, col2 = st.columns(2)
        with col1 :
            fig1 = Image.open("Medias/mgss_coverage_png/" + option + "_subspecies_coverage_boxplot.png") 
            st.image(fig1)

        with col2 :
            fig2 = Image.open("Medias/mgss_coverage_png/" + option + "_subspecies_coverage_by_NoGenes.png")
            st.image(fig2)

with tab2 :

        st.subheader("Clustering parameters")
        df2 = clust[clust['species'] == option].drop('Number_of_unassigned', axis = 1)
        st.dataframe(df2)
        
        st.subheader("Subspecies stats")
        st.dataframe(mgss_coverage[mgss_coverage['sub_species'].apply(lambda x : x.split(".")[0]) == option])

with tab3 :
        image = Image.open("Medias/heatmap_presence_absence/_" + option + "_heatmap_presence_absence.png")
        st.image(image)
