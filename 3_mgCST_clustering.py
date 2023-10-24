import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
import subprocess
import base64  

# sns.set_theme(style="whitegrid", palette="cubehelix")

subspecies_with_colors = pd.read_csv("Data/subspecies_with_colors.csv")
mgcsts_samples = pd.read_csv("Data/mgCST_samples_color.csv")
mgCSTs_sort = pd.read_csv("Data/mgCST_sort_color.csv")
projects = pd.read_csv("Data/VIRGO2_projects.csv")


st.sidebar.subheader("Parameters")
deepsplit = st.sidebar.slider(label="deepSplit", min_value=0, max_value=4)
mincluster = st.sidebar.slider(label="minClusterSize", min_value=10, max_value=50)

# Save parameters into a csv file for the execution of Rscript to built the heatmap
parameters = pd.DataFrame({"minClusterSize" : [mincluster], "deepsplit" : [deepsplit]})
parameters.to_csv("R_scripts/mgCSTs_parameters_streamlit.csv", index = False)

data = mgcsts_samples[(mgcsts_samples['deepSplit'] == deepsplit) & (mgcsts_samples['minClusterSize'] == mincluster)]
data = data.reset_index(drop = True)
data = data.merge(projects, on = "sampleID", how = "left")

data2 = mgCSTs_sort[(mgCSTs_sort['deepSplit'] == deepsplit) & (mgCSTs_sort['minClusterSize'] == mincluster)]
data2 = data2.reset_index(drop = True)

count_sample = []
for element in data2['dtc'].values :
    count_sample.append(data.groupby(['dtc']).count()['sampleID'][element])
data2['count_sample'] = count_sample

st.subheader("Clustering parameters")
col1, col2 = st.columns(2)

with col1 :
    g = sns.barplot(x = 'mgCST', y = 'count_sample', data = data2 , legend = True, hue = 'mgCST', palette=list(data2['color_mgCST']))
    fig1 = g.figure
    plt.xlabel("mgCSTs")
    plt.ylabel("Number of samples")
    g.tick_params(axis='x', which='major', labelsize= 8, labelrotation=70)
    # g.legend(title = "Dominant Taxa", loc='upper center', bbox_to_anchor=(0.5, -0.2), fancybox=True, shadow=True, ncol = 2)
    new_legend = data2['color_mgCST'].apply(lambda x : mcolors.to_rgba(x)).values

    new_patch = []
    for i,j in zip(new_legend, data2['domTaxa'].values) :
        new_patch.append(mpatches.Patch(color = i, label = j))

    g.legend(handles=new_patch, title = 'Dominant taxa',loc='upper center', bbox_to_anchor=(0.5, -0.2), fancybox=True, shadow=True, ncol = 2)
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
with st.expander("See table"):
    st.dataframe(data2[['mgCST','domTaxa','meanRelabund','color_mgCST']])

heatmap = st.button(label = "Generate heatmap (1min30)")
if heatmap :
    st.subheader("MgCSTs Heatmap")

    process = subprocess.Popen(["Rscript", "mgCSTs_heatmap.R"])
    result = process.communicate()

    def displayPDF(file):
        # Opening file from file path
        with open(file, "rb") as f:
            base64_pdf = base64.b64encode(f.read()).decode('utf-8')

        # Embedding PDF in HTML
        pdf_display = F'<iframe src="data:application/pdf;base64,{base64_pdf}" width="700" height="1000" type="application/pdf"></iframe>'

        # Displaying File
        st.markdown(pdf_display, unsafe_allow_html=True)

    displayPDF("Medias/mgCST_heatmap.pdf")