import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import subprocess
import base64  
import plotly.express as px
import numpy as np

st.set_page_config(layout="wide")

# subspecies_with_colors = pd.read_csv("Data/subspecies_with_colors.csv")
mgcsts_samples_df = pd.read_csv("Data/mgCSTs.samples.df2.csv")
mgCSTs_sort = pd.read_csv("Data/mgCST_sort_color2.csv")
projects = pd.read_csv("Data/VIRGO2_projects.csv")


# Slider for parameters variation
st.sidebar.subheader("Parameters")
parameters = pd.read_csv("Data/mgCSTs_parameters_streamlit.csv")
value_ds = parameters['deepsplit'].values[0]
value_mcs = parameters['minClusterSize'].values[0]
deepsplit = st.sidebar.slider(label="deepSplit", min_value=0, max_value=4, value = value_ds)
minclustersize = st.sidebar.slider(label="minClusterSize", min_value=10, max_value=50, value = value_mcs)
parameters = pd.DataFrame({"minClusterSize" : [minclustersize], "deepsplit" : [deepsplit]})
parameters.to_csv("Data/mgCSTs_parameters_streamlit.csv", index = False)

mgcsts_samples = mgcsts_samples_df[(mgcsts_samples_df['deepSplit'] == deepsplit) & (mgcsts_samples_df['minClusterSize'] == minclustersize)]
mgcsts_samples = mgcsts_samples.reset_index(drop = True)
mgcsts_samples = mgcsts_samples.merge(projects, on = "sampleID", how = "left")

mgcsts = mgCSTs_sort[(mgCSTs_sort['deepSplit'] == deepsplit) & (mgCSTs_sort['minClusterSize'] == minclustersize)]
mgcsts = mgcsts.reset_index(drop = True)

count_sample = []
for element in mgcsts['dtc'].values :
    count_sample.append(mgcsts_samples.groupby(['dtc']).count()['sampleID'][element])
mgcsts['count_sample'] = count_sample

st.container()
st.subheader("Clustering parameters")

mgcsts['mgCST'] = mgcsts['mgCST'].astype(str)

col1, col2, col3 = st.columns(3)
with col1 :
    # fig1,g = plt.subplots()
    # g = sns.barplot(x = 'mgCST', y = 'count_sample', data = mgcsts , legend = True, hue = 'mgCST', palette=list(mgcsts['color_mgCST']),edgecolor='black', linewidth=0.5)
    # plt.xlabel("mgCST")
    # plt.ylabel("Number of samples")
    # plt.tick_params(axis='x', which='major', labelsize= 8, labelrotation=70)
    # g.legend(title = "mgCSTs", loc='upper center', bbox_to_anchor=(0.5, -0.2), fancybox=True, shadow=True, ncol = 5)
    # g.grid(False)
    # st.pyplot(fig1)
    fig = px.bar(mgcsts, x='mgCST', y='count_sample', color='mgCST',
                 hover_data=['domTaxa'],
                 color_discrete_sequence=list(mgcsts['color_mgCST']), 
                 labels={"count_sample" : "Number of samples"},
                 title = "Number of samples on each MgCSTs - colored by mgCST")
    subtitle_text = f"MinClusterSize = {minclustersize} <br> DeepSplit = {deepsplit}"
    fig.add_annotation(
        dict(text=subtitle_text,
            xref='paper', yref='paper',
            x=0, y=1.05,
            showarrow=False,
            font=dict(size=12)
            ))

    fig.update_xaxes(tickmode='array', tickvals=mgcsts['mgCST'])
    fig.update_layout(
    xaxis=dict(
        tickangle=0,  # Change the angle of the tick labels
        tickfont=dict(size=10)  # Change the size of the tick labels
        ))
    st.plotly_chart(fig, use_container_width=True)

with col2 :
    # Project Colors
    proj_cols = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#d1ba36', '#a65628', '#f781bf', '#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e']
    df2 = mgcsts_samples.groupby(["Project", "mgCST"]).size().reset_index().pivot(columns='Project', index = 'mgCST', values =0).reset_index()
    # df2 = df2.reset_index()
    # st.dataframe(df2)
    # st.write(df2.columns)
    # mgcsts_samples
#     # Plot the figure
#     h = df2.plot(kind='bar', stacked=True, color=proj_cols, edgecolor='black', linewidth=0.5)
#     fig2 = h.figure
#     plt.xlabel("mgCST")
#     plt.ylabel("Number of samples")
#     h.legend(title = "Project",loc='upper right', fontsize="small")
#     h.tick_params(axis='x', which='major', labelsize= 8, labelrotation=70)
#     h.legend(title = "Project",loc='upper center', bbox_to_anchor=(0.5, -0.2), fancybox=True, shadow=True, ncol = 3)
#     h.grid(False)
#     st.pyplot(fig2)
    fig = px.bar(df2, x='mgCST', y=df2.drop('mgCST', axis=1).columns.values, labels={'value' : 'Number of samples'},
                 title = f"Number of samples on each MgCSTs - colored by projects")
    subtitle_text = f"MinClusterSize = {minclustersize} <br> DeepSplit = {deepsplit}"
    fig.add_annotation(
        dict(text=subtitle_text,
            xref='paper', yref='paper',
            x=0, y=1.05,
            showarrow=False,
            font=dict(size=12)
            )
    )
    fig.update_xaxes(tickmode='array', tickvals=mgcsts['mgCST'])
    fig.update_layout(
    xaxis=dict(
        tickangle=0,  # Change the angle of the tick labels
        tickfont=dict(size=10)  # Change the size of the tick labels
        ))
    st.plotly_chart(fig, use_container_width=True)

with col3 :
    File_S6 = pd.read_excel('Data/File_S6_clean.xlsx')
    old_mgCST = File_S6[['mapID', 'mgCST']]
    old_mgCST = old_mgCST.rename(columns={'mapID':'sampleID'})
    new_mgCST = mgcsts_samples[(mgcsts_samples['deepSplit'] == deepsplit) & (mgcsts_samples['minClusterSize'] == minclustersize)][['sampleID', 'mgCST']]
    new_vs_old = pd.merge(new_mgCST, old_mgCST, on='sampleID', how='inner')
    bubble_data = new_vs_old.groupby(['mgCST_x', 'mgCST_y']).size().reset_index(name='count')
    bubble_data = bubble_data.rename(columns={"mgCST_x":"mgCST", "mgCST_y":"old_mgCST"})

    bubble_color = []
    mgcsts['mgCST'] = mgcsts['mgCST'].astype(int)
    for i in sorted(bubble_data['mgCST'].unique()):
        bubble_color.append(mgcsts[mgcsts['mgCST'] == i]['color_mgCST'].values[0])
        
    # # Plot the figure with seaborn
    # fig3, f = plt.subplots()
    # f = sns.scatterplot(data=bubble_data, x = 'mgCST', y = 'old_mgCST', hue='mgCST', size='count', edgecolor='black', palette=list(bubble_color))
    # plt.xlabel("mgCST")
    # plt.ylabel("old_mgCST")
    # plt.grid(True)
    # plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), fancybox=True, shadow=True, ncol = 5)
    # st.pyplot(fig3)
    

    bubble_data['mgCST'] = bubble_data['mgCST'].astype(str)
    fig = px.scatter(bubble_data, x='mgCST', y = 'old_mgCST', color='mgCST',color_discrete_sequence=list(bubble_color), size='count',
                     title = "Number of samples in previous MgCST vs New MgCST")
    subtitle_text = f"MinClusterSize = {minclustersize} <br> DeepSplit = {deepsplit}"
    fig.add_annotation(
        dict(text=subtitle_text,
            xref='paper', yref='paper',
            x=0, y=1.05,
            showarrow=False,
            font=dict(size=12)
            )
    )
    fig.update_xaxes(tickmode='array', tickvals=mgcsts['mgCST'])
    fig.update_yaxes(tickmode='array', tickvals=np.arange(1,28))
    fig.update_layout(
    xaxis=dict(
        tickangle=0,  # Change the angle of the tick labels
        tickfont=dict(size=10)  # Change the size of the tick labels
        ),
    yaxis=dict(
        tickangle=0,  # Change the angle of the tick labels
        tickfont=dict(size=10)  # Change the size of the tick labels
        ),
    )
    st.plotly_chart(fig, use_container_width=True)


st.subheader("Most abund species per mgCSTs")
with st.expander("See table"):
    st.dataframe(mgcsts[['mgCST','domTaxa','meanRelabund','color_mgCST']])

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