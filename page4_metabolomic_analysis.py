import streamlit as st
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
from sklearn.decomposition import PCA
from page3_mgCST_clustering import data2

@st.cache_data
def read_csv(url):
     return pd.read_csv(url)
color_mgCST = read_csv("Data/mgCST_sort_color.csv")
mgcsts_samples_df = read_csv("Data/mgCSTs.samples.df.csv")
project = read_csv("Data/VIRGO2_projects.csv")
mgCSTs_sort = read_csv("Data/mgCST_sort_color.csv")
data = read_csv("Data/mgCST_samples_color.csv")

# Slider for parameters variation
st.sidebar.subheader("Parameters")
parameters = pd.read_csv("R_scripts/mgCSTs_parameters_streamlit.csv")
value_ds = parameters['deepsplit'].values[0]
value_mcs = parameters['minClusterSize'].values[0]
deepsplit = st.sidebar.slider(label="deepSplit", min_value=0, max_value=4, value = value_ds)
minclustersize = st.sidebar.slider(label="minClusterSize", min_value=10, max_value=50, value = value_mcs)
parameters = pd.DataFrame({"minClusterSize" : [minclustersize], "deepsplit" : [deepsplit]})
parameters.to_csv("R_scripts/mgCSTs_parameters_streamlit.csv", index = False)

@st.cache_data
def read_pickle(path):
    with open(path, 'rb') as file:
        df = pickle.load(file)
    return df
metabolomics = read_pickle("/Users/amaros/Desktop/mgss2/log_norm.pkl")
pca_model =  read_pickle("Data/pca_model.pkl")

@st.cache_data
def pca_model_2_compo_load(minclustersize, deepsplit):
    return pca_model[minclustersize][deepsplit][0]
principal_components = pca_model_2_compo_load(minclustersize, deepsplit)['principal_components']
explained_variance = pca_model_2_compo_load(minclustersize, deepsplit)['explained_var_ration']
sampleID  = pca_model_2_compo_load(minclustersize, deepsplit)['sampleID']
mgCST = pca_model_2_compo_load(minclustersize, deepsplit)['mgCST']

@st.cache_data
def filter_df(df, minclustersize, deepsplit):
     return df[(df['minClusterSize'] == minclustersize) & (df['deepSplit'] == deepsplit)]
mgcsts_samples = filter_df(mgcsts_samples_df, minclustersize, deepsplit)
color_mgCST = filter_df(color_mgCST, minclustersize, deepsplit)
data = mgcsts_samples[(mgcsts_samples['deepSplit'] == deepsplit) & (mgcsts_samples['minClusterSize'] == minclustersize)]
data2 = mgCSTs_sort[(mgCSTs_sort['deepSplit'] == deepsplit) & (mgCSTs_sort['minClusterSize'] == minclustersize)]

data2 = data2.reset_index(drop = True)
count_sample = []
for element in data2['dtc'].values :
    count_sample.append(data.groupby(['dtc']).count()['sampleID'][element])
data2['count_sample'] = count_sample

# Assign mgCST color for each sample
color_mgCST = color_mgCST[['mgCST', 'color_mgCST']].reset_index(drop = True)
color_mgCST = color_mgCST[color_mgCST['mgCST'].isin(mgCST)]

st.subheader("PCA analysis - all samples")
# Plot the figure
col1, col2 = st.columns(2)
pca_df = pd.DataFrame(data=principal_components, columns=['PC1', 'PC2'])
new_legend = color_mgCST['color_mgCST'].apply(lambda x : mcolors.to_rgba(x)).values
proj = pd.merge(sampleID, project, on='sampleID', how='inner')['Project']
st.container()
with col1 :
    fig,f = plt.subplots()
    f = sns.scatterplot(pca_df, x='PC1', y='PC2', hue=mgCST, palette=new_legend)
    plt.xlabel('PC1 : ' + str(round(explained_variance[0]*100,2)) + "%")
    plt.ylabel('PC2 : ' + str(round(explained_variance[1]*100,2)) + "%")
    new_legend = color_mgCST['color_mgCST'].apply(lambda x : mcolors.to_rgba(x)).values
    new_patch = []
    for i,j in zip(new_legend, color_mgCST['mgCST'].values) :
        new_patch.append(mpatches.Patch(color = i, label = j))
    f.legend(handles=new_patch, title = 'mgCSTs',loc='center',bbox_to_anchor=(1.2, 0.5), fontsize = "7",fancybox=True, shadow=True, ncol = 2)
    st.pyplot(fig)

with  col2 :   
    fig, g = plt.subplots()
    g = sns.scatterplot(data = pca_df, x='PC1', y='PC2', hue=proj)
    plt.xlabel('PC1 : ' + str(round(explained_variance[0]*100,2)) + "%")
    plt.ylabel('PC2 : ' + str(round(explained_variance[1]*100,2)) + "%")
    g.legend(title = 'Projects',loc='center',bbox_to_anchor=(1.2, 0.5), fontsize = "7",fancybox=True, shadow=True, ncol = 1)
    st.pyplot(fig)

st.container()
st.write("Metabolome Data : 1280 samples")
st.dataframe(pd.DataFrame(mgCST.value_counts().sort_index().rename(index = "count_sample")).transpose())
st.write("All samples : 2528")
st.dataframe(data2[['mgCST','count_sample']].set_index('mgCST').transpose())
st.write('Minclustersize :', minclustersize)
st.write('Deepsplit : ', deepsplit)




# # Comparison of 2 groups
st.subheader("Comparison between 2 groups")

st.container()
# PCA comparison between 2 different mgCST - same figure
col1, col2 = st.columns(2)
with col1 :
    # mgCST1 = st.slider("GroupA", 1, color_mgCST.shape[0],(1,color_mgCST.shape[0]), key='mgCST1')
    mgCST1 = st.slider("GroupA", 1, 41,(1,41), key='mgCST1')
with col2 :
    # mgCST2 = st.slider("GroupB", 1,color_mgCST.shape[0],(1,color_mgCST.shape[0]), key='mgCST2')
    mgCST2 = st.slider("GroupB", 1, 41,(1,41), key='mgCST2')

# st.container()
@st.cache_data
def data_pca(m1,m2):
    df1 = mgcsts_samples[mgcsts_samples['mgCST'].isin(range(m1[0],m1[1]+1))]
    df1.loc[:,'label'] = "GroupA"
    df2 = mgcsts_samples[mgcsts_samples['mgCST'].isin(range(m2[0],m2[1]+1))]
    df2.loc[:,'label'] = "GroupB"
    df = pd.concat([df1,df2], axis = 0)
    data1 = pd.merge(df, metabolomics, on = "sampleID", how = "inner")
    return data1

data1 = data_pca(mgCST1, mgCST2)


mgCSTs = data1['mgCST']
data1 = data1.drop(['dtc','domTaxa','relabund','minClusterSize','deepSplit', 'sampleID', 'mgCST', 'label'], axis = 1)
new_colors = color_mgCST[color_mgCST['mgCST'].isin(mgCSTs)]


run_pca = st.button('Run PCA', key='run_pca')

if run_pca :
    # Create and fit your PCA models
    pca = PCA(n_components=3)
    principal_components = pca.fit_transform(data1)
    explained_variance = pca.explained_variance_ratio_

    pca_df = pd.DataFrame(data = principal_components, columns=['PC1', 'PC2', 'PC3'])
    new_legend = new_colors['color_mgCST'].apply(lambda x : mcolors.to_rgba(x)).values

    @st.cache_data
    def display_pca(df, pc0, pc1, a ,b):
        fig,f = plt.subplots()
        f = sns.scatterplot(x=df[pc0], y=df[pc1], hue=mgCSTs, palette=list(new_legend))
        plt.xlabel(pc0 + ' : ' + str(round(explained_variance[a]*100,2)) + "%")
        plt.ylabel(pc1 + ' : ' + str(round(explained_variance[b]*100,2)) + "%")
        new_patch = []
        for i,j in zip(new_legend, new_colors['mgCST'].values) :
            new_patch.append(mpatches.Patch(color = i, label = j))
        f.legend(handles=new_patch, title = 'mgCSTs',loc='center',bbox_to_anchor=(1.2, 0.5), fontsize = "7",fancybox=True, shadow=True, ncol = 2)
        return fig
    
    tab1, tab2 = st.tabs(['PC1 PC2', 'PC2 PC3'])
    with tab1 :    
        st.pyplot(display_pca(pca_df, 'PC1', 'PC2', 0, 1))
    with tab2 :
        st.pyplot(display_pca(pca_df, 'PC2', 'PC3', 1, 2))


