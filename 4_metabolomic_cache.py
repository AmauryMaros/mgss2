import streamlit as st
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
from sklearn.decomposition import PCA


@st.cache_data
def read_csv(url):
     return pd.read_csv(url)
color_mgCST = read_csv("Data/mgCST_sort_color.csv")
mgcsts_samples_df = read_csv("Data/mgCSTs.samples.df.csv")
project = read_csv("Data/VIRGO2_projects.csv")

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

# Assign mgCST color for each sample
color_mgCST = color_mgCST[['mgCST', 'color_mgCST']].reset_index(drop = True)

st.subheader("PCA analysis - all samples")
# Plot the figure
col1, col2 = st.columns(2)
pca_df = pd.DataFrame(data=principal_components, columns=['PC1', 'PC2'])
new_legend = color_mgCST['color_mgCST'].apply(lambda x : mcolors.to_rgba(x)).values
proj = pd.merge(sampleID, project, on='sampleID', how='inner')['Project']

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
st.dataframe(pd.DataFrame(mgCST.value_counts().sort_index().rename(index = "Number of samples")).transpose())
st.write('Minclustersize :', minclustersize)
st.write('Deepsplit : ', deepsplit)


# # Comparison of 2 groups
st.subheader("Comparison between 2 groups")

st.container()
# PCA comparison between 2 different mgCST - same figure
col1, col2 = st.columns(2)
with col1 :
    mgCST1 = st.slider("GroupA", 1, color_mgCST.shape[0],(1,color_mgCST.shape[0]), key='mgCST1')
with col2 :
    mgCST2 = st.slider("GroupB", 1,color_mgCST.shape[0],(1,color_mgCST.shape[0]), key='mgCST2')

# st.container()
@st.cache_data
def data_pca(m1,m2):
    df1 = mgcsts_samples[mgcsts_samples['mgCST'].isin(range(m1[0],m1[1]+1))]
    df1.loc[:,'label'] = "GroupA"
    df2 = mgcsts_samples[mgcsts_samples['mgCST'].isin(range(m2[0],m2[1]+1))]
    df2.loc[:,'label'] = "GroupB"
    df = pd.concat([df1,df2], axis = 0)
    data1 = pd.merge(df, metabolomics, on = "sampleID", how = "inner")
    mgCSTs = data1['mgCST']
    label = data1['mgCST']
    data1 = data1.drop(['dtc','domTaxa','relabund','minClusterSize','deepSplit', 'sampleID', 'mgCST', 'label'], axis = 1)
    return data1, mgCSTs, label

data1, mgCSTs, groups = data_pca(mgCST1, mgCST2)

@st.cache_data
def pca_group(data):
    pca = PCA(n_components=3)
    df = pd.DataFrame(data = pca.fit_transform(data), columns = ['PC1','PC2','PC3'])
    explained_variance = pca.explained_variance_ratio_
    return df, explained_variance
pca = pca_group(data_pca(mgCST1, mgCST2)[0])

new_colors = color_mgCST[color_mgCST['mgCST'].isin(mgCSTs)]
new_legend = new_colors['color_mgCST'].apply(lambda x : mcolors.to_rgba(x)).values

tab1, tab2 = st.tabs(['PC1 PC2', 'PC2 PC3'])


      
@st.cache_data
def display_pca(df, pc0, pc1, explained_variance, a ,b):
        fig,f = plt.subplots()
        f = sns.scatterplot(x=df[pc0], y=df[pc1], hue=mgCSTs, palette=list(new_legend), style = groups)
        plt.xlabel(pc0 + ' : ' + str(round(explained_variance[a]*100,2)) + "%")
        plt.ylabel(pc1 + ' : ' + str(round(explained_variance[b]*100,2)) + "%")
        plt.legend(loc='best')
        # new_patch = []
        # for i,j in zip(new_legend, new_colors['mgCST'].values) :
        #     new_patch.append(mpatches.Patch(color = i, label = j))
        # f.legend(handles=new_patch, title = 'mgCSTs',loc='center',bbox_to_anchor=(1.2, 0.5), fontsize = "7",fancybox=True, shadow=True, ncol = 2)
        return fig

with tab1 :
    run_pca12 = st.button('Run PC1 PC2', key='pc12')
    if run_pca12 :
        st.pyplot(display_pca(pca[0], 'PC1', 'PC2',pca[1], 0, 1))

with tab2 :
    run_pca23 = st.button('Run PC2 PC3', key='pc23')
    if run_pca23:        
        st.pyplot(display_pca(pca[0], 'PC2', 'PC3', pca[1], 1, 2))


