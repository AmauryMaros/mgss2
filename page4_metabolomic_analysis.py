import streamlit as st
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import matplotlib.colors as mcolors
from sklearn.decomposition import PCA
import numpy as np
import plotly.express as px

st.set_page_config(layout="wide")


# Slider for parameters variation
st.sidebar.subheader("Parameters")
parameters = pd.read_csv("Data/mgCSTs_parameters_streamlit.csv")
value_ds = parameters['deepsplit'].values[0]
value_mcs = parameters['minClusterSize'].values[0]
deepsplit = st.sidebar.slider(label="deepSplit", min_value=0, max_value=4, value = value_ds)
minclustersize = st.sidebar.slider(label="minClusterSize", min_value=10, max_value=50, value = value_mcs)
parameters = pd.DataFrame({"minClusterSize" : [minclustersize], "deepsplit" : [deepsplit]})
parameters.to_csv("Data/mgCSTs_parameters_streamlit.csv", index = False)

# Load Data
@st.cache_data
def read_csv(url):
     return pd.read_csv(url)
mgcsts_samples_df = read_csv("Data/mgCSTs.samples.df.csv")
project = read_csv("Data/VIRGO2_projects.csv")
mgCSTs_sort = read_csv("Data/mgCST_sort_color.csv")

@st.cache_data
def read_pickle(path):
    with open(path, 'rb') as file:
        df = pickle.load(file)
    return df
metabolomics = read_pickle("/Users/amaros/Desktop/mgss2/log_norm.pkl")      # change this path
pca_model =  read_pickle("Data/pca_model.pkl")

@st.cache_data
def pca_model_data(minclustersize, deepsplit):
    return pca_model[minclustersize][deepsplit][0]
principal_components = pca_model_data(minclustersize, deepsplit)['principal_components']
explained_variance = pca_model_data(minclustersize, deepsplit)['explained_var_ratio']
sampleID  = pca_model_data(minclustersize, deepsplit)['sampleID']
mgCST = pca_model_data(minclustersize, deepsplit)['mgCST']

@st.cache_data
def filter_df(df, minclustersize, deepsplit):
     return df[(df['minClusterSize'] == minclustersize) & (df['deepSplit'] == deepsplit)]
mgcsts_samples = filter_df(mgcsts_samples_df, minclustersize, deepsplit)
mgcsts = mgCSTs_sort[(mgCSTs_sort['deepSplit'] == deepsplit) & (mgCSTs_sort['minClusterSize'] == minclustersize)]

mgcsts = mgcsts.reset_index(drop = True)
count_sample = []
for element in mgcsts['dtc'].values :
    count_sample.append(mgcsts_samples.groupby(['dtc']).count()['sampleID'][element])
mgcsts['count_sample'] = count_sample

# Assign mgCST color for each sample
color_mgCST = mgcsts[['mgCST', 'color_mgCST']].reset_index(drop = True)
color_mgCST = color_mgCST[color_mgCST['mgCST'].isin(mgCST)]

st.header("PCA - all samples", divider='gray')

# Create a dataframe with PCA datas for visualization
principal_components = pd.DataFrame(data=principal_components, columns = ['PC1','PC2','PC3'])
pca_df = pd.concat([sampleID, principal_components,mgCST], axis=1)
pca_df = pd.merge(pca_df, project, on='sampleID', how='inner') 
pca_df = pd.merge(pca_df, color_mgCST, on='mgCST', how='inner')      
pca_df = pca_df.sort_values(by='mgCST', ascending=True)
pca_df = pca_df.sort_values('mgCST', ascending=True).reset_index(drop=True)
pca_df['mgCST'] = pca_df['mgCST'].astype(str)   # to be interpret as categorical column for plotly

st.container()
col1, col2 = st.columns(2)
with col1 :
    # fig,f = plt.subplots()
    # f = sns.scatterplot(pca_df, x='PC1', y='PC2', hue='mgCST', palette=list(color_mgCST['color_mgCST'].values), edgecolor = 'black')
    # plt.xlabel('PC1 : ' + str(round(explained_variance[0]*100,2)) + "%")
    # plt.ylabel('PC2 : ' + str(round(explained_variance[1]*100,2)) + "%")
    # plt.legend(title = 'mgCSTs',loc='center',bbox_to_anchor=(1.2, 0.5), fontsize = "7",fancybox=True, shadow=True, ncol = 2)
    # st.pyplot(fig)
    # df = pca_df
    # # df = df.sort_values('mgCST', ascending=True)
    # df['mgCST'] = df['mgCST'].astype(str)

    fig = px.scatter(pca_df, x='PC1', y='PC2', color='mgCST',
                     color_discrete_sequence=color_mgCST['color_mgCST'].values,
                     title="PCA - colored by mgCSTs",
                     labels = {'PC1' : "PC1 : "+str(round(explained_variance[0]*100,2))+"%",
                               'PC2' : "PC2 : "+str(round(explained_variance[1]*100,2)) + "%"})
    st.plotly_chart(fig, use_container_width=True)



with  col2 :   
    # fig, g = plt.subplots()
    # g = sns.scatterplot(data = pca_df, x='PC1', y='PC2', hue='Project')
    # plt.xlabel('PC1 : ' + str(round(explained_variance[0]*100,2)) + "%")
    # plt.ylabel('PC2 : ' + str(round(explained_variance[1]*100,2)) + "%")
    # g.legend(title = 'Projects',loc='center',bbox_to_anchor=(1.2, 0.5), fontsize = "7",fancybox=True, shadow=True, ncol = 1)
    # st.pyplot(fig)
    fig = px.scatter(pca_df, x='PC1', y='PC2', color='Project',
                     title="PCA - colored by projects",
                     labels = {'PC1' : "PC1 : "+str(round(explained_variance[0]*100,2))+"%",
                               'PC2' : "PC2 : "+str(round(explained_variance[1]*100,2)) + "%"})
    st.plotly_chart(fig, use_container_width=True)

st.container()
col1, col2 = st.columns(2)
with col1 :
    st.write("Metabolome Data : 1280 samples")
    st.dataframe(pd.DataFrame(mgCST.value_counts().sort_index().rename(index = "count_sample")).transpose())
    st.write("All samples : 2528")
    st.dataframe(mgcsts[['mgCST','count_sample']].set_index('mgCST').transpose())
with col2 :
    st.write("Clustering parameters")
    st.write('Minclustersize :', minclustersize)
    st.write('Deepsplit : ', deepsplit)


# # Comparison of 2 groups
st.header("Comparison between 2 groups", divider='gray')

st.container()
# PCA comparison between 2 different mgCST - same figure
col1, col2 = st.columns(2)
with col1 :
    mgCST1 = st.slider("GroupA", 1, 41,(1,41), key='mgCST1')    # this return a tuple (a, b)
with col2 :
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
groups = data1['label']
data1 = data1.drop(['dtc','domTaxa','relabund','minClusterSize','deepSplit', 'sampleID', 'mgCST', 'label'], axis = 1)
new_colors = color_mgCST[color_mgCST['mgCST'].isin(mgCSTs)]


run_pca = st.button('Run PCA', key='run_pca')

import matplotlib.lines as mlines

if run_pca :
    # Create and fit your PCA models
    pca = PCA(n_components=6)
    principal_components = pca.fit_transform(data1)
    explained_variance = pca.explained_variance_ratio_
    components_compo = pca.components_

    pca_df = pd.DataFrame(data = principal_components, columns=['PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6'])
    pca_df = pd.concat([pca_df, mgCSTs], axis=1)
    pca_df = pca_df.sort_values(by='mgCST', ascending=True)
    new_legend = new_colors['color_mgCST']#.apply(lambda x : mcolors.to_rgba(x)).values
    pca_df['mgCST'] = pca_df['mgCST'].astype(str)
    
    @st.cache_data
    def display_pca(df, pc0, pc1, a ,b):
        # fig,f = plt.subplots()
        # f = sns.scatterplot(x=df[pc0], y=df[pc1], hue=mgCSTs, palette=list(new_legend), style=mgCSTs, edgecolor="black")
        # plt.xlabel(pc0 + ' : ' + str(round(explained_variance[a]*100,2)) + "%")
        # plt.ylabel(pc1 + ' : ' + str(round(explained_variance[b]*100,2)) + "%")
        # plt.legend(loc='center', bbox_to_anchor=(1.2, 0.5), fontsize = "7",fancybox=True, shadow=True, ncol = 1)
        # fig = px.scatter(df, x=df[pc0], y=df[pc1], color=mgCSTs)
        fig = px.scatter(df, x=pc0, y=pc1, color='mgCST', symbol='mgCST',
                         color_discrete_sequence=new_legend.values,
                         title=pc0+" - "+pc1+" representation",
                         labels = {pc0 : pc0 + " : " + str(round(explained_variance[a]*100,2))+"%",
                                   pc1 : pc1 + " : " + str(round(explained_variance[b]*100,2)) + "%"})
        return fig
    
    cols = metabolomics.columns
    
    @st.cache_data
    def top_features(pc, n) :
        arr = pca.components_[pc]
        n_largest = n
        max_indices = np.argpartition(-arr, n_largest)[:n_largest]
        min_indices = np.argpartition(arr, n_largest)[:n_largest]
        pc_df_pos = pd.DataFrame({'Features':[cols[i] for i in max_indices],
                              'explained_variance' : [arr[i] for i in max_indices]}).sort_values(by='explained_variance', ascending=False)
        pc_df_neg = pd.DataFrame({'Features':[cols[i] for i in min_indices],
                              'explained_variance' : [arr[i] for i in min_indices]}).sort_values(by='explained_variance', ascending=True)
        return pc_df_pos.reset_index(drop=True), pc_df_neg.reset_index(drop=True)
    
    st.subheader("PCA Visualization")
    # tab1, tab2, tab3 = st.tabs(['PC1 PC2', 'PC3 PC4', 'PC5 PC6'])
    col1, col2, col3 = st.columns(3)
    with col1 :
        # st.pyplot(display_pca(pca_df, 'PC1', 'PC2', 0, 1))
        st.plotly_chart(display_pca(pca_df, 'PC1', 'PC2', 0, 1), use_container_width=True)
    with col2 :
        # st.pyplot(display_pca(pca_df, 'PC3', 'PC4', 2, 3))
        st.plotly_chart(display_pca(pca_df, 'PC3', 'PC4', 2, 3), use_container_width=True)
    with col3 :
        # st.pyplot(display_pca(pca_df, 'PC5', 'PC6', 4, 5))
        st.plotly_chart(display_pca(pca_df, 'PC5', 'PC6', 4, 5), use_container_width=True)

    st.subheader("Features importances (explained variance)")
    # tab1, tab2, tab3, tab4, tab5, tab6 = st.tabs(['PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6'])
    # tab = [tab1, tab2, tab3, tab4, tab5, tab6]
    # for i in range(6):
    #     with tab[i] :
    #         col1, col2 = st.columns(2)
    #         with col1 :
    #             st.write("Top5 PC" + str(i+1) + " pos")
    #             st.dataframe(top_features(i,5)[0])
    #         with col2 :
    #             st.write("Top5 PC" + str(i+1) + " neg")
    #             st.dataframe(top_features(i,5)[1])

    # Number of features to consider
    n = 5
    PCs = ['PC1'] * n + ['PC2'] * n + ['PC3'] * n + ['PC4'] * n + ['PC5'] * n + ['PC6'] * n
    plotly_df_pos = pd.concat([top_features(0,n)[0],top_features(1,n)[0],top_features(2,n)[0],top_features(3,n)[0],top_features(4,n)[0],top_features(5,n)[0]], axis = 0)        
    plotly_df_pos['PCs'] = PCs
    plotly_df_pos = plotly_df_pos.reset_index(drop=True)
    plotly_df_pos['feats_rank'] = [1,2,3,4,5]*6
    
    plotly_df_neg = pd.concat([top_features(0,n)[1],top_features(1,n)[1],top_features(2,n)[1],top_features(3,n)[1],top_features(4,n)[1],top_features(5,n)[1]], axis = 0)        
    plotly_df_neg['PCs'] = PCs
    plotly_df_neg = plotly_df_neg.reset_index(drop=True)
    plotly_df_neg['feats_rank'] = [1,2,3,4,5]*6

    col1, col2 = st.columns(2)
    with col1 :
        df_pos = plotly_df_pos.sort_values(['feats_rank', 'PCs']).reset_index(drop=True)
        df_pos['feats_rank'] = df_pos['feats_rank'].astype(str)
        fig = px.bar(df_pos, x='PCs', y='explained_variance', color='feats_rank', text='Features', barmode='group', title='5 most contributing features - Positive correlation')
        st.plotly_chart(fig, theme='streamlit', use_container_width=True)
    with col2 :
        df_neg = plotly_df_neg.sort_values(['feats_rank', 'PCs']).reset_index(drop=True)
        df_neg['feats_rank'] = df_neg['feats_rank'].astype(str)
        fig = px.bar(df_neg, x='PCs', y='explained_variance', color='feats_rank', text='Features', barmode='group', title='5 most contributing features - Negative correlation')
        st.plotly_chart(fig, theme='streamlit', use_container_width=True)

