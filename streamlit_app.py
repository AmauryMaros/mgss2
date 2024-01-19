from st_pages import Page, show_pages
import streamlit as st


show_pages(
    [   
        Page("streamlit_app.py", "Home"),
        # Page("page1_mgss_clustering.py", "MgSs clustering"),
        Page("page2_mgss_analysis.py", "MgSS Analysis"),
        Page("page3_mgCST_clustering.py", "MgCSTs Analysis"),
        Page("page4_metabolomic_analysis.py", "Metabolomic Analysis"),
        Page("page5_mgcst_classifier.py", "MgST Classifier")

    ] 
)

st.image("Medias/Holm_Lab_Logo.png")

st.title("MgCSTs project")

st.markdown('<div style="text-align: justify;">A Lactobacillus-dominated vaginal microbiome provides the first line of defense against adverse genital tract health outcomes.\
            However, there is limited understanding of the mechanisms by which the vaginal microbiome modulates protection, \
            as prior work mostly described its composition through morphologic assessment and marker gene sequencing methods that do not capture functional information. \
            To address this gap, we developed metagenomic community state types (mgCSTs) \
            which use metagenomic sequences to describe and define vaginal microbiomes based on both composition and functional potential.</div>', unsafe_allow_html = True)

st.write("https://www.jbholmlab.org/")
