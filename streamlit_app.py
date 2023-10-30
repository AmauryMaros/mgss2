from st_pages import Page, show_pages
import streamlit as st


show_pages(
    [
        Page("page1_mgss_clustering.py", "MgSs clustering"),
        Page("page2_mgss_analysis.py", "MgSS Analysis"),
        Page("page3_mgCST_clustering.py", "MgCSTs Analysis"),
        Page("page4_metabolomic_analysis.py", "Metabolomic Analysis")

    ] 
)

st.title("MgCSTs project")

