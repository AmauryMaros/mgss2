from st_pages import Page, show_pages
import streamlit as st


show_pages(
    [
        Page("1_clustering.py", "MgSs clustering"),
        Page("2_mgss.py", "MgSS Analysis"),
        Page("3_mgCST_clustering.py", "MgCSTs Analysis"),
        Page("4_metabolomic_analysis.py", "Metabolomic Analysis")

    ] 
)

st.title("MgCSTs project")

