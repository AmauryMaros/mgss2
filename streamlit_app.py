from st_pages import Page, show_pages, Section
import streamlit as st


show_pages(
    [
        Page("pages/1_clustering.py", "MgSs clustering"),
        Page("pages/2_mgss.py", "MgSS Analysis"),
        Page("pages/3_mgCST_clustering.py", "MgCSTs Analysis")

    ] 
)

st.title("MgCSTs project")

