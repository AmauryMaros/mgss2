from st_pages import Page, show_pages, Section
import streamlit as st


show_pages(
    [
        Page("./streamlit_app.py", "Home", "🏠"),
        Section(name="MgSs", icon=":beginner:"),
        Page("pages/1_clustering.py", "MgSs clustering"),
        Page("pages/2_mgss.py", "MgSS Analysis"),
        Section(name="MCSTs", icon=":beginner:"),
        Page("pages/3_mgCST_clustering.py", "MgCSTs clustering"),
        Page("pages/4_MgCSTs.py", "MgCSTs Analysis")

    ] 
)

st.title("MgCSTs project")

