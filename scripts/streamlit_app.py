import os
from pathlib import Path
import tempfile

from loguru import logger
from PIL import Image
import streamlit as st

import grodecoder as gd


if __name__ == "__main__":
    logger.remove()
    logger.add(st.write, format="{time:YYYY-MM-DD HH:mm:ss} | {message}", level="INFO")
    logger.add("test.log", format="{time:YYYY-MM-DD HH:mm:ss} | {level} | {message}")
    st.markdown(
        """
        <style>
        p {
            margin-bottom: 0px;
        }
        </style>
        """,
        unsafe_allow_html=True,
    )
    logo = Image.open("assets/grodecoder_logo.png")
    st.sidebar.image(logo)
    st.sidebar.markdown(
        f"""
    **GroDecoder** extracts and identifies
    the molecular components of a structure file (PDB or GRO)
    issued from a molecular dynamics simulation.

    ---
    [Source code](https://github.com/pierrepo/grodecoder)  

    Last commit date: `{gd.get_git_last_commit_date()}`  
    Last commit hash: `{gd.get_git_last_commit_hash()[:7]}`
    """
    )


    st.markdown("# GroDecoder 📦")

    uploaded_file = st.file_uploader("Choose a structure file", type=["pdb", "gro"])

    st.markdown(
        """
    Examples:
    <a href="https://raw.githubusercontent.com/pierrepo/grodecoder/main/data/examples/barstar.gro" download target="_blank">Barstar</a>
    <br /><br />
    """,
        unsafe_allow_html=True,
    )

    st.write("Options:")
    check_connectivity_value = st.checkbox(
        "Check connectivity (add degree and number of edge in the fingerprint) ",
        value=False,
    )
    bond_threshold_value = st.text_input(
        "Threshold value ('auto' or any positive value):", value="auto"
    )

    if (uploaded_file is not None) and st.button("Run analysis"):
        temp_dir = tempfile.mkdtemp()
        path = Path(temp_dir) / uploaded_file.name
        with open(path, "wb") as f:
            f.write(uploaded_file.getvalue())
        st.write(f"Saved to {path}")
        with st.spinner("Running analysis..."):
            JSON_filepath = gd.main(
                path,
                check_connectivity=check_connectivity_value,
                bond_threshold=bond_threshold_value,
                query_pdb=True,
            )

        with open(JSON_filepath, "r") as file:
            file_content = file.read()
        st.toast("Analysis completed successfully 🎉", icon="🥳")
        st.download_button(
            label="Download molecular inventory",
            data=file_content,
            file_name=JSON_filepath,
            mime="text/plain",
        )