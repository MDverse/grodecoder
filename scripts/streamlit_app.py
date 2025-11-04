import tempfile
from pathlib import Path
from textwrap import fill

import streamlit as st
from PIL import Image

import grodecoder as gd



def load_image(image_path):
    return Image.open(image_path)


def show_segments(decoded):
    segments = decoded.inventory.segments
    if segments:
        st.markdown("### 🧩 Segments (Proteins/Nucleic Acids)")
        st.info(f"🔬 Found **{len(segments)}** biological segment{'s' if len(segments) != 1 else ''}")

        for i, seg in enumerate(segments):
            num_atoms = seg.number_of_atoms
            num_residues = seg.number_of_residues
            sequence = seg.sequence

            st.markdown(f"#### Segment #{i + 1}: {seg.molecular_type.value.capitalize()}")

            # Metrics in columns
            col1, col2 = st.columns(2)
            col1.markdown("⚛️ Atoms")
            col1.write(f"{num_atoms:,}")

            col2.markdown("📏 Residues")
            col2.write(f"{num_residues:,}")

            with st.expander("🧬 View Sequence", expanded=False):
                formatted_seq = fill(sequence, width=60)
                st.code(formatted_seq, language="text")

    else:
        st.warning("🔍 No segments detected in this structure.")


def show_small_molecules(decoded):
    small_molecules = decoded.inventory.small_molecules
    if small_molecules:
        st.markdown("### 🧪 Small Molecules, Ions, Solvents, Lipids")
        mol_data = [
            {
                "Name": mol.name,
                "Type": mol.molecular_type,
                "Description": mol.description,
                "#Atoms": mol.number_of_atoms,
                "#Residues": mol.number_of_residues,
            }
            for mol in small_molecules
        ]
        st.dataframe(mol_data, width="content")
        st.caption(f"Total small molecules: {len(small_molecules)}")
    else:
        st.info("No small molecules detected.")


def show_inventory(decoded):
    st.subheader("🗂️ Decoded Structure Inventory")

    show_segments(decoded)
    st.markdown("---")
    show_small_molecules(decoded)


st.set_page_config(
    page_title="GroDecoder",
    page_icon="🧬",
    layout="wide",
)

"""
# 🧬 GroDecoder  Structure Explorer
"""

"""
""" # Adds a little vertical space

"""
Upload a structure file (`.gro`, `.pdb`, `.coor`, `.crd`) and explore its decoded molecular structure.
"""

# ========================================================
# Sidebar
# ========================================================
st.sidebar.image(load_image("assets/grodecoder_logo.png"), use_container_width=True)
st.sidebar.title("Grodecoder")
st.sidebar.markdown("[Source code](https://github.com/pierrepo/grodecoder)")


# ========================================================
# Main Interface
# ========================================================
uploaded_file = st.file_uploader(
    "Choose a structure file",
    type=["gro", "pdb", "coor", "crd"],
    help="Supported formats: .gro, .pdb, .coor, .crd",
)

bond_threshold = st.slider(
    "Interchain bond detection threshold (Å)",
    min_value=1.0,
    max_value=10.0,
    value=5.0,
    step=0.1
)

compact_serialization = st.toggle("Compact serialization (no atom indices)", value=True)

if uploaded_file:
    with tempfile.NamedTemporaryFile(delete=False, suffix=Path(uploaded_file.name).suffix) as tmp:
        tmp.write(uploaded_file.read())
        tmp_path = Path(tmp.name)

    st.info(f"Processing file: `{uploaded_file.name}`")
    try:
        decoded = gd.decode_structure(tmp_path, bond_threshold=bond_threshold)
    except Exception as e:
        st.error(f"Error decoding file: {e}")
    else:
        serialization_mode = "compact" if compact_serialization else "full"
        json_data = decoded.model_dump_json(indent=2, context={"serialization_mode": serialization_mode})

        st.success("Decoding successful! 🎉")
        st.download_button(
            label="Download JSON result",
            data=json_data,
            file_name=f"{Path(uploaded_file.name).stem}_decoded.json",
            mime="application/json",
        )

        show_inventory(decoded)


    finally:
        tmp_path.unlink(missing_ok=True)
else:
    st.warning("Please upload a structure file to begin.")
