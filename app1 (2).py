import streamlit as st
from Bio import SeqIO
import numpy as np
import seaborn as sns 
import matplotlib.pyplot as plt
import io
from scipy.cluster.hierarchy import linkage, dendrogram
from alignment_logic import global_alignment, local_alignment, get_symbol

st.set_page_config(page_title="Bioinformatics Alignment", layout="centered")
st.title("Sequence Alignment Tool")

def run_batch_alignment(sequences, match, mismatch, gap):
    n = len(sequences)
    matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            _, _, score = global_alignment(str(sequences[i].seq), str(sequences[j].seq), match, mismatch, gap)
            matrix[i, j] = matrix[j, i] = score

    return matrix

with st.sidebar:
    st.header("Parameters")
    match = st.number_input("Match Score", value=5)
    mismatch = st.number_input("Mismatch Score", value=-1)
    gap = st.number_input("Gap Penalty", value=-3)

    mode = st.radio("Alignment Mode", ("Pairwise", "Batch All-vs-All"))

input_method = st.radio("Choose input method:", ("Manual Text", "Upload FASTA File"))

def get_sequence(method, label, key):
    if method == "Manual Text":

        seq = st.text_input(f"enter {label}:",key=key)
        if seq:
            return seq.strip().upper()
        
    else:
        uploaded_file = st.file_uploader(f"Upload {label} file", type=["fasta", "fa"], key=key)
        if uploaded_file :
            
            stringio = io.StringIO(uploaded_file.getvalue().decode("utf-8"))
            try:

                record = SeqIO.read(stringio, "fasta")
                return str(record.seq).upper()
            except Exception as e:
                st.error(f"Error reading FASTA file: {e}")
                return None
    return ""

if mode == "Pairwise":
    seq1 = get_sequence(input_method, "first sequence", "seq1")
    seq2 = get_sequence(input_method, "second sequence", "seq2")
    alignment_type = st.radio("Type", ("Global (Needleman-Wunsch)", "Local (Smith-Waterman)"))
    
    if st.button("Run Alignment"):
        if not seq1 or not seq2:
            st.error("Please enter both sequences!")
        else:
            if "Global" in alignment_type:
                align1, align2, score = global_alignment(seq1, seq2, match, mismatch, gap)
            else:
                align1, align2, score = local_alignment(seq1, seq2, match, mismatch, gap)
            
            st.success(f"Alignment Score: {score}")
            symbol = get_symbol(align1, align2)
            line_length = 80
            st.subheader("Alignment Result:")
            for i in range(0, len(align1), line_length):
                st.code(f"{align1[i : i + line_length]}\n{symbol[i : i + line_length]}\n{align2[i : i + line_length]}", language="text")

elif mode == "Batch All-vs-All":
    st.subheader("Batch Alignment & Phylogeny")
    batch_file = st.file_uploader("Upload Multi-FASTA file", type=["fasta", "fa"])
    
    if batch_file:
        stringio = io.StringIO(batch_file.getvalue().decode("utf-8"))
        sequences = list(SeqIO.parse(stringio, "fasta"))
        
        if st.button("Run Batch Alignment"):
            if len(sequences) < 2:
                st.warning("Please upload a file with at least 2 sequences.")
            else:
                with st.spinner("Aligning sequences and building matrix..."):
                    matrix = run_batch_alignment(sequences, match, mismatch, gap)
                    
                    st.subheader("Similarity Matrix Heatmap")
                    fig, ax = plt.subplots(figsize=(8, 6))
                    sns.heatmap(matrix, annot=True, xticklabels=[s.id for s in sequences], 
                                yticklabels=[s.id for s in sequences], cmap="YlGnBu", ax=ax)
                    st.pyplot(fig)
                    
                    st.subheader("Phylogenetic Quick-View (Tree)")
                    try:
                        max_score = np.max(matrix)
                        dist_matrix = max_score - matrix
                        
                        from scipy.spatial.distance import squareform
                        condensed_dist = squareform(dist_matrix)
                        
                        linked = linkage(condensed_dist, method='average')
                        
                        fig_tree, ax_tree = plt.subplots(figsize=(8, 5))
                        dendrogram(linked, labels=[s.id for s in sequences], ax=ax_tree, orientation='top')
                        plt.xticks(rotation=45, ha='right')
                        st.pyplot(fig_tree)
                    except Exception as tree_error:
                        st.error(f"Could not generate tree: {tree_error}")























