import pandas as pd
import streamlit as st
import altair as alt
from PIL import Image
from Bio.SeqUtils import ProtParam
from Bio.Seq import Seq
from Bio.SeqUtils import molecular_weight
from Bio import pairwise2

# Mapping of amino acid codes to full names
AA_NAMES = {
    'A': 'Alanine', 'R': 'Arginine', 'N': 'Asparagine', 'D': 'Aspartic Acid',
    'C': 'Cysteine', 'E': 'Glutamic Acid', 'Q': 'Glutamine', 'G': 'Glycine',
    'H': 'Histidine', 'I': 'Isoleucine', 'L': 'Leucine', 'K': 'Lysine',
    'M': 'Methionine', 'F': 'Phenylalanine', 'P': 'Proline', 'S': 'Serine',
    'T': 'Threonine', 'W': 'Tryptophan', 'Y': 'Tyrosine', 'V': 'Valine'
}

# Page Title and Image
st.title('Bioinformatics Tool Suite')

image = Image.open('dna-logo.jpg')
st.image(image, use_column_width=True)

st.write("""
# Bioinformatics Tool Suite

This app performs various bioinformatics analyses including DNA nucleotide count, protein secondary structure prediction, codon optimization, and additional functionalities.

***
""")

# Input Text Box
sequence_input = ">DNA Query 2\nGAACACGTGGAGGCAAACAGGAAGGTGAAGAAGAACTTATCCTATCAGGACGGAAGGTCCTGTGCTCGGG\nATCTTCCAGACGTCGCGACTCTAAATTGCCCCCTCTGAGGTCAAGGAACACAAGATGGTTTTGGAAATGC\nTGAACCCGATACATTATAACATCACCAGCATCGTGCCTGAAGCCATGCCTGCTGCCACCATGCCAGTCCT"

sequence_type = st.radio("Select:", (
"DNA Sequence Counter", "Protein Secondary Structure", "Codon Optimization", "Global Alignment", "Local Alignment", "Isoelectric Point"))

sequence = st.text_area("Sequence input", sequence_input, height=250)

if sequence_type == "DNA Sequence Counter":
    sequence = sequence.splitlines()[1:]
    sequence = ''.join(sequence)
    st.write("""
    ***
    """)
    st.header('INPUT (DNA Query)')
    st.write(sequence)

    st.header('OUTPUT (DNA Nucleotide Count)')

    def DNA_nucleotide_count(seq):
        d = dict([
            ('A', seq.count('A')),
            ('T', seq.count('T')),
            ('G', seq.count('G')),
            ('C', seq.count('C'))
        ])
        return d

    X = DNA_nucleotide_count(sequence)
    st.write(X)

    # Add GC content calculation
    gc_content = (X['G'] + X['C']) / sum(X.values())
    st.write(f"GC Content: {gc_content:.2f}")

    # Add Translate DNA to Protein functionality for DNA sequence
    if st.checkbox("Translate DNA to Protein"):
        seq_obj = Seq(sequence)
        protein_seq = seq_obj.translate()
        st.header('OUTPUT (Translated Protein Sequence)')
        st.write(str(protein_seq))

    # Display DataFrame
    st.subheader('DNA Nucleotide Count DataFrame')
    df_dna = pd.DataFrame.from_dict(X, orient='index', columns=['Count'])

    # Display the nucleotide counts as a bar chart
    chart_dna = alt.Chart(df_dna.reset_index()).mark_bar().encode(
        x='index:O',
        y='Count:Q'
    ).properties(
        width=alt.Step(80)  # controls width of bar.
    )
    st.altair_chart(chart_dna)

    # Add Sequence Reverse and Complement functionality
    if st.checkbox("Reverse and Complement Sequence"):
        reversed_complemented_seq = sequence[::-1].translate(str.maketrans('ATGC', 'TACG'))
        st.header('OUTPUT (Reverse Complemented Sequence)')
        st.write(reversed_complemented_seq)

elif sequence_type == "Protein Secondary Structure":
    st.write("""
    ***
    """)
    st.header('INPUT (Protein Query)')
    st.write(sequence)

    st.header('OUTPUT (Protein Secondary Structure Prediction)')
    seq_obj = Seq(sequence)
    prot_param = ProtParam.ProteinAnalysis(str(seq_obj))
    prot_seq_length = len(sequence)
    aa_count = prot_param.count_amino_acids()
    st.write(f"Protein Sequence Length: {prot_seq_length} amino acids")
    st.write("Amino Acid Counts:")
    aa_count_with_names = {AA_NAMES.get(aa, aa): count for aa, count in aa_count.items()}
    st.write(aa_count_with_names)

    # Calculate amino acid frequencies
    total_aa = sum(aa_count.values())
    aa_freq = {AA_NAMES.get(aa, aa): count / total_aa for aa, count in aa_count.items()}
    st.write("Amino Acid Frequencies:")
    st.write(aa_freq)

    # Calculate and display molecular weight
    mol_weight = molecular_weight(sequence)
    st.write(f"Molecular Weight: {mol_weight:.2f} Da")

    # Perform protein secondary structure prediction (example using Bio.SeqUtils.ProtParam)
    helix, turn, sheet = prot_param.secondary_structure_fraction()
    st.write(f"Predicted Secondary Structure:")
    st.write(f"Helix: {helix:.2f}, Turn: {turn:.2f}, Sheet: {sheet:.2f}")

    # Create a DataFrame for secondary structure fractions
    df_secondary_structure = pd.DataFrame({
        'Structure': ['Helix', 'Turn', 'Sheet'],
        'Fraction': [helix, turn, sheet]
    })

    # Display the secondary structure fractions as a bar chart
    st.subheader('Secondary Structure Prediction')
    chart_structure = alt.Chart(df_secondary_structure).mark_bar().encode(
        x='Structure',
        y='Fraction'
    ).properties(
        width=alt.Step(80)  # controls width of bar.
    )
    st.altair_chart(chart_structure)

elif sequence_type == "Codon Optimization":
    st.header('Codon Optimization')
    expression_system = st.radio("Select Expression System", ("Bacteria", "Eukaryotes"))
    if st.button("Optimize Codons"):
        seq_obj = Seq(sequence)
        if expression_system == "Bacteria":
            optimized_seq = seq_obj.reverse_complement()
        else:
            optimized_seq = seq_obj.transcribe().reverse_complement()
        st.header('OUTPUT (Optimized Sequence)')
        st.write(str(optimized_seq))

elif sequence_type == "Global Alignment" or sequence_type == "Local Alignment":
    st.header(sequence_type)
    st.subheader("Enter Second Sequence:")
    sequence2 = st.text_area("Second Sequence", height=250)
    if st.button("Perform Alignment"):
        seq1 = Seq(sequence)
        seq2 = Seq(sequence2)

        if sequence_type == "Global Alignment":
            alignments = pairwise2.align.globalxx(seq1, seq2)
        elif sequence_type == "Local Alignment":
            alignments = pairwise2.align.localxx(seq1, seq2)

        aligned_seq1 = alignments[0][0]
        aligned_seq2 = alignments[0][1]
        alignment_score = alignments[0][2]

        st.write("Aligned Sequences:")
        st.write(aligned_seq1)
        st.write(aligned_seq2)
        st.write(f"Alignment Score: {alignment_score}")

        # Visualize alignment
        st.subheader("Alignment Visualization")
        alignment_df = pd.DataFrame({
            'Aligned Sequence 1': list(aligned_seq1),
            'Aligned Sequence 2': list(aligned_seq2)
        })
        st.write(alignment_df)

elif sequence_type == "Isoelectric Point":
    st.header('Isoelectric Point Calculation')
    seq_obj = Seq(sequence)
    prot_param = ProtParam.ProteinAnalysis(str(seq_obj))
    pI = prot_param.isoelectric_point()
    st.write(f"Isoelectric Point (pI): {pI:.2f}")
