import math
from io import BytesIO, StringIO

import logomaker
import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st
from Bio import SeqIO


st.set_page_config(
    page_title="Protein Sequence Logo Generator",
    page_icon="🧬",
    layout="wide",
)

st.title("🧬 Protein Sequence Logo Generator")
st.markdown(
    "Analyze aligned protein sequences and visualize amino acid "
    "conservation using sequence logos."
)

# Sample FASTA file users can download to test the app.
sample_fasta = """>Protein1
MKTLLLTLVV
>Protein2
MKTLLLALVV
>Protein3
MKTLLLTLVA
>Protein4
MKTLLLTLVV
"""

st.download_button(
    label="Download Sample FASTA",
    data=sample_fasta,
    file_name="sample_proteins.fasta",
    mime="text/plain",
)

uploaded_file = st.file_uploader(
    "Upload your protein FASTA file",
    type=["fasta", "fa", "txt"],
)

if uploaded_file is not None:
    try:
        fasta_text = uploaded_file.getvalue().decode("utf-8")
        records = list(SeqIO.parse(StringIO(fasta_text), "fasta"))

        sequences = [
            str(record.seq).upper()
            for record in records
        ]

        valid_amino_acids = set("ACDEFGHIKLMNPQRSTVWY")

        if not sequences:
            st.error("No FASTA sequences found.")

        elif any(
            not sequence or set(sequence) - valid_amino_acids
            for sequence in sequences
        ):
            st.error(
                "Invalid or empty protein sequence detected. "
                "Use the 20 standard amino-acid letters only."
            )

        elif len({len(sequence) for sequence in sequences}) != 1:
            st.error(
                "All sequences must be aligned and have equal lengths."
            )

        else:
            number_of_sequences = len(sequences)
            sequence_length = len(sequences[0])
            amino_acids = "ACDEFGHIKLMNPQRSTVWY"

            st.success(
                f"Successfully processed {number_of_sequences} sequences!"
            )

            # Display basic sequence statistics.
            col1, col2, col3 = st.columns(3)
            col1.metric("Number of Sequences", number_of_sequences)
            col2.metric("Sequence Length", sequence_length)
            col3.metric("Amino Acids", 20)

            # Calculate amino-acid frequency at every position.
            frequency_matrix = pd.DataFrame(
                0,
                index=range(sequence_length),
                columns=list(amino_acids),
                dtype=int,
            )

            for position in range(sequence_length):
                for sequence in sequences:
                    amino_acid = sequence[position]
                    frequency_matrix.loc[position, amino_acid] += 1

            st.subheader("Amino Acid Frequency")

            display_frequency = frequency_matrix.copy()
            display_frequency.index = range(1, sequence_length + 1)
            display_frequency.index.name = "Position"
            st.dataframe(display_frequency)

            # Convert frequencies into probabilities.
            probabilities = frequency_matrix.div(
                frequency_matrix.sum(axis=1),
                axis=0,
            )

            # Calculate Shannon entropy for each position.
            entropy = []
            for _, row in probabilities.iterrows():
                position_entropy = -sum(
                    probability * math.log2(probability)
                    for probability in row
                    if probability > 0
                )
                entropy.append(position_entropy)

            # Information content = maximum entropy - observed entropy.
            maximum_entropy = math.log2(20)
            information = [
                max(0, maximum_entropy - value)
                for value in entropy
            ]

            st.subheader("Information Content")
            information_df = pd.DataFrame(
                {
                    "Position": range(1, sequence_length + 1),
                    "Information Content (bits)": information,
                }
            )
            st.dataframe(information_df, hide_index=True)

            # Each letter's height is its probability multiplied by
            # the information content at that position.
            logo_matrix = probabilities.mul(information, axis=0)

            st.subheader("Generated Protein Sequence Logo")

            fig, ax = plt.subplots(
                figsize=(max(12, sequence_length * 0.8), 6)
            )

            amino_acid_colors = {
                # Hydrophobic
                "A": "#808080",
                "V": "#808080",
                "L": "#808080",
                "I": "#808080",
                "M": "#808080",
                "F": "#808080",
                "W": "#808080",
                "P": "#808080",
                # Polar
                "S": "#2CA02C",
                "T": "#2CA02C",
                "N": "#2CA02C",
                "Q": "#2CA02C",
                "C": "#2CA02C",
                "Y": "#2CA02C",
                "G": "#2CA02C",
                # Positively charged
                "K": "#1F77B4",
                "R": "#1F77B4",
                "H": "#1F77B4",
                # Negatively charged
                "D": "#D62728",
                "E": "#D62728",
            }

            logomaker.Logo(
                logo_matrix,
                ax=ax,
                color_scheme=amino_acid_colors,
            )

            ax.set_xlabel("Amino Acid Position")
            ax.set_ylabel("Information Content (bits)")
            ax.set_title("Protein Sequence Logo")
            ax.set_xticks(range(sequence_length))
            ax.set_xticklabels(range(1, sequence_length + 1))
            ax.set_ylim(0, maximum_entropy + 0.2)
            ax.grid(axis="y", alpha=0.2)

            plt.tight_layout()
            st.pyplot(fig)

            # Convert the generated plot into an in-memory PNG.
            image_buffer = BytesIO()
            fig.savefig(
                image_buffer,
                format="png",
                dpi=300,
                bbox_inches="tight",
            )
            image_buffer.seek(0)

            st.download_button(
                label="Download Sequence Logo",
                data=image_buffer.getvalue(),
                file_name="protein_sequence_logo.png",
                mime="image/png",
            )

            # Allow users to download the frequency matrix as CSV.
            st.download_button(
                label="Download Frequency Table",
                data=frequency_matrix.to_csv(
                    index_label="Position"
                ).encode("utf-8"),
                file_name="amino_acid_frequency.csv",
                mime="text/csv",
            )

            plt.close(fig)

    except UnicodeDecodeError:
        st.error("The uploaded file could not be read as UTF-8 text.")
    except Exception as error:
        st.error(f"Error processing FASTA file: {error}")
