import streamlit as st
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction, molecular_weight

st.markdown(
    """
    <style>
    /* Main Layout Background and Font */
    .stApp {
        background: linear-gradient(135deg, #0F172A 0%, #1E293B 100%);
        font-family: 'Inter', system-ui, sans-serif;
    }

    /* Global Typography overrides */
    h1, h2, h3, p, span, label {
        color: #F8FAFC !important;
    }

    /* Premium Title Styling */
    .main-title {
        font-size: 2.8rem !important;
        font-weight: 800 !important;
        background: linear-gradient(90deg, #38BDF8, #818CF8);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
        margin-bottom: 0.2rem;
    }

    .subtitle {
        color: #94A3B8 !important;
        font-size: 1.1rem;
        margin-bottom: 2rem;
    }

    /* Glassmorphism Text Area Input Container */
    div[data-testid="stTextArea"] textarea {
        background-color: rgba(30, 41, 59, 0.7) !important;
        border: 1px solid rgba(255, 255, 255, 0.1) !important;
        border-radius: 12px !important;
        color: #E2E8F0 !important;
        box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.1);
        transition: all 0.3s ease;
    }
    div[data-testid="stTextArea"] textarea:focus {
        border-color: #6366F1 !important;
        box-shadow: 0 0 0 2px rgba(99, 102, 241, 0.2) !important;
    }

    /* Glowing Primary Button Styling */
    div[data-testid="stButton"] button {
        background: linear-gradient(90deg, #6366F1, #4F46E5) !important;
        color: white !important;
        border: none !important;
        border-radius: 10px !important;
        padding: 0.6rem 2.5rem !important;
        font-weight: 600 !important;
        font-size: 1rem !important;
        box-shadow: 0 4px 12px rgba(99, 102, 241, 0.3);
        transition: transform 0.2s, box-shadow 0.2s !important;
        width: 100%;
    }
    div[data-testid="stButton"] button:hover {
        transform: translateY(-2px);
        box-shadow: 0 6px 20px rgba(99, 102, 241, 0.4);
    }

    /* Beautiful Dashboard Metrics Cards */
    div[data-testid="stMetricValue"] {
        font-size: 1.8rem !important;
        font-weight: 700 !important;
        color: #38BDF8 !important;
    }
    div[data-testid="stMetricLabel"] p {
        color: #94A3B8 !important;
        text-transform: uppercase;
        font-size: 0.75rem !important;
        letter-spacing: 0.05em;
    }
    </style>
    """,
    unsafe_allow_html=True
)

# Apps title
st.title("🧬 DNA Sequence Analyzer")
st.markdown("Paste your DNA sequence below (eg. ATGC.....) to analyze it")

# Text input box
sequence_input = st.text_area("Paste DNA Sequence", height=200, placeholder="ATGC.....")

# --- 🚀 ADDED INTERACTIVE TRANSLATION CONFIGURATION 🚀 ---
st.markdown("### ⚙️ Translation Settings")
translation_mode = st.radio(
    "Choose how you want to handle translation for long sequences:",
    [
        "Complete Sequence (Show stop codons as '*')",
        "First Coding Region Only (Find ATG Start → Translate until Stop Codon)"
    ]
)

# Clean the input (remove spaces and newlines)
clean_seq = "".join(sequence_input.split()).upper()

# Run Button at the bottom
if st.button("Run Analysis"):
    if clean_seq:
        # Input Validation check
        valid_bases = set("ATGCN")
        if not set(clean_seq).issubset(valid_bases):
            st.error("❌ Invalid sequence! Please use only standard A, T, G, and C characters.")
        else:
            dna_seq = Seq(clean_seq)

            # Calculations
            seq_len = len(dna_seq)
            counts = {base: dna_seq.count(base) for base in "ATGC"}
            gc_content = gc_fraction(dna_seq) * 100
            rna_seq = dna_seq.transcribe()
            mol_wt = molecular_weight(dna_seq, seq_type="DNA")

            # Handle Translation Strategies based on user choice
            if "Complete Sequence" in translation_mode:
                # Pad to avoid incomplete codon terminal warnings
                pad_len = (3 - (len(dna_seq) % 3)) % 3
                padded_seq = dna_seq + ("N" * pad_len)
                # Keep to_stop=False to show entire sequence length with asterisks
                prot_seq = padded_seq.translate(to_stop=False)
                caption_text = "Complete translated amino acid sequence. '*' represents structural Stop Codons."
            else:
                # Find the actual start index of coding gene frame (ATG)
                start_idx = dna_seq.find("ATG")
                if start_idx != -1:
                    coding_dna = dna_seq[start_idx:]
                    prot_seq = coding_dna.translate(to_stop=True)
                    caption_text = f"Targeted translation starting at codon index position {start_idx} (ATG) ending at the first biological stop marker."
                else:
                    prot_seq = "⚠️ No standard ATG Start Codon found in this sequence sequence string!"
                    caption_text = "Unable to extract a standard target reading frame."

            # Display results
            st.header("Analysis Results")

            # Row 1: Key Metrics
            m1, m2, m3 = st.columns(3)
            m1.metric("Sequence Length", f"{seq_len:,} bp")
            m2.metric("GC Content", f"{gc_content:.2f}%")
            m3.metric("Molecular Weight", f"{mol_wt:,.2f} Da")

            col1, col2 = st.columns(2)
            with col1:
                st.subheader("Nucleotide counts")
                st.write(f"**A:** {counts['A']:,}  |  **G:** {counts['G']:,}")
                st.write(f"**T:** {counts['T']:,}  |  **C:** {counts['C']:,}")

            with col2:
                st.subheader("Composition Chart")
                chart_data = pd.DataFrame(list(counts.items()), columns=["Base", "Count"])
                st.bar_chart(chart_data, x="Base", y="Count", color="#4CAF50")

            # Final Downstream Sequencing Blocks
            st.markdown("#### 🔄 Downstream Expressions")

            tab1, tab2 = st.tabs(["🧬 Transcription (RNA)", "🧪 Translation (Protein)"])

            with tab1:
                st.caption("RNA transcript generated by substituting thymine (T) with uracil (U).")
                st.html(f"""
                    <div style="
                        background-color: #0B0F19; 
                        color: #4ADE80; 
                        padding: 16px; 
                        border-radius: 10px; 
                        font-family: 'Courier New', monospace; 
                        font-size: 14px;
                        white-space: pre-wrap; 
                        word-wrap: break-word;
                        border: 1px solid rgba(255,255,255,0.05);
                        letter-spacing: 1px;
                        max-height: 300px;
                        overflow-y: auto;
                    ">{str(rna_seq)}</div>
                """)

            with tab2:
                st.caption(caption_text)
                st.html(f"""
                    <div style="
                        background-color: #0B0F19; 
                        color: #F43F5E; 
                        padding: 16px; 
                        border-radius: 10px; 
                        font-family: 'Courier New', monospace; 
                        font-size: 14px;
                        white-space: pre-wrap; 
                        word-wrap: break-word;
                        border: 1px solid rgba(255,255,255,0.05);
                        letter-spacing: 1px;
                        max-height: 300px;
                        overflow-y: auto;
                    ">{str(prot_seq)}</div>
                """)
    else:
        st.warning("⚠️ Please provide a valid DNA sequence block first.")
