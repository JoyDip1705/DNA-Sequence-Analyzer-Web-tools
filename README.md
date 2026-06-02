🧬 DNA Sequence Analyzer

A sleek, premium web application built with Streamlit and Biopython designed for rapid DNA sequence analysis, transcription, and interactive translation.The application features a modern dark-mode user interface with a custom glassmorphism design, providing researchers and students with instant genetic metrics, interactive nucleotide visualization charts, and customizable downstream expression readouts.

✨ Key Features

    • Instant Genetic Metrics: Computes sequence length, exact molecular weight (Da), and GC/AT content percentage instantly.

    • Interactive Charting: Generates responsive, clean bar charts mapping nucleotide frequencies using native Streamlit graphing engines.

    • Smart Transcription Engine: Simulates biological transcription by processing DNA into RNA strings while handling edge cases.

    • Dual-Mode Translation Configuration:
       • Complete Sequence Mode: Translates the entire string while highlighting structural Stop Codons as *.
       • First Coding Region Mode: Scans for the biological initiation site (ATG Start Codon) and stops cleanly at the first downstream stop marker.

    • Input Sanitization & Validation: Automatically handles whitespaces, numeric debris, and case issues while flagging invalid non-nucleotide characters.

    • Premium Dark Mode UI: Custom-styled CSS containing animated elements, glowing metric markers, and clean, readable code blocks for sequence streams.

🛠️ Tech Stack

    • Core Framework: Streamlit
    • Bioinformatics Engine: Biopython (Bio.Seq)
    • Data Processing: Pandas
