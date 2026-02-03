![Status](https://img.shields.io/badge/Surveillance-Active-green)
![Python](https://img.shields.io/badge/Python-3.8%2B-blue)
![Bioinformatics](https://img.shields.io/badge/Bio-BioPython-orange)
![License](https://img.shields.io/badge/license-MIT-lightgrey)

## 🧬 Executive Summary
**Pf-GVP** is a modular bioinformatics pipeline designed to automate the detection of drug-resistant mutations in *Plasmodium falciparum*. Developed to support molecular surveillance initiatives (such as the **MSMT2** project), this tool bridges the gap between raw next-generation sequencing (NGS) data and actionable public health reports.

The pipeline processes raw read data, performs quality control and alignment, identifies Single Nucleotide Polymorphisms (SNPs), and cross-references them against a curated database of known resistance markers (e.g., **Kelch13**, **PfCRT**, **PfMDR1**).

## 🚀 Key Features

* **Automated Quality Control**: Filters raw FASTQ reads based on Phred quality scores and read length to ensure high-fidelity analysis.
* **Variant Calling Engine**: Identifies SNPs and Indels in target genes relative to the *P. falciparum* 3D7 reference genome.
* **Drug Resistance Classifier**: Automatically flags mutations associated with resistance to:
    * **Artemisinin** (e.g., PfK13 C580Y)
    * **Chloroquine** (e.g., PfCRT K76T)
    * **Multi-drug resistance** (PfMDR1)
* **Structured SQL Storage**: Archives all processed samples and variants in a relational database (SQLite) for longitudinal surveillance and trend analysis.
* **Visual Reporting**: Includes a Jupyter Notebook dashboard for visualizing mutation frequency and resistance prevalence.

## 🔬 Scientific Background
Molecular surveillance is critical for tracking the spread of partial Artemisinin resistance in Tanzania. This pipeline focuses on the **Kelch13 (K13)** propeller domain.

* **Input**: Raw sequencing reads (FASTQ) from clinical isolates.
* **Logic**: The pipeline translates nucleotide changes (e.g., `G` to `A` at pos 1739) into amino acid changes (e.g., `C580Y`), which determines the resistance phenotype.

## 🛠️ Tech Stack

| Component | Technology |
| :--- | :--- |
| **Language** | Python 3.9 |
| **Bioinformatics** | BioPython (SeqIO, Align) |
| **Data Processing** | Pandas, NumPy |
| **Database** | SQLite3 (Relational Data Management) |
| **Visualization** | Matplotlib, Seaborn |
| **Workflow** | Modular Python Scripts |

## 🏗️ Pipeline Architecture

```mermaid
graph TD;
    A[Raw FASTQ Data] -->|preprocessing.py| B(QC & Filtering);
    B -->|alignment.py| C(Reference Mapping);
    C -->|variant_calling.py| D(SNP Detection);
    D -->|annotation.py| E(Functional Annotation);
    E -->|resistance_detection.py| F{Check Resistance DB};
    F -->|Yes| G[Flagged: RESISTANT];
    F -->|No| H[Flagged: WILD TYPE];
    G --> I[SQL Database];
    H --> I[SQL Database];
    I --> J[Dashboard & Reports];
⚡ Getting Started
Prerequisites
Python 3.8 or higher

Jupyter Notebook (for the demo)

Installation
Clone the repository

Bash
git clone [https://github.com/Mike-Sanga/genomic-variant-pipeline.git](https://github.com/Mike-Sanga/genomic-variant-pipeline.git)
cd genomic-variant-pipeline
Install dependencies

Bash
pip install -r requirements.txt
Usage
1. Run the Full Pipeline We have included a synthetic data generator to simulate a sequencer run containing the C580Y mutation.

Bash
python run_pipeline.py
Expected Output:

Plaintext
[1] Generating synthetic data containing C580Y mutation...
[2] Running Quality Control...
[3] Aligning reads to Pf Kelch13 Reference...
[4] Calling Variants...
[5] Annotating and Checking Drug Resistance Database...
[6] Storing results in SQL Database...

Detected Resistance Markers:
 -> C580Y: Resistant to Artemisinin
2. Visual Analysis Open the notebook to view the molecular epidemiology dashboard.

Bash
jupyter notebook notebooks/analysis_demo.ipynb
📂 Project Structure
Plaintext
genomic-variant-pipeline/
├── data/
│   ├── raw/                  # FASTQ files
│   └── reference/            # Pf3D7 Reference Genome
├── src/
│   ├── alignment.py          # Sequence alignment algorithms
│   ├── variant_calling.py    # SNP detection logic
│   ├── resistance_detection.py # Database of WHO markers
│   └── database.py           # SQL Interface
├── notebooks/                # Visualization & Reporting
├── run_pipeline.py           # Main execution script
└── requirements.txt
🛡️ Note on Implementation
This project uses a pure-Python implementation of alignment and calling algorithms for portability and demonstration purposes. In a production high-throughput environment, this pipeline would wrap standard binaries like BWA-MEM, GATK, and SnpEff, orchestrated via this Python framework.

👤 Author
Mike Levison Sanga

Role: Laboratory Scientist & Data Scientist

Specialization: Molecular Biology, Genomics, ISO 15189, Python/R

Email: mykiie85@gmail.com

GitHub: mykiie85

LinkedIn: Mike Sanga