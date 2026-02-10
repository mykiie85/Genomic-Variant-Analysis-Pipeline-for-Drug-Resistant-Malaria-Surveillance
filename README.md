# Plasmodium falciparum Genomic Variant Pipeline

## Molecular Surveillance for Antimalarial Drug Resistance

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Bioinformatics](https://img.shields.io/badge/field-bioinformatics-green.svg)]()

A comprehensive bioinformatics pipeline for detecting antimalarial drug resistance mutations in *Plasmodium falciparum* from next-generation sequencing data. Designed for use by national malaria reference laboratories and molecular surveillance programs.

---

## Table of Contents

- [Overview](#overview)
- [Scientific Background](#scientific-background)
- [Features](#features)
- [Workflow Diagram](#workflow-diagram)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Usage](#usage)
- [Pipeline Steps](#pipeline-steps)
- [Output Files](#output-files)
- [Database Schema](#database-schema)
- [Example Output](#example-output)
- [Data Sources](#data-sources)
- [Contributing](#contributing)
- [License](#license)
- [References](#references)

---

## Overview

This pipeline provides an end-to-end solution for molecular surveillance of antimalarial drug resistance. It processes raw sequencing data (FASTQ) through quality control, alignment, variant calling, annotation, and resistance classification to produce actionable reports for clinical and public health decision-making.

### Key Capabilities

- **Quality Control**: Comprehensive FASTQ quality assessment and filtering
- **Variant Detection**: SNP and INDEL calling in drug resistance genes
- **Resistance Classification**: Automated detection of WHO-validated resistance markers
- **Data Management**: SQLite database for sample tracking and surveillance data
- **Reporting**: Publication-ready visualizations and clinical reports

### Target Users

- National Malaria Reference Laboratories
- Molecular Surveillance Programs
- Research Institutions
- Public Health Organizations

---

## Scientific Background

### Malaria Drug Resistance

*Plasmodium falciparum* is responsible for the majority of malaria deaths worldwide. The emergence and spread of antimalarial drug resistance poses a major threat to malaria control efforts. Molecular surveillance enables early detection of resistance, guiding treatment policy decisions.

### Key Resistance Genes

| Gene | Chromosome | Drug(s) Affected | Key Mutations |
|------|------------|------------------|---------------|
| **pfcrt** | 7 | Chloroquine, Amodiaquine | K76T (primary) |
| **pfmdr1** | 5 | Lumefantrine, Mefloquine | N86Y, Y184F, D1246Y |
| **kelch13** | 13 | Artemisinin | C580Y, R539T, Y493H |
| **dhfr** | 4 | Pyrimethamine | N51I, C59R, S108N |
| **dhps** | 8 | Sulfadoxine | A437G, K540E, A581G |

### WHO-Validated Markers

This pipeline implements the WHO-validated molecular markers for antimalarial drug resistance surveillance, as defined in the [WHO Malaria Threat Maps](https://www.who.int/teams/global-malaria-programme/surveillance).

---

## Features

### Core Functionality

- ✅ FASTQ quality control and filtering
- ✅ Simulated alignment to Pf3D7 reference genome
- ✅ SNP/INDEL variant calling (VCF 4.2 format)
- ✅ Functional variant annotation
- ✅ Drug resistance marker detection
- ✅ SQLite database storage
- ✅ Automated report generation
- ✅ Publication-ready visualizations

### Drugs Monitored

| Drug | Class | Status |
|------|-------|--------|
| Chloroquine | Aminoquinoline | ✅ Implemented |
| Artemisinin | Artemisinin | ✅ Implemented |
| Lumefantrine | Aryl-amino alcohol | ✅ Implemented |
| Mefloquine | Aryl-amino alcohol | ✅ Implemented |
| Amodiaquine | Aminoquinoline | ✅ Implemented |
| Pyrimethamine | Antifolate | ✅ Implemented |
| Sulfadoxine | Antifolate | ✅ Implemented |
| SP (Combo) | Antifolate | ✅ Implemented |

---

## Workflow Diagram

```
┌─────────────────────────────────────────────────────────────────┐
│                    MOLECULAR SURVEILLANCE PIPELINE               │
└─────────────────────────────────────────────────────────────────┘

    ┌──────────────┐
    │   FASTQ      │  Raw sequencing data
    │   Input      │  (Illumina, ONT, etc.)
    └──────┬───────┘
           │
           ▼
    ┌──────────────┐
    │  Step 1:     │  • Read length distribution
    │  Quality     │  • GC content analysis
    │  Control     │  • Q20/Q30 statistics
    └──────┬───────┘  • Low-quality filtering
           │
           ▼
    ┌──────────────┐
    │  Step 2:     │  • Map to Pf3D7 reference
    │  Alignment   │  • Generate SAM/BAM
    │              │  • Mapping statistics
    └──────┬───────┘
           │
           ▼
    ┌──────────────┐
    │  Step 3:     │  • SNP detection
    │  Variant     │  • INDEL detection
    │  Calling     │  • VCF generation
    └──────┬───────┘  • Quality filtering
           │
           ▼
    ┌──────────────┐
    │  Step 4:     │  • Gene annotation
    │  Annotation  │  • Amino acid changes
    │              │  • Effect prediction
    └──────┬───────┘
           │
           ▼
    ┌──────────────┐
    │  Step 5:     │  • Known marker lookup
    │  Resistance  │  • Drug classification
    │  Detection   │  • Confidence scoring
    └──────┬───────┘
           │
           ├─────────────────────────────┐
           ▼                             ▼
    ┌──────────────┐              ┌──────────────┐
    │  Step 6:     │              │  Step 7:     │
    │  Database    │              │  Reports &   │
    │  Storage     │              │  Visuals     │
    └──────────────┘              └──────────────┘
           │                             │
           └──────────────┬──────────────┘
                          ▼
               ┌──────────────────┐
               │  Clinical Report │
               │  Surveillance    │
               │  Data Export     │
               └──────────────────┘
```

---

## Installation

### Prerequisites

- Python 3.8 or higher
- pip package manager
- 4GB RAM minimum
- 10GB disk space

### Step 1: Clone Repository

```bash
git clone https://github.com/mykiie85/Genomic-Variant-Analysis-Pipeline-for-Drug-Resistant-Malaria-Surveillance.git
cd genomic-variant-pipeline
```

### Step 2: Create Virtual Environment (Recommended)

```bash
python -m venv venv
source venv/bin/activate  # Linux/macOS
# or
venv\Scripts\activate     # Windows
```

### Step 3: Install Dependencies

```bash
pip install -r requirements.txt
```

### Step 4: Verify Installation

```bash
python -c "import pandas; import Bio; print('Installation successful!')"
```

### Dependencies

```
biopython>=1.79
pysam>=0.19.0
pandas>=1.3.0
numpy>=1.21.0
matplotlib>=3.4.0
seaborn>=0.11.0
jupyter>=1.0.0
sqlalchemy>=1.4.0
```

---

## Quick Start

### Run with Demo Data

```bash
# Generate synthetic test data and run pipeline
python run_pipeline.py --generate-test-data -o ./output

# Or use provided demo data
python run_pipeline.py -i data/raw/demo_sample.fastq -o ./output
```

### View Results

```bash
# Check output directory
ls -la output/reports/

# View resistance report
cat output/reports/*_resistance_report.txt
```

---

## Usage

### Basic Usage

```bash
python run_pipeline.py -i <input.fastq> -o <output_dir>
```

### Command Line Options

| Option | Short | Description | Default |
|--------|-------|-------------|---------|
| `--input` | `-i` | Input FASTQ file | Required |
| `--output` | `-o` | Output directory | `./output` |
| `--sample-id` | `-s` | Sample identifier | Auto-generated |
| `--database` | `-d` | SQLite database path | `<output>/surveillance.db` |
| `--min-quality` | `-q` | Minimum base quality | 20 |
| `--min-length` | `-l` | Minimum read length | 50 |
| `--verbose` | `-v` | Enable verbose output | False |
| `--generate-test-data` | | Generate synthetic data | False |
| `--num-reads` | | Reads for synthetic data | 1000 |

### Examples

```bash
# Standard analysis
python run_pipeline.py -i sample.fastq -o results/

# High-stringency filtering
python run_pipeline.py -i sample.fastq -o results/ -q 30 -l 100

# Verbose output with custom sample ID
python run_pipeline.py -i sample.fastq -o results/ -s "SAMPLE_001" -v

# Generate 5000 synthetic reads for testing
python run_pipeline.py --generate-test-data --num-reads 5000 -o test_output/
```

### Using as Python Library

```python
from src.preprocessing import FASTQProcessor
from src.variant_calling import VariantCaller
from src.resistance_detection import ResistanceDetector

# Process FASTQ
processor = FASTQProcessor("sample.fastq")
stats = processor.get_quality_stats()
print(f"Total reads: {stats.total_reads}")

# Detect resistance
detector = ResistanceDetector()
profile = detector.analyze_vcf("variants.vcf", "SAMPLE_001")

if profile.has_any_resistance:
    print(f"Resistant to: {', '.join(profile.resistant_drugs)}")
```

---

## Pipeline Steps

### Step 1: FASTQ Quality Control

Analyzes raw sequencing data quality and filters low-quality reads.

**Metrics calculated:**
- Total reads and bases
- Read length distribution
- GC content (P. falciparum ~19%)
- Per-base quality scores
- Q20/Q30 percentages

**Filtering criteria:**
- Minimum quality score (default: Q20)
- Minimum read length (default: 50 bp)
- Maximum N content (default: 10%)

### Step 2: Sequence Alignment

Maps filtered reads to the *P. falciparum* 3D7 reference genome.

**Output:**
- SAM/BAM alignment files
- Mapping statistics
- Coverage per gene

**Note:** In demo mode, alignment is simulated. For production use, integrate with BWA-MEM or similar aligner.

### Step 3: Variant Calling

Identifies SNPs and INDELs in aligned reads.

**Output:**
- VCF 4.2 format
- Quality-filtered variants
- Allele frequencies

**Quality metrics:**
- Depth of coverage
- Variant quality score
- Strand bias

### Step 4: Variant Annotation

Annotates variants with functional information.

**Annotations include:**
- Gene name and ID
- Amino acid change
- Codon change
- Effect type (missense, synonymous, etc.)
- Known resistance marker flag

### Step 5: Drug Resistance Detection

Classifies samples by resistance status for each drug.

**Classification levels:**
- **SENSITIVE**: No resistance markers detected
- **RESISTANT**: Primary resistance marker present
- **PARTIAL**: Secondary markers only
- **UNKNOWN**: Insufficient data

**Confidence scoring:**
- Based on marker validation level
- Allele frequency consideration
- Multiple marker agreement

### Step 6: Database Storage

Stores all results in SQLite database for querying and surveillance.

**Tables:**
- `samples`: Sample metadata
- `variants`: Called variants
- `annotations`: Functional annotations
- `resistance_calls`: Drug classifications
- `gene_coverage`: Coverage statistics

### Step 7: Report Generation

Produces clinical and surveillance reports.

**Reports include:**
- Sample summary
- QC metrics
- Variant statistics
- Drug-by-drug resistance status
- Visualizations

---

## Output Files

```
output/
├── qc/
│   ├── SAMPLE_qc_report.txt        # Quality control report
│   └── SAMPLE_filtered.fastq       # Filtered reads
├── alignment/
│   ├── SAMPLE.sam                  # Alignment file
│   └── SAMPLE.bam                  # Binary alignment
├── variants/
│   ├── SAMPLE.vcf                  # Raw variants
│   └── SAMPLE_filtered.vcf         # Quality-filtered
├── annotation/
│   └── SAMPLE_annotated.vcf        # Annotated variants
├── reports/
│   ├── SAMPLE_resistance_report.txt  # Resistance report
│   └── SAMPLE_summary.txt            # Full summary
├── surveillance.db                   # SQLite database
└── pipeline.log                      # Execution log
```

---

## Database Schema

### Entity Relationship

```
┌─────────────┐       ┌─────────────┐       ┌──────────────┐
│   samples   │───┬───│  variants   │───────│ annotations  │
└─────────────┘   │   └─────────────┘       └──────────────┘
                  │
                  ├───│ resistance_ │
                  │   │   calls     │
                  │   └─────────────┘
                  │
                  └───│ gene_       │
                      │  coverage   │
                      └─────────────┘
```

### Key Tables

**samples**
```sql
sample_id, collection_date, location, country, 
sequencing_platform, coverage, processing_status
```

**variants**
```sql
chromosome, position, ref_allele, alt_allele, 
quality, depth, allele_frequency, genotype
```

**resistance_calls**
```sql
drug, status, confidence, mutations_found, 
primary_mutation, recommendations
```

### Example Queries

```sql
-- Get all resistant samples
SELECT s.sample_id, r.drug, r.mutations_found
FROM samples s
JOIN resistance_calls r ON s.sample_id = r.sample_id
WHERE r.status = 'resistant';

-- Mutation frequency by country
SELECT s.country, a.aa_change, COUNT(*) as count
FROM samples s
JOIN annotations a ON s.sample_id = a.sample_id
WHERE a.is_resistance_mutation = 1
GROUP BY s.country, a.aa_change;
```

---

## Example Output

### Resistance Report

```
======================================================================
PLASMODIUM FALCIPARUM DRUG RESISTANCE REPORT
======================================================================

Sample ID:     SYNTHETIC_TEST_001
Analysis Date: 2024-01-15 10:30:00

----------------------------------------------------------------------
DRUG RESISTANCE SUMMARY
----------------------------------------------------------------------

Drug                 Status       Confidence   Mutations
----------------------------------------------------------------------
Chloroquine          RESISTANT    0.95         K76T, M74I
Artemisinin          SENSITIVE    0.90         None
Lumefantrine         PARTIAL      0.75         Y184F
Mefloquine           SENSITIVE    0.90         None
Amodiaquine          RESISTANT    0.85         N86Y
Pyrimethamine        RESISTANT    0.95         N51I, C59R, S108N
Sulfadoxine          RESISTANT    0.90         A437G, K540E
SP (Combo)           RESISTANT    0.95         Quintuple mutant

----------------------------------------------------------------------
DETECTED MUTATIONS
----------------------------------------------------------------------

Gene      Position  Mutation  AF      Drug(s) Affected
----------------------------------------------------------------------
pfcrt     76        K76T      0.95    Chloroquine
pfcrt     74        M74I      0.82    Chloroquine
pfmdr1    86        N86Y      0.78    Amodiaquine
pfmdr1    184       Y184F     0.65    Lumefantrine
dhfr      51        N51I      0.92    Pyrimethamine
dhfr      59        C59R      0.88    Pyrimethamine
dhfr      108       S108N     0.95    Pyrimethamine
dhps      437       A437G     0.90    Sulfadoxine
dhps      540       K540E     0.85    Sulfadoxine

----------------------------------------------------------------------
INTERPRETATION
----------------------------------------------------------------------

⚠️  CHLOROQUINE RESISTANCE: K76T detected - primary marker

⚠️  SP RESISTANCE: Quintuple mutant detected (dhfr triple + dhps double)
    Consider alternative treatment options.

✓  ARTEMISININ SENSITIVE: No kelch13 propeller mutations detected

======================================================================
```

### Sample Visualization

![Mutation Frequency Plot](docs/images/mutation_frequency.png)
*Example mutation frequency across samples*

![Resistance Distribution](docs/images/resistance_pie.png)
*Drug resistance distribution*

---

## Data Sources

### Reference Genome

- **Pf3D7 Reference**: [PlasmoDB](https://plasmodb.org)
- **Version**: Release 3

### Public Datasets

The pipeline is compatible with data from:

- **[Pf3k Project](https://www.malariagen.net/pf3k)**: >3,000 P. falciparum genomes
- **[MalariaGEN](https://www.malariagen.net/)**: Global genomic surveillance data
- **[NCBI SRA](https://www.ncbi.nlm.nih.gov/sra)**: Raw sequencing data repository

### Downloading Pf3k Data

```bash
# Example: Download sample from ENA
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR000/ERR000001/ERR000001.fastq.gz
gunzip ERR000001.fastq.gz

# Run pipeline
python run_pipeline.py -i ERR000001.fastq -o results/
```

---

## Jupyter Notebook

An interactive tutorial notebook is available at `notebooks/analysis_demo.ipynb`.

### Launch Notebook

```bash
jupyter notebook notebooks/analysis_demo.ipynb
```

### Contents

1. **Introduction**: Pipeline overview and objectives
2. **Data Loading**: Working with FASTQ files
3. **Quality Control**: QC analysis and visualization
4. **Variant Analysis**: Understanding variant calls
5. **Resistance Detection**: Interpreting resistance results
6. **Database Queries**: Surveillance data analysis
7. **Visualization**: Creating publication figures

---

## How This Supports Molecular Surveillance

### National Programs

This pipeline supports national malaria control programs by:

1. **Standardized Analysis**: Consistent methodology across sites
2. **Rapid Turnaround**: Automated processing reduces delays
3. **Data Integration**: SQLite database enables trend analysis
4. **Policy Support**: Evidence for treatment guideline updates

### WHO Alignment

The pipeline implements WHO-recommended molecular markers and follows WHO guidelines for antimalarial drug efficacy surveillance.

### Integration Options

- **DHIS2**: Export data for national health information systems
- **REDCap**: Clinical data capture integration
- **GenBank**: Sequence submission formatting

---

## Troubleshooting

### Common Issues

**Memory Error during alignment**
```bash
# Reduce batch size
python run_pipeline.py -i sample.fastq -o output/ --batch-size 10000
```

**Database locked**
```bash
# Close other connections or copy database
cp surveillance.db surveillance_backup.db
```

**Low mapping rate**
- Check FASTQ quality scores
- Verify organism (P. falciparum expected)
- Consider mixed infection

---

## Contributing

We welcome contributions! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

### Development Setup

```bash
git clone https://github.com/your-org/genomic-variant-pipeline.git
cd genomic-variant-pipeline
pip install -e ".[dev]"
pytest tests/
```

---

## License

This project is licensed under the MIT License - see [LICENSE](LICENSE) file.

---

## References

1. WHO. *Report on antimalarial drug efficacy, resistance and response*. 2020.
2. Ariey F, et al. A molecular marker of artemisinin-resistant *Plasmodium falciparum* malaria. *Nature*. 2014.
3. Fidock DA, et al. Mutations in the P. falciparum digestive vacuole transmembrane protein PfCRT and evidence for their role in chloroquine resistance. *Mol Cell*. 2000.
4. MalariaGEN Plasmodium falciparum Community Project. *Nature*. 2021.

---

## Contact


- **Email**: mykiie85@gmail.com


---

*Built for malaria elimination.*
