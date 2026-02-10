-- ============================================================================
-- Plasmodium falciparum Molecular Surveillance Database Schema
-- ============================================================================
-- 
-- This schema supports storage and querying of:
-- - Sample metadata
-- - Genomic variants (SNPs, INDELs)
-- - Functional annotations
-- - Drug resistance classifications
-- - Gene coverage statistics
--
-- Designed for SQLite but compatible with PostgreSQL with minor modifications.
--
-- Author: Molecular Surveillance Laboratory
-- Version: 1.0.0
-- ============================================================================

-- Enable foreign key constraints (SQLite)
PRAGMA foreign_keys = ON;

-- ============================================================================
-- TABLE: samples
-- ============================================================================
-- Stores metadata for each sequenced sample
-- Primary table for linking all analysis results

CREATE TABLE IF NOT EXISTS samples (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    
    -- Sample identification
    sample_id TEXT UNIQUE NOT NULL,
    external_id TEXT,                    -- External reference (e.g., lab ID)
    
    -- Collection information
    collection_date DATE,
    collection_site TEXT,                -- Specific collection site/facility
    location TEXT,                       -- Region/district
    country TEXT,
    latitude REAL,
    longitude REAL,
    
    -- Patient information (de-identified)
    patient_age INTEGER,                 -- Age in years
    patient_sex TEXT CHECK(patient_sex IN ('M', 'F', 'U')),
    
    -- Clinical information
    parasitemia REAL,                    -- Parasites per microliter
    clinical_status TEXT,                -- e.g., 'uncomplicated', 'severe'
    treatment_history TEXT,              -- Previous antimalarial treatments
    
    -- Sequencing information
    sequencing_platform TEXT,            -- e.g., 'Illumina MiSeq'
    sequencing_date DATE,
    library_prep TEXT,                   -- Library preparation method
    
    -- Quality metrics
    total_reads INTEGER,
    reads_mapped INTEGER,
    mean_quality REAL,
    coverage REAL,                       -- Mean genome coverage
    
    -- Processing status
    processing_status TEXT DEFAULT 'pending' 
        CHECK(processing_status IN ('pending', 'processing', 'completed', 'failed')),
    analysis_date TIMESTAMP,
    
    -- Metadata
    notes TEXT,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

-- ============================================================================
-- TABLE: variants
-- ============================================================================
-- Stores called genomic variants (SNPs and INDELs)

CREATE TABLE IF NOT EXISTS variants (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    
    -- Link to sample
    sample_id TEXT NOT NULL,
    
    -- Genomic position
    chromosome TEXT NOT NULL,            -- e.g., 'Pf3D7_07_v3'
    position INTEGER NOT NULL,           -- 1-based position
    
    -- Variant alleles
    ref_allele TEXT NOT NULL,            -- Reference allele
    alt_allele TEXT NOT NULL,            -- Alternate allele
    
    -- Quality metrics
    quality REAL,                        -- Phred-scaled quality score
    depth INTEGER,                       -- Read depth at position
    allele_frequency REAL,               -- Variant allele frequency
    
    -- Genotype
    genotype TEXT,                       -- e.g., '0/1', '1/1'
    genotype_quality INTEGER,
    
    -- Variant type
    variant_type TEXT CHECK(variant_type IN ('SNP', 'INS', 'DEL', 'MNP', 'COMPLEX')),
    
    -- Filter status
    filter_status TEXT DEFAULT 'PASS',   -- PASS or filter reason
    
    -- Additional INFO fields
    strand_bias REAL,
    mapping_quality REAL,
    
    -- Metadata
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    
    -- Foreign key constraint
    FOREIGN KEY (sample_id) REFERENCES samples(sample_id) ON DELETE CASCADE,
    
    -- Ensure unique variants per sample
    UNIQUE(sample_id, chromosome, position, ref_allele, alt_allele)
);

-- ============================================================================
-- TABLE: annotations
-- ============================================================================
-- Stores functional annotations for variants

CREATE TABLE IF NOT EXISTS annotations (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    
    -- Link to variant
    variant_id INTEGER NOT NULL,
    sample_id TEXT NOT NULL,
    
    -- Gene information
    gene TEXT,                           -- Gene symbol (e.g., 'pfcrt')
    gene_id TEXT,                        -- Gene ID (e.g., 'PF3D7_0709000')
    transcript_id TEXT,
    
    -- Functional effect
    effect TEXT,                         -- e.g., 'missense', 'synonymous'
    impact TEXT CHECK(impact IN ('HIGH', 'MODERATE', 'LOW', 'MODIFIER')),
    
    -- Amino acid change
    aa_position INTEGER,                 -- Position in protein
    aa_change TEXT,                      -- e.g., 'K76T'
    codon_change TEXT,                   -- e.g., 'AAA>ACA'
    
    -- Coding information
    cds_position INTEGER,
    protein_length INTEGER,
    
    -- Resistance annotation
    is_resistance_mutation BOOLEAN DEFAULT FALSE,
    resistance_drugs TEXT,               -- Comma-separated drug names
    resistance_evidence TEXT,            -- e.g., 'WHO validated', 'candidate'
    
    -- External references
    dbsnp_id TEXT,
    pubmed_ids TEXT,                     -- Comma-separated PubMed IDs
    
    -- Metadata
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    
    -- Foreign key constraints
    FOREIGN KEY (variant_id) REFERENCES variants(id) ON DELETE CASCADE,
    FOREIGN KEY (sample_id) REFERENCES samples(sample_id) ON DELETE CASCADE
);

-- ============================================================================
-- TABLE: resistance_calls
-- ============================================================================
-- Stores drug resistance classifications per sample

CREATE TABLE IF NOT EXISTS resistance_calls (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    
    -- Link to sample
    sample_id TEXT NOT NULL,
    
    -- Drug information
    drug TEXT NOT NULL,                  -- Drug name
    drug_class TEXT,                     -- e.g., 'AMINOQUINOLINE', 'ARTEMISININ'
    
    -- Resistance classification
    status TEXT NOT NULL 
        CHECK(status IN ('sensitive', 'resistant', 'partial', 'unknown')),
    
    -- Confidence and evidence
    confidence REAL,                     -- Confidence score (0-1)
    mutations_found TEXT,                -- Comma-separated list of mutations
    primary_mutation TEXT,               -- Primary resistance marker
    
    -- Additional information
    interpretation TEXT,                 -- Clinical interpretation
    recommendations TEXT,                -- Treatment recommendations
    
    -- Metadata
    analysis_method TEXT DEFAULT 'marker-based',
    reference_version TEXT,              -- Resistance marker database version
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    
    -- Foreign key constraint
    FOREIGN KEY (sample_id) REFERENCES samples(sample_id) ON DELETE CASCADE,
    
    -- Ensure unique drug calls per sample
    UNIQUE(sample_id, drug)
);

-- ============================================================================
-- TABLE: gene_coverage
-- ============================================================================
-- Stores per-gene coverage statistics

CREATE TABLE IF NOT EXISTS gene_coverage (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    
    -- Link to sample
    sample_id TEXT NOT NULL,
    
    -- Gene information
    gene TEXT NOT NULL,
    chromosome TEXT,
    start_pos INTEGER,
    end_pos INTEGER,
    
    -- Coverage metrics
    mean_coverage REAL,
    median_coverage REAL,
    min_coverage INTEGER,
    max_coverage INTEGER,
    
    -- Breadth metrics
    bases_covered INTEGER,               -- Bases with >= 1x coverage
    pct_covered_1x REAL,                 -- % of gene with >= 1x
    pct_covered_10x REAL,                -- % of gene with >= 10x
    pct_covered_30x REAL,                -- % of gene with >= 30x
    
    -- Quality metrics
    mean_mapping_quality REAL,
    
    -- Metadata
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    
    -- Foreign key constraint
    FOREIGN KEY (sample_id) REFERENCES samples(sample_id) ON DELETE CASCADE,
    
    -- Ensure unique gene coverage per sample
    UNIQUE(sample_id, gene)
);

-- ============================================================================
-- TABLE: reference_markers
-- ============================================================================
-- Reference database of known resistance markers

CREATE TABLE IF NOT EXISTS reference_markers (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    
    -- Gene and mutation
    gene TEXT NOT NULL,
    mutation TEXT NOT NULL,              -- e.g., 'K76T'
    aa_position INTEGER,
    
    -- Genomic location
    chromosome TEXT,
    genomic_position INTEGER,
    ref_allele TEXT,
    alt_allele TEXT,
    
    -- Drug association
    drug TEXT NOT NULL,
    resistance_level TEXT 
        CHECK(resistance_level IN ('high', 'moderate', 'low', 'candidate')),
    is_primary_marker BOOLEAN DEFAULT FALSE,
    
    -- Evidence
    evidence_level TEXT,                 -- e.g., 'WHO validated', 'in vitro'
    geographic_distribution TEXT,        -- Where marker is prevalent
    first_reported_year INTEGER,
    
    -- References
    pubmed_ids TEXT,
    who_reference TEXT,
    
    -- Metadata
    notes TEXT,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    
    -- Ensure unique markers
    UNIQUE(gene, mutation, drug)
);

-- ============================================================================
-- TABLE: analysis_runs
-- ============================================================================
-- Tracks pipeline execution history

CREATE TABLE IF NOT EXISTS analysis_runs (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    
    -- Run identification
    run_id TEXT UNIQUE NOT NULL,
    
    -- Samples processed
    sample_ids TEXT,                     -- Comma-separated sample IDs
    batch_name TEXT,
    
    -- Pipeline information
    pipeline_version TEXT,
    reference_genome TEXT DEFAULT 'Pf3D7',
    reference_version TEXT DEFAULT '3',
    
    -- Run parameters
    parameters TEXT,                     -- JSON of run parameters
    
    -- Status and timing
    status TEXT CHECK(status IN ('running', 'completed', 'failed')),
    start_time TIMESTAMP,
    end_time TIMESTAMP,
    runtime_seconds REAL,
    
    -- Error tracking
    error_message TEXT,
    
    -- Metadata
    created_by TEXT,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

-- ============================================================================
-- INDEXES
-- ============================================================================
-- Optimize common query patterns

-- Sample lookups
CREATE INDEX IF NOT EXISTS idx_samples_sample_id ON samples(sample_id);
CREATE INDEX IF NOT EXISTS idx_samples_country ON samples(country);
CREATE INDEX IF NOT EXISTS idx_samples_collection_date ON samples(collection_date);
CREATE INDEX IF NOT EXISTS idx_samples_location ON samples(location);

-- Variant queries
CREATE INDEX IF NOT EXISTS idx_variants_sample_id ON variants(sample_id);
CREATE INDEX IF NOT EXISTS idx_variants_chromosome_position ON variants(chromosome, position);
CREATE INDEX IF NOT EXISTS idx_variants_gene_position ON variants(chromosome, position);

-- Annotation queries
CREATE INDEX IF NOT EXISTS idx_annotations_variant_id ON annotations(variant_id);
CREATE INDEX IF NOT EXISTS idx_annotations_sample_id ON annotations(sample_id);
CREATE INDEX IF NOT EXISTS idx_annotations_gene ON annotations(gene);
CREATE INDEX IF NOT EXISTS idx_annotations_resistance ON annotations(is_resistance_mutation);

-- Resistance queries
CREATE INDEX IF NOT EXISTS idx_resistance_sample_id ON resistance_calls(sample_id);
CREATE INDEX IF NOT EXISTS idx_resistance_drug_status ON resistance_calls(drug, status);
CREATE INDEX IF NOT EXISTS idx_resistance_status ON resistance_calls(status);

-- Gene coverage queries
CREATE INDEX IF NOT EXISTS idx_gene_coverage_sample ON gene_coverage(sample_id);
CREATE INDEX IF NOT EXISTS idx_gene_coverage_gene ON gene_coverage(gene);

-- Reference markers
CREATE INDEX IF NOT EXISTS idx_markers_gene ON reference_markers(gene);
CREATE INDEX IF NOT EXISTS idx_markers_drug ON reference_markers(drug);

-- ============================================================================
-- VIEWS
-- ============================================================================
-- Commonly used query views

-- View: Sample resistance summary
CREATE VIEW IF NOT EXISTS v_sample_resistance_summary AS
SELECT 
    s.sample_id,
    s.collection_date,
    s.country,
    s.location,
    GROUP_CONCAT(DISTINCT CASE WHEN rc.status = 'resistant' THEN rc.drug END) as resistant_drugs,
    GROUP_CONCAT(DISTINCT CASE WHEN rc.status = 'sensitive' THEN rc.drug END) as sensitive_drugs,
    COUNT(DISTINCT CASE WHEN rc.status = 'resistant' THEN rc.drug END) as resistant_count,
    COUNT(DISTINCT CASE WHEN rc.status = 'sensitive' THEN rc.drug END) as sensitive_count
FROM samples s
LEFT JOIN resistance_calls rc ON s.sample_id = rc.sample_id
GROUP BY s.sample_id;

-- View: Resistance mutation frequencies
CREATE VIEW IF NOT EXISTS v_mutation_frequencies AS
SELECT 
    a.gene,
    a.aa_change,
    a.resistance_drugs,
    COUNT(DISTINCT a.sample_id) as sample_count,
    ROUND(100.0 * COUNT(DISTINCT a.sample_id) / 
          (SELECT COUNT(DISTINCT sample_id) FROM samples), 2) as frequency_pct
FROM annotations a
WHERE a.is_resistance_mutation = 1
GROUP BY a.gene, a.aa_change, a.resistance_drugs
ORDER BY sample_count DESC;

-- View: Geographic resistance distribution
CREATE VIEW IF NOT EXISTS v_geographic_resistance AS
SELECT 
    s.country,
    s.location,
    rc.drug,
    COUNT(DISTINCT s.sample_id) as total_samples,
    SUM(CASE WHEN rc.status = 'resistant' THEN 1 ELSE 0 END) as resistant_count,
    ROUND(100.0 * SUM(CASE WHEN rc.status = 'resistant' THEN 1 ELSE 0 END) / 
          COUNT(DISTINCT s.sample_id), 2) as resistance_pct
FROM samples s
JOIN resistance_calls rc ON s.sample_id = rc.sample_id
GROUP BY s.country, s.location, rc.drug
ORDER BY s.country, rc.drug;

-- View: Variant summary per gene
CREATE VIEW IF NOT EXISTS v_gene_variant_summary AS
SELECT 
    a.gene,
    COUNT(DISTINCT v.id) as total_variants,
    SUM(CASE WHEN v.variant_type = 'SNP' THEN 1 ELSE 0 END) as snps,
    SUM(CASE WHEN v.variant_type IN ('INS', 'DEL') THEN 1 ELSE 0 END) as indels,
    SUM(CASE WHEN a.effect = 'missense' THEN 1 ELSE 0 END) as missense,
    SUM(CASE WHEN a.is_resistance_mutation = 1 THEN 1 ELSE 0 END) as resistance_mutations
FROM variants v
JOIN annotations a ON v.id = a.variant_id
WHERE a.gene IS NOT NULL
GROUP BY a.gene
ORDER BY total_variants DESC;

-- ============================================================================
-- SEED DATA: Known Resistance Markers
-- ============================================================================
-- WHO-validated and candidate resistance markers

INSERT OR IGNORE INTO reference_markers (gene, mutation, aa_position, drug, resistance_level, is_primary_marker, evidence_level) VALUES
-- pfcrt (Chloroquine resistance)
('pfcrt', 'K76T', 76, 'Chloroquine', 'high', 1, 'WHO validated'),
('pfcrt', 'M74I', 74, 'Chloroquine', 'moderate', 0, 'WHO validated'),
('pfcrt', 'N75E', 75, 'Chloroquine', 'moderate', 0, 'WHO validated'),
('pfcrt', 'A220S', 220, 'Chloroquine', 'moderate', 0, 'Associated'),
('pfcrt', 'Q271E', 271, 'Chloroquine', 'low', 0, 'Associated'),
('pfcrt', 'N326S', 326, 'Chloroquine', 'low', 0, 'Associated'),
('pfcrt', 'I356T', 356, 'Chloroquine', 'low', 0, 'Associated'),

-- pfmdr1 (Multidrug resistance)
('pfmdr1', 'N86Y', 86, 'Chloroquine', 'moderate', 0, 'Associated'),
('pfmdr1', 'N86Y', 86, 'Amodiaquine', 'moderate', 1, 'Associated'),
('pfmdr1', 'Y184F', 184, 'Lumefantrine', 'moderate', 0, 'Associated'),
('pfmdr1', 'D1246Y', 1246, 'Mefloquine', 'moderate', 0, 'Associated'),
('pfmdr1', 'S1034C', 1034, 'Mefloquine', 'low', 0, 'Candidate'),
('pfmdr1', 'N1042D', 1042, 'Mefloquine', 'low', 0, 'Candidate'),

-- kelch13 (Artemisinin resistance) - WHO validated markers
('kelch13', 'C580Y', 580, 'Artemisinin', 'high', 1, 'WHO validated'),
('kelch13', 'R539T', 539, 'Artemisinin', 'high', 1, 'WHO validated'),
('kelch13', 'Y493H', 493, 'Artemisinin', 'high', 1, 'WHO validated'),
('kelch13', 'I543T', 543, 'Artemisinin', 'high', 1, 'WHO validated'),
('kelch13', 'R561H', 561, 'Artemisinin', 'high', 1, 'WHO validated'),
('kelch13', 'P574L', 574, 'Artemisinin', 'high', 1, 'WHO validated'),
('kelch13', 'M476I', 476, 'Artemisinin', 'high', 1, 'WHO validated'),
('kelch13', 'F446I', 446, 'Artemisinin', 'moderate', 0, 'WHO validated'),
('kelch13', 'N458Y', 458, 'Artemisinin', 'moderate', 0, 'WHO validated'),
('kelch13', 'P553L', 553, 'Artemisinin', 'moderate', 0, 'WHO validated'),
('kelch13', 'R622I', 622, 'Artemisinin', 'moderate', 0, 'WHO validated'),
('kelch13', 'A675V', 675, 'Artemisinin', 'moderate', 0, 'Associated'),

-- dhfr (Pyrimethamine resistance)
('dhfr', 'N51I', 51, 'Pyrimethamine', 'high', 0, 'WHO validated'),
('dhfr', 'C59R', 59, 'Pyrimethamine', 'high', 0, 'WHO validated'),
('dhfr', 'S108N', 108, 'Pyrimethamine', 'high', 1, 'WHO validated'),
('dhfr', 'I164L', 164, 'Pyrimethamine', 'high', 0, 'WHO validated'),

-- dhps (Sulfadoxine resistance)
('dhps', 'S436A', 436, 'Sulfadoxine', 'moderate', 0, 'WHO validated'),
('dhps', 'A437G', 437, 'Sulfadoxine', 'high', 1, 'WHO validated'),
('dhps', 'K540E', 540, 'Sulfadoxine', 'high', 0, 'WHO validated'),
('dhps', 'A581G', 581, 'Sulfadoxine', 'high', 0, 'WHO validated'),
('dhps', 'A613S', 613, 'Sulfadoxine', 'moderate', 0, 'WHO validated');

-- ============================================================================
-- TRIGGERS
-- ============================================================================
-- Automatic timestamp updates

CREATE TRIGGER IF NOT EXISTS update_samples_timestamp 
AFTER UPDATE ON samples
BEGIN
    UPDATE samples SET updated_at = CURRENT_TIMESTAMP WHERE id = NEW.id;
END;

CREATE TRIGGER IF NOT EXISTS update_reference_markers_timestamp 
AFTER UPDATE ON reference_markers
BEGIN
    UPDATE reference_markers SET updated_at = CURRENT_TIMESTAMP WHERE id = NEW.id;
END;

-- ============================================================================
-- End of Schema
-- ============================================================================
