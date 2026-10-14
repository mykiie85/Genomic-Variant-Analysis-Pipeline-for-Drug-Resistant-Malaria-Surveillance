"""
database.py - Database Management Module

This module handles storage and retrieval of genomic surveillance data
in a SQLite database. It provides methods for storing samples, variants,
annotations, and resistance calls.

Author: Molecular Surveillance Team
Version: 1.0.0

Database Schema:
- samples: Sample metadata
- variants: Called variants
- annotations: Variant annotations
- resistance_calls: Drug resistance classifications
"""

import os
import sqlite3
from typing import Dict, List, Optional, Tuple, Any
from datetime import datetime
from contextlib import contextmanager

from src.utils import setup_logging, ensure_directory
from src.resistance_detection import SampleResistanceProfile, ResistanceStatus


# Columns of the ``samples`` table (see initialize_database / get_schema_sql).
# Column names cannot be bound as SQL parameters, so any caller-supplied
# field name must be validated against this whitelist before being placed
# into a query string.
SAMPLE_COLUMNS = frozenset({
    'sample_id',
    'collection_date',
    'location',
    'country',
    'region',
    'latitude',
    'longitude',
    'patient_age',
    'patient_sex',
    'parasitemia',
    'sequencing_date',
    'sequencing_platform',
    'read_count',
    'mean_coverage',
    'created_at',
    'updated_at',
})


# =============================================================================
# DATABASE MANAGER CLASS
# =============================================================================

class DatabaseManager:
    """
    Manage SQLite database for genomic surveillance data.
    
    This class provides methods for creating, populating, and querying
    the surveillance database.
    
    Parameters
    ----------
    db_path : str
        Path to SQLite database file
    
    Example
    -------
    >>> db = DatabaseManager('database/surveillance.db')
    >>> db.initialize_database()
    >>> db.insert_sample('SAMPLE001', 'Tanzania', '2024-01-15')
    """
    
    def __init__(self, db_path: str):
        """Initialize database manager."""
        self.db_path = db_path
        self.logger = setup_logging()
        
        # Ensure database directory exists
        ensure_directory(os.path.dirname(db_path) or '.')
    
    @contextmanager
    def get_connection(self):
        """
        Context manager for database connections.
        
        Yields
        ------
        sqlite3.Connection
            Database connection
        """
        conn = sqlite3.connect(self.db_path)
        conn.row_factory = sqlite3.Row
        try:
            yield conn
            conn.commit()
        except Exception as e:
            conn.rollback()
            raise e
        finally:
            conn.close()
    
    def initialize_database(self) -> None:
        """
        Create database tables if they don't exist.
        
        Creates tables for:
        - samples
        - variants
        - annotations
        - resistance_calls
        - gene_coverage
        """
        self.logger.info(f"Initializing database: {self.db_path}")
        
        with self.get_connection() as conn:
            cursor = conn.cursor()
            
            # Samples table
            cursor.execute('''
                CREATE TABLE IF NOT EXISTS samples (
                    sample_id TEXT PRIMARY KEY,
                    collection_date DATE,
                    location TEXT,
                    country TEXT,
                    region TEXT,
                    latitude REAL,
                    longitude REAL,
                    patient_age INTEGER,
                    patient_sex TEXT,
                    parasitemia REAL,
                    sequencing_date DATE,
                    sequencing_platform TEXT,
                    read_count INTEGER,
                    mean_coverage REAL,
                    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
                )
            ''')
            
            # Variants table
            cursor.execute('''
                CREATE TABLE IF NOT EXISTS variants (
                    variant_id INTEGER PRIMARY KEY AUTOINCREMENT,
                    sample_id TEXT NOT NULL,
                    chromosome TEXT NOT NULL,
                    position INTEGER NOT NULL,
                    ref_allele TEXT NOT NULL,
                    alt_allele TEXT NOT NULL,
                    quality REAL,
                    filter TEXT,
                    depth INTEGER,
                    allele_frequency REAL,
                    genotype TEXT,
                    variant_type TEXT,
                    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                    FOREIGN KEY (sample_id) REFERENCES samples(sample_id),
                    UNIQUE(sample_id, chromosome, position, alt_allele)
                )
            ''')
            
            # Annotations table
            cursor.execute('''
                CREATE TABLE IF NOT EXISTS annotations (
                    annotation_id INTEGER PRIMARY KEY AUTOINCREMENT,
                    variant_id INTEGER NOT NULL,
                    gene TEXT,
                    gene_id TEXT,
                    effect TEXT,
                    codon_change TEXT,
                    amino_acid_change TEXT,
                    amino_acid_position INTEGER,
                    is_resistance_mutation BOOLEAN,
                    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                    FOREIGN KEY (variant_id) REFERENCES variants(variant_id)
                )
            ''')
            
            # Resistance calls table
            cursor.execute('''
                CREATE TABLE IF NOT EXISTS resistance_calls (
                    call_id INTEGER PRIMARY KEY AUTOINCREMENT,
                    sample_id TEXT NOT NULL,
                    drug TEXT NOT NULL,
                    status TEXT NOT NULL,
                    confidence REAL,
                    mutations_found TEXT,
                    notes TEXT,
                    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                    FOREIGN KEY (sample_id) REFERENCES samples(sample_id),
                    UNIQUE(sample_id, drug)
                )
            ''')
            
            # Gene coverage table
            cursor.execute('''
                CREATE TABLE IF NOT EXISTS gene_coverage (
                    coverage_id INTEGER PRIMARY KEY AUTOINCREMENT,
                    sample_id TEXT NOT NULL,
                    gene TEXT NOT NULL,
                    mean_coverage REAL,
                    coverage_breadth REAL,
                    variant_count INTEGER,
                    FOREIGN KEY (sample_id) REFERENCES samples(sample_id),
                    UNIQUE(sample_id, gene)
                )
            ''')
            
            # Create indexes for common queries
            cursor.execute('''
                CREATE INDEX IF NOT EXISTS idx_variants_sample 
                ON variants(sample_id)
            ''')
            cursor.execute('''
                CREATE INDEX IF NOT EXISTS idx_variants_gene 
                ON annotations(gene)
            ''')
            cursor.execute('''
                CREATE INDEX IF NOT EXISTS idx_resistance_drug 
                ON resistance_calls(drug, status)
            ''')
            cursor.execute('''
                CREATE INDEX IF NOT EXISTS idx_samples_location 
                ON samples(country, region)
            ''')
        
        self.logger.info("Database initialized successfully")
    
    def insert_sample(
        self,
        sample_id: str,
        location: str = None,
        collection_date: str = None,
        **kwargs
    ) -> None:
        """
        Insert a new sample record.
        
        Parameters
        ----------
        sample_id : str
            Unique sample identifier
        location : str, optional
            Collection location
        collection_date : str, optional
            Collection date (YYYY-MM-DD)
        **kwargs : dict
            Additional sample metadata. Keys must be columns of the
            ``samples`` table (see ``SAMPLE_COLUMNS``).
        
        Raises
        ------
        ValueError
            If a keyword argument is not a valid ``samples`` column.
        """
        with self.get_connection() as conn:
            cursor = conn.cursor()
            
            # Build insert query dynamically based on provided fields
            fields = ['sample_id']
            values = [sample_id]
            placeholders = ['?']
            
            if location:
                fields.append('location')
                values.append(location)
                placeholders.append('?')
            
            if collection_date:
                fields.append('collection_date')
                values.append(collection_date)
                placeholders.append('?')
            
            for key, value in kwargs.items():
                if value is not None:
                    # Column names are interpolated into the SQL text, so
                    # only allow known columns from the samples schema.
                    if key not in SAMPLE_COLUMNS or key in fields:
                        raise ValueError(f"Invalid sample field: {key!r}")
                    fields.append(key)
                    values.append(value)
                    placeholders.append('?')
            
            # All identifiers in `fields` come from SAMPLE_COLUMNS; all
            # values are passed separately as bound parameters.
            query = f'''
                INSERT OR REPLACE INTO samples ({', '.join(fields)})
                VALUES ({', '.join(placeholders)})
            '''
            
            cursor.execute(query, values)
        
        self.logger.info(f"Inserted sample: {sample_id}")
    
    def insert_variant(
        self,
        sample_id: str,
        chromosome: str,
        position: int,
        ref_allele: str,
        alt_allele: str,
        **kwargs
    ) -> int:
        """
        Insert a variant record.
        
        Parameters
        ----------
        sample_id : str
            Sample identifier
        chromosome : str
            Chromosome name
        position : int
            Genomic position
        ref_allele : str
            Reference allele
        alt_allele : str
            Alternative allele
        **kwargs : dict
            Additional variant fields
        
        Returns
        -------
        int
            Inserted variant ID
        """
        with self.get_connection() as conn:
            cursor = conn.cursor()
            
            cursor.execute('''
                INSERT OR REPLACE INTO variants 
                (sample_id, chromosome, position, ref_allele, alt_allele,
                 quality, filter, depth, allele_frequency, genotype, variant_type)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ''', (
                sample_id, chromosome, position, ref_allele, alt_allele,
                kwargs.get('quality'),
                kwargs.get('filter'),
                kwargs.get('depth'),
                kwargs.get('allele_frequency'),
                kwargs.get('genotype'),
                kwargs.get('variant_type')
            ))
            
            return cursor.lastrowid
    
    def insert_annotation(
        self,
        variant_id: int,
        gene: str = None,
        effect: str = None,
        amino_acid_change: str = None,
        **kwargs
    ) -> None:
        """
        Insert a variant annotation.
        
        Parameters
        ----------
        variant_id : int
            Variant ID from variants table
        gene : str, optional
            Gene name
        effect : str, optional
            Variant effect
        amino_acid_change : str, optional
            Amino acid change notation
        """
        with self.get_connection() as conn:
            cursor = conn.cursor()
            
            cursor.execute('''
                INSERT INTO annotations
                (variant_id, gene, gene_id, effect, codon_change,
                 amino_acid_change, amino_acid_position, is_resistance_mutation)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?)
            ''', (
                variant_id,
                gene,
                kwargs.get('gene_id'),
                effect,
                kwargs.get('codon_change'),
                amino_acid_change,
                kwargs.get('amino_acid_position'),
                kwargs.get('is_resistance_mutation', False)
            ))
    
    def insert_resistance_call(
        self,
        sample_id: str,
        drug: str,
        status: str,
        confidence: float = None,
        mutations_found: List[str] = None,
        notes: str = None
    ) -> None:
        """
        Insert a resistance classification.
        
        Parameters
        ----------
        sample_id : str
            Sample identifier
        drug : str
            Drug name
        status : str
            Resistance status
        confidence : float, optional
            Classification confidence
        mutations_found : list, optional
            List of mutations found
        notes : str, optional
            Additional notes
        """
        with self.get_connection() as conn:
            cursor = conn.cursor()
            
            mutations_str = ','.join(mutations_found) if mutations_found else None
            
            cursor.execute('''
                INSERT OR REPLACE INTO resistance_calls
                (sample_id, drug, status, confidence, mutations_found, notes)
                VALUES (?, ?, ?, ?, ?, ?)
            ''', (sample_id, drug, status, confidence, mutations_str, notes))
    
    def store_resistance_profile(
        self,
        profile: SampleResistanceProfile
    ) -> None:
        """
        Store a complete resistance profile.
        
        Parameters
        ----------
        profile : SampleResistanceProfile
            Resistance profile object
        """
        for drug, result in profile.drug_results.items():
            mutations = [m.mutation for m in result.mutations_found]
            
            self.insert_resistance_call(
                sample_id=profile.sample_id,
                drug=drug,
                status=str(result.status),
                confidence=result.confidence,
                mutations_found=mutations,
                notes=result.notes
            )
        
        self.logger.info(f"Stored resistance profile for {profile.sample_id}")
    
    def import_vcf(
        self,
        vcf_path: str,
        sample_id: str = None,
        **sample_metadata
    ) -> Dict[str, int]:
        """
        Import variants from VCF file into database.
        
        Parameters
        ----------
        vcf_path : str
            Path to VCF file
        sample_id : str, optional
            Sample identifier (derived from VCF if not provided)
        **sample_metadata : dict
            Additional sample metadata
        
        Returns
        -------
        dict
            Import statistics
        """
        if sample_id is None:
            sample_id = os.path.splitext(os.path.basename(vcf_path))[0]
            sample_id = sample_id.replace('_annotated', '').replace('_filtered', '')
        
        self.logger.info(f"Importing VCF: {vcf_path}")
        
        stats = {
            'variants_imported': 0,
            'annotations_imported': 0,
            'errors': 0
        }
        
        # Insert sample
        self.insert_sample(sample_id, **sample_metadata)
        
        with open(vcf_path, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                try:
                    fields = line.strip().split('\t')
                    
                    chrom = fields[0]
                    pos = int(fields[1])
                    ref = fields[3]
                    alt = fields[4]
                    qual = float(fields[5]) if fields[5] != '.' else None
                    filt = fields[6]
                    info = fields[7]
                    
                    # Parse INFO field
                    info_dict = {}
                    for item in info.split(';'):
                        if '=' in item:
                            key, value = item.split('=', 1)
                            info_dict[key] = value
                    
                    # Parse FORMAT/SAMPLE
                    genotype = '.'
                    if len(fields) > 9:
                        format_keys = fields[8].split(':')
                        sample_values = fields[9].split(':')
                        sample_dict = dict(zip(format_keys, sample_values))
                        genotype = sample_dict.get('GT', '.')
                    
                    # Determine variant type
                    if len(ref) == 1 and len(alt) == 1:
                        var_type = 'SNP'
                    elif len(ref) < len(alt):
                        var_type = 'INS'
                    else:
                        var_type = 'DEL'
                    
                    # Insert variant
                    variant_id = self.insert_variant(
                        sample_id=sample_id,
                        chromosome=chrom,
                        position=pos,
                        ref_allele=ref,
                        alt_allele=alt,
                        quality=qual,
                        filter=filt,
                        depth=int(info_dict.get('DP', 0)) if info_dict.get('DP') else None,
                        allele_frequency=float(info_dict.get('AF', 0)) if info_dict.get('AF') else None,
                        genotype=genotype,
                        variant_type=var_type
                    )
                    
                    stats['variants_imported'] += 1
                    
                    # Insert annotation if present
                    gene = info_dict.get('GENE')
                    if gene and gene != 'intergenic':
                        self.insert_annotation(
                            variant_id=variant_id,
                            gene=gene,
                            effect=info_dict.get('EFFECT'),
                            amino_acid_change=info_dict.get('AA_CHANGE'),
                            is_resistance_mutation=(info_dict.get('RESISTANCE') == 'YES')
                        )
                        stats['annotations_imported'] += 1
                
                except Exception as e:
                    self.logger.warning(f"Error importing variant: {e}")
                    stats['errors'] += 1
        
        self.logger.info(
            f"VCF import complete: {stats['variants_imported']} variants, "
            f"{stats['annotations_imported']} annotations"
        )
        
        return stats
    
    # =========================================================================
    # QUERY METHODS
    # =========================================================================
    
    def get_sample(self, sample_id: str) -> Optional[Dict]:
        """
        Retrieve sample information.
        
        Parameters
        ----------
        sample_id : str
            Sample identifier
        
        Returns
        -------
        dict or None
            Sample record
        """
        with self.get_connection() as conn:
            cursor = conn.cursor()
            cursor.execute(
                'SELECT * FROM samples WHERE sample_id = ?',
                (sample_id,)
            )
            row = cursor.fetchone()
            return dict(row) if row else None
    
    def get_sample_variants(
        self,
        sample_id: str,
        gene: str = None
    ) -> List[Dict]:
        """
        Get variants for a sample.
        
        Parameters
        ----------
        sample_id : str
            Sample identifier
        gene : str, optional
            Filter by gene name
        
        Returns
        -------
        list of dict
            Variant records with annotations
        """
        with self.get_connection() as conn:
            cursor = conn.cursor()
            
            query = '''
                SELECT v.*, a.gene, a.effect, a.amino_acid_change, a.is_resistance_mutation
                FROM variants v
                LEFT JOIN annotations a ON v.variant_id = a.variant_id
                WHERE v.sample_id = ?
            '''
            params = [sample_id]
            
            if gene:
                query += ' AND a.gene = ?'
                params.append(gene)
            
            cursor.execute(query, params)
            return [dict(row) for row in cursor.fetchall()]
    
    def get_resistance_summary(
        self,
        drug: str = None,
        country: str = None
    ) -> List[Dict]:
        """
        Get resistance summary statistics.
        
        Parameters
        ----------
        drug : str, optional
            Filter by drug
        country : str, optional
            Filter by country
        
        Returns
        -------
        list of dict
            Aggregated resistance data
        """
        with self.get_connection() as conn:
            cursor = conn.cursor()
            
            query = '''
                SELECT 
                    r.drug,
                    r.status,
                    COUNT(*) as count,
                    s.country
                FROM resistance_calls r
                JOIN samples s ON r.sample_id = s.sample_id
                WHERE 1=1
            '''
            params = []
            
            if drug:
                query += ' AND r.drug = ?'
                params.append(drug)
            
            if country:
                query += ' AND s.country = ?'
                params.append(country)
            
            query += ' GROUP BY r.drug, r.status, s.country'
            
            cursor.execute(query, params)
            return [dict(row) for row in cursor.fetchall()]
    
    def get_mutation_frequencies(
        self,
        gene: str = None
    ) -> List[Dict]:
        """
        Get mutation frequency statistics.
        
        Parameters
        ----------
        gene : str, optional
            Filter by gene
        
        Returns
        -------
        list of dict
            Mutation frequencies
        """
        with self.get_connection() as conn:
            cursor = conn.cursor()
            
            query = '''
                SELECT 
                    a.gene,
                    a.amino_acid_change,
                    COUNT(*) as count,
                    COUNT(*) * 100.0 / (SELECT COUNT(DISTINCT sample_id) FROM samples) as frequency
                FROM annotations a
                WHERE a.is_resistance_mutation = 1
            '''
            params = []
            
            if gene:
                query += ' AND a.gene = ?'
                params.append(gene)
            
            query += ' GROUP BY a.gene, a.amino_acid_change ORDER BY count DESC'
            
            cursor.execute(query, params)
            return [dict(row) for row in cursor.fetchall()]
    
    def get_all_samples(self) -> List[Dict]:
        """Get all sample records."""
        with self.get_connection() as conn:
            cursor = conn.cursor()
            cursor.execute('SELECT * FROM samples ORDER BY collection_date DESC')
            return [dict(row) for row in cursor.fetchall()]
    
    def get_statistics(self) -> Dict:
        """
        Get database statistics.
        
        Returns
        -------
        dict
            Database statistics
        """
        with self.get_connection() as conn:
            cursor = conn.cursor()
            
            stats = {}
            
            cursor.execute('SELECT COUNT(*) FROM samples')
            stats['total_samples'] = cursor.fetchone()[0]
            
            cursor.execute('SELECT COUNT(*) FROM variants')
            stats['total_variants'] = cursor.fetchone()[0]
            
            cursor.execute('SELECT COUNT(*) FROM annotations WHERE is_resistance_mutation = 1')
            stats['resistance_mutations'] = cursor.fetchone()[0]
            
            cursor.execute('SELECT COUNT(*) FROM resistance_calls WHERE status = "Resistant"')
            stats['resistant_samples'] = cursor.fetchone()[0]
            
            cursor.execute('SELECT DISTINCT drug FROM resistance_calls')
            stats['drugs_tracked'] = cursor.fetchall()
            
            return stats


# =============================================================================
# CONVENIENCE FUNCTIONS
# =============================================================================

def create_database(db_path: str) -> DatabaseManager:
    """
    Create and initialize a new database.
    
    Parameters
    ----------
    db_path : str
        Path for database file
    
    Returns
    -------
    DatabaseManager
        Initialized database manager
    """
    db = DatabaseManager(db_path)
    db.initialize_database()
    return db


def get_schema_sql() -> str:
    """
    Get SQL schema as string.
    
    Returns
    -------
    str
        SQL CREATE TABLE statements
    """
    return '''
-- Genomic Variant Pipeline Database Schema
-- Version 1.0.0

-- Samples table: stores sample metadata
CREATE TABLE IF NOT EXISTS samples (
    sample_id TEXT PRIMARY KEY,
    collection_date DATE,
    location TEXT,
    country TEXT,
    region TEXT,
    latitude REAL,
    longitude REAL,
    patient_age INTEGER,
    patient_sex TEXT,
    parasitemia REAL,
    sequencing_date DATE,
    sequencing_platform TEXT,
    read_count INTEGER,
    mean_coverage REAL,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

-- Variants table: stores called variants
CREATE TABLE IF NOT EXISTS variants (
    variant_id INTEGER PRIMARY KEY AUTOINCREMENT,
    sample_id TEXT NOT NULL,
    chromosome TEXT NOT NULL,
    position INTEGER NOT NULL,
    ref_allele TEXT NOT NULL,
    alt_allele TEXT NOT NULL,
    quality REAL,
    filter TEXT,
    depth INTEGER,
    allele_frequency REAL,
    genotype TEXT,
    variant_type TEXT,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (sample_id) REFERENCES samples(sample_id),
    UNIQUE(sample_id, chromosome, position, alt_allele)
);

-- Annotations table: stores variant annotations
CREATE TABLE IF NOT EXISTS annotations (
    annotation_id INTEGER PRIMARY KEY AUTOINCREMENT,
    variant_id INTEGER NOT NULL,
    gene TEXT,
    gene_id TEXT,
    effect TEXT,
    codon_change TEXT,
    amino_acid_change TEXT,
    amino_acid_position INTEGER,
    is_resistance_mutation BOOLEAN,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (variant_id) REFERENCES variants(variant_id)
);

-- Resistance calls table: stores drug resistance classifications
CREATE TABLE IF NOT EXISTS resistance_calls (
    call_id INTEGER PRIMARY KEY AUTOINCREMENT,
    sample_id TEXT NOT NULL,
    drug TEXT NOT NULL,
    status TEXT NOT NULL,
    confidence REAL,
    mutations_found TEXT,
    notes TEXT,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (sample_id) REFERENCES samples(sample_id),
    UNIQUE(sample_id, drug)
);

-- Gene coverage table
CREATE TABLE IF NOT EXISTS gene_coverage (
    coverage_id INTEGER PRIMARY KEY AUTOINCREMENT,
    sample_id TEXT NOT NULL,
    gene TEXT NOT NULL,
    mean_coverage REAL,
    coverage_breadth REAL,
    variant_count INTEGER,
    FOREIGN KEY (sample_id) REFERENCES samples(sample_id),
    UNIQUE(sample_id, gene)
);

-- Indexes for performance
CREATE INDEX IF NOT EXISTS idx_variants_sample ON variants(sample_id);
CREATE INDEX IF NOT EXISTS idx_variants_gene ON annotations(gene);
CREATE INDEX IF NOT EXISTS idx_resistance_drug ON resistance_calls(drug, status);
CREATE INDEX IF NOT EXISTS idx_samples_location ON samples(country, region);
'''


if __name__ == '__main__':
    print("Database Management Module")
    print("-" * 50)
    
    # Print schema
    print("\nDatabase Schema:")
    print(get_schema_sql())
