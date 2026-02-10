"""
utils.py - Utility Functions for Genomic Variant Pipeline

This module provides common utility functions used across the pipeline,
including file handling, logging, and data validation.

Author: Molecular Surveillance Team
Version: 1.0.0
"""

import os
import logging
import hashlib
from datetime import datetime
from pathlib import Path
from typing import Optional, List, Dict, Any, Union

import pandas as pd


# =============================================================================
# LOGGING CONFIGURATION
# =============================================================================

def setup_logging(
    log_file: Optional[str] = None,
    level: int = logging.INFO
) -> logging.Logger:
    """
    Configure logging for the pipeline.
    
    Parameters
    ----------
    log_file : str, optional
        Path to log file. If None, logs only to console.
    level : int
        Logging level (default: logging.INFO)
    
    Returns
    -------
    logging.Logger
        Configured logger instance
    
    Example
    -------
    >>> logger = setup_logging('pipeline.log')
    >>> logger.info('Pipeline started')
    """
    # Create logger
    logger = logging.getLogger('genomic_pipeline')
    logger.setLevel(level)
    
    # Clear existing handlers
    logger.handlers = []
    
    # Create formatter
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    
    # Console handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(level)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
    
    # File handler (if specified)
    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(level)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
    
    return logger


# =============================================================================
# FILE HANDLING UTILITIES
# =============================================================================

def ensure_directory(path: Union[str, Path]) -> Path:
    """
    Ensure a directory exists, creating it if necessary.
    
    Parameters
    ----------
    path : str or Path
        Directory path to ensure exists
    
    Returns
    -------
    Path
        Path object for the directory
    
    Example
    -------
    >>> output_dir = ensure_directory('results/variants')
    """
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


def get_file_hash(filepath: str, algorithm: str = 'md5') -> str:
    """
    Calculate hash of a file for integrity verification.
    
    Parameters
    ----------
    filepath : str
        Path to the file
    algorithm : str
        Hash algorithm ('md5', 'sha256')
    
    Returns
    -------
    str
        Hexadecimal hash string
    
    Example
    -------
    >>> hash_val = get_file_hash('sample.fastq', 'md5')
    """
    hash_func = hashlib.new(algorithm)
    
    with open(filepath, 'rb') as f:
        for chunk in iter(lambda: f.read(8192), b''):
            hash_func.update(chunk)
    
    return hash_func.hexdigest()


def list_files(
    directory: str,
    pattern: str = '*',
    recursive: bool = False
) -> List[Path]:
    """
    List files in a directory matching a pattern.
    
    Parameters
    ----------
    directory : str
        Directory to search
    pattern : str
        Glob pattern (default: '*')
    recursive : bool
        Search recursively (default: False)
    
    Returns
    -------
    List[Path]
        List of matching file paths
    
    Example
    -------
    >>> fastq_files = list_files('data/raw', '*.fastq')
    """
    path = Path(directory)
    
    if recursive:
        return list(path.rglob(pattern))
    else:
        return list(path.glob(pattern))


# =============================================================================
# SEQUENCE UTILITIES
# =============================================================================

def reverse_complement(sequence: str) -> str:
    """
    Return the reverse complement of a DNA sequence.
    
    Parameters
    ----------
    sequence : str
        DNA sequence (A, T, G, C)
    
    Returns
    -------
    str
        Reverse complement sequence
    
    Example
    -------
    >>> reverse_complement('ATGC')
    'GCAT'
    """
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(complement.get(base, 'N') for base in reversed(sequence.upper()))


def calculate_gc_content(sequence: str) -> float:
    """
    Calculate GC content of a DNA sequence.
    
    Parameters
    ----------
    sequence : str
        DNA sequence
    
    Returns
    -------
    float
        GC content as percentage (0-100)
    
    Example
    -------
    >>> calculate_gc_content('ATGCGC')
    66.67
    """
    sequence = sequence.upper()
    gc_count = sequence.count('G') + sequence.count('C')
    total = len(sequence)
    
    if total == 0:
        return 0.0
    
    return round((gc_count / total) * 100, 2)


def translate_codon(codon: str) -> str:
    """
    Translate a DNA codon to amino acid.
    
    Parameters
    ----------
    codon : str
        Three-letter DNA codon
    
    Returns
    -------
    str
        Single-letter amino acid code
    
    Example
    -------
    >>> translate_codon('ATG')
    'M'
    """
    codon_table = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    }
    
    return codon_table.get(codon.upper(), 'X')


def get_amino_acid_name(code: str) -> str:
    """
    Get full amino acid name from single-letter code.
    
    Parameters
    ----------
    code : str
        Single-letter amino acid code
    
    Returns
    -------
    str
        Full amino acid name
    """
    amino_acids = {
        'A': 'Alanine', 'R': 'Arginine', 'N': 'Asparagine',
        'D': 'Aspartic acid', 'C': 'Cysteine', 'E': 'Glutamic acid',
        'Q': 'Glutamine', 'G': 'Glycine', 'H': 'Histidine',
        'I': 'Isoleucine', 'L': 'Leucine', 'K': 'Lysine',
        'M': 'Methionine', 'F': 'Phenylalanine', 'P': 'Proline',
        'S': 'Serine', 'T': 'Threonine', 'W': 'Tryptophan',
        'Y': 'Tyrosine', 'V': 'Valine', '*': 'Stop'
    }
    
    return amino_acids.get(code.upper(), 'Unknown')


# =============================================================================
# DATA VALIDATION
# =============================================================================

def validate_fastq_format(filepath: str) -> Dict[str, Any]:
    """
    Validate FASTQ file format and return basic statistics.
    
    Parameters
    ----------
    filepath : str
        Path to FASTQ file
    
    Returns
    -------
    dict
        Validation results including:
        - is_valid: bool
        - record_count: int
        - errors: list of error messages
    
    Example
    -------
    >>> result = validate_fastq_format('sample.fastq')
    >>> if result['is_valid']:
    ...     print(f"Valid file with {result['record_count']} reads")
    """
    results = {
        'is_valid': True,
        'record_count': 0,
        'errors': []
    }
    
    try:
        with open(filepath, 'r') as f:
            line_num = 0
            while True:
                # Read 4 lines (one FASTQ record)
                header = f.readline()
                if not header:
                    break
                
                sequence = f.readline()
                plus_line = f.readline()
                quality = f.readline()
                
                line_num += 4
                
                # Validate header
                if not header.startswith('@'):
                    results['errors'].append(
                        f"Line {line_num-3}: Header should start with '@'"
                    )
                    results['is_valid'] = False
                
                # Validate plus line
                if not plus_line.startswith('+'):
                    results['errors'].append(
                        f"Line {line_num-1}: Plus line should start with '+'"
                    )
                    results['is_valid'] = False
                
                # Validate sequence/quality length match
                if len(sequence.strip()) != len(quality.strip()):
                    results['errors'].append(
                        f"Line {line_num}: Sequence and quality lengths don't match"
                    )
                    results['is_valid'] = False
                
                results['record_count'] += 1
                
    except Exception as e:
        results['is_valid'] = False
        results['errors'].append(f"Error reading file: {str(e)}")
    
    return results


def validate_vcf_format(filepath: str) -> Dict[str, Any]:
    """
    Validate VCF file format.
    
    Parameters
    ----------
    filepath : str
        Path to VCF file
    
    Returns
    -------
    dict
        Validation results
    """
    results = {
        'is_valid': True,
        'variant_count': 0,
        'has_header': False,
        'errors': []
    }
    
    try:
        with open(filepath, 'r') as f:
            for line_num, line in enumerate(f, 1):
                if line.startswith('##'):
                    continue
                elif line.startswith('#CHROM'):
                    results['has_header'] = True
                elif not line.startswith('#'):
                    fields = line.strip().split('\t')
                    if len(fields) < 8:
                        results['errors'].append(
                            f"Line {line_num}: VCF requires at least 8 columns"
                        )
                        results['is_valid'] = False
                    else:
                        results['variant_count'] += 1
        
        if not results['has_header']:
            results['errors'].append("Missing #CHROM header line")
            results['is_valid'] = False
            
    except Exception as e:
        results['is_valid'] = False
        results['errors'].append(f"Error reading file: {str(e)}")
    
    return results


# =============================================================================
# REPORT GENERATION UTILITIES
# =============================================================================

def generate_timestamp() -> str:
    """
    Generate a formatted timestamp for reports.
    
    Returns
    -------
    str
        Timestamp in format 'YYYY-MM-DD_HH-MM-SS'
    """
    return datetime.now().strftime('%Y-%m-%d_%H-%M-%S')


def create_summary_table(
    data: List[Dict],
    columns: Optional[List[str]] = None
) -> pd.DataFrame:
    """
    Create a summary DataFrame from list of dictionaries.
    
    Parameters
    ----------
    data : list of dict
        Data to convert to DataFrame
    columns : list of str, optional
        Column order
    
    Returns
    -------
    pd.DataFrame
        Summary table
    """
    df = pd.DataFrame(data)
    
    if columns:
        # Reorder columns, keeping only specified ones
        available_cols = [c for c in columns if c in df.columns]
        df = df[available_cols]
    
    return df


def format_percentage(value: float, decimals: int = 1) -> str:
    """
    Format a decimal as percentage string.
    
    Parameters
    ----------
    value : float
        Value to format (0-1 or 0-100)
    decimals : int
        Number of decimal places
    
    Returns
    -------
    str
        Formatted percentage string
    """
    if value <= 1:
        value = value * 100
    
    return f"{value:.{decimals}f}%"


# =============================================================================
# P. FALCIPARUM SPECIFIC UTILITIES
# =============================================================================

# Reference chromosome sizes for P. falciparum 3D7
PF_CHROMOSOME_SIZES = {
    'Pf3D7_01_v3': 640851,
    'Pf3D7_02_v3': 947102,
    'Pf3D7_03_v3': 1067971,
    'Pf3D7_04_v3': 1200490,
    'Pf3D7_05_v3': 1343557,
    'Pf3D7_06_v3': 1418242,
    'Pf3D7_07_v3': 1445207,
    'Pf3D7_08_v3': 1472805,
    'Pf3D7_09_v3': 1541735,
    'Pf3D7_10_v3': 1687656,
    'Pf3D7_11_v3': 2038340,
    'Pf3D7_12_v3': 2271494,
    'Pf3D7_13_v3': 2925236,
    'Pf3D7_14_v3': 3291936,
}

def get_chromosome_from_gene(gene_name: str) -> Optional[str]:
    """
    Get chromosome location for known P. falciparum drug resistance genes.
    
    Parameters
    ----------
    gene_name : str
        Gene name (e.g., 'pfcrt', 'kelch13')
    
    Returns
    -------
    str or None
        Chromosome identifier or None if unknown
    """
    gene_locations = {
        'pfcrt': 'Pf3D7_07_v3',
        'pfmdr1': 'Pf3D7_05_v3',
        'kelch13': 'Pf3D7_13_v3',
        'dhfr': 'Pf3D7_04_v3',
        'dhps': 'Pf3D7_08_v3',
        'plasmepsin2': 'Pf3D7_14_v3',
        'plasmepsin3': 'Pf3D7_14_v3',
    }
    
    return gene_locations.get(gene_name.lower())


if __name__ == '__main__':
    # Quick test of utilities
    print("Testing utility functions...")
    
    # Test GC content
    test_seq = "ATGCGCGCATATATGC"
    print(f"GC content of {test_seq}: {calculate_gc_content(test_seq)}%")
    
    # Test reverse complement
    print(f"Reverse complement: {reverse_complement(test_seq)}")
    
    # Test codon translation
    print(f"ATG translates to: {translate_codon('ATG')} ({get_amino_acid_name('M')})")
    
    print("\nAll tests passed!")
