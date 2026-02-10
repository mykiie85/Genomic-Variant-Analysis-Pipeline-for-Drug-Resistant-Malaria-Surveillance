"""
preprocessing.py - FASTQ Quality Control Module

This module handles quality control and preprocessing of raw FASTQ sequencing data.
It includes functions for reading sequences, calculating quality metrics,
and filtering low-quality reads.

Author: Molecular Surveillance Team
Version: 1.0.0

Typical Usage:
    from preprocessing import FASTQProcessor
    
    processor = FASTQProcessor('sample.fastq')
    stats = processor.get_quality_stats()
    processor.filter_reads(min_quality=20, output_path='filtered.fastq')
"""

import os
from typing import Dict, List, Tuple, Optional, Generator, Any
from dataclasses import dataclass, field
from collections import defaultdict

import numpy as np
import pandas as pd

from src.utils import (
    calculate_gc_content,
    setup_logging,
    ensure_directory
)


# =============================================================================
# DATA CLASSES
# =============================================================================

@dataclass
class FASTQRead:
    """
    Represents a single FASTQ read.
    
    Attributes
    ----------
    header : str
        Read identifier (starts with @)
    sequence : str
        DNA sequence
    quality : str
        Quality scores (Phred+33 encoded)
    """
    header: str
    sequence: str
    quality: str
    
    @property
    def read_id(self) -> str:
        """Extract read ID from header."""
        return self.header.split()[0].lstrip('@')
    
    @property
    def length(self) -> int:
        """Return read length."""
        return len(self.sequence)
    
    @property
    def gc_content(self) -> float:
        """Calculate GC content."""
        return calculate_gc_content(self.sequence)
    
    @property
    def mean_quality(self) -> float:
        """Calculate mean Phred quality score."""
        scores = [ord(c) - 33 for c in self.quality]
        return np.mean(scores) if scores else 0.0
    
    @property
    def min_quality(self) -> int:
        """Return minimum quality score."""
        scores = [ord(c) - 33 for c in self.quality]
        return min(scores) if scores else 0
    
    def to_fastq_string(self) -> str:
        """Convert back to FASTQ format string."""
        return f"{self.header}\n{self.sequence}\n+\n{self.quality}\n"


@dataclass
class QualityStats:
    """
    Quality statistics for a FASTQ file.
    
    Attributes
    ----------
    total_reads : int
        Total number of reads
    total_bases : int
        Total number of bases
    mean_read_length : float
        Average read length
    mean_quality : float
        Average quality score
    gc_content : float
        Overall GC content percentage
    quality_distribution : dict
        Distribution of quality scores
    length_distribution : dict
        Distribution of read lengths
    """
    total_reads: int = 0
    total_bases: int = 0
    mean_read_length: float = 0.0
    mean_quality: float = 0.0
    gc_content: float = 0.0
    min_length: int = 0
    max_length: int = 0
    q20_bases: int = 0  # Bases with quality >= 20
    q30_bases: int = 0  # Bases with quality >= 30
    quality_distribution: Dict[int, int] = field(default_factory=dict)
    length_distribution: Dict[int, int] = field(default_factory=dict)
    per_position_quality: List[float] = field(default_factory=list)
    
    @property
    def q20_percentage(self) -> float:
        """Percentage of bases with Q>=20."""
        return (self.q20_bases / self.total_bases * 100) if self.total_bases > 0 else 0.0
    
    @property
    def q30_percentage(self) -> float:
        """Percentage of bases with Q>=30."""
        return (self.q30_bases / self.total_bases * 100) if self.total_bases > 0 else 0.0
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for reporting."""
        return {
            'total_reads': self.total_reads,
            'total_bases': self.total_bases,
            'mean_read_length': round(self.mean_read_length, 2),
            'mean_quality': round(self.mean_quality, 2),
            'gc_content': round(self.gc_content, 2),
            'min_length': self.min_length,
            'max_length': self.max_length,
            'q20_percentage': round(self.q20_percentage, 2),
            'q30_percentage': round(self.q30_percentage, 2)
        }


# =============================================================================
# FASTQ PROCESSOR CLASS
# =============================================================================

class FASTQProcessor:
    """
    Process and analyze FASTQ sequencing data.
    
    This class provides methods for reading FASTQ files, calculating
    quality metrics, and filtering reads based on quality thresholds.
    
    Parameters
    ----------
    filepath : str
        Path to input FASTQ file
    
    Attributes
    ----------
    filepath : str
        Path to the FASTQ file
    stats : QualityStats
        Quality statistics (populated after calling get_quality_stats)
    
    Example
    -------
    >>> processor = FASTQProcessor('sample.fastq')
    >>> stats = processor.get_quality_stats()
    >>> print(f"Total reads: {stats.total_reads}")
    >>> processor.filter_reads(min_quality=20, output_path='filtered.fastq')
    """
    
    def __init__(self, filepath: str):
        """Initialize FASTQ processor."""
        self.filepath = filepath
        self.stats: Optional[QualityStats] = None
        self.logger = setup_logging()
        
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"FASTQ file not found: {filepath}")
    
    def read_fastq(self) -> Generator[FASTQRead, None, None]:
        """
        Read FASTQ file and yield reads one at a time.
        
        Yields
        ------
        FASTQRead
            Individual FASTQ read objects
        
        Example
        -------
        >>> for read in processor.read_fastq():
        ...     print(f"Read {read.read_id}: {read.length} bp")
        """
        with open(self.filepath, 'r') as f:
            while True:
                header = f.readline().strip()
                if not header:
                    break
                
                sequence = f.readline().strip()
                _ = f.readline()  # Plus line
                quality = f.readline().strip()
                
                yield FASTQRead(
                    header=header,
                    sequence=sequence,
                    quality=quality
                )
    
    def get_quality_stats(self) -> QualityStats:
        """
        Calculate comprehensive quality statistics for the FASTQ file.
        
        Returns
        -------
        QualityStats
            Object containing all quality metrics
        
        Example
        -------
        >>> stats = processor.get_quality_stats()
        >>> print(f"Mean quality: {stats.mean_quality}")
        >>> print(f"GC content: {stats.gc_content}%")
        """
        self.logger.info(f"Calculating quality statistics for {self.filepath}")
        
        stats = QualityStats()
        
        all_qualities = []
        all_lengths = []
        all_gc = []
        position_qualities = defaultdict(list)
        
        for read in self.read_fastq():
            stats.total_reads += 1
            stats.total_bases += read.length
            
            all_lengths.append(read.length)
            all_gc.append(read.gc_content)
            
            # Quality scores for this read
            for i, char in enumerate(read.quality):
                q_score = ord(char) - 33
                all_qualities.append(q_score)
                position_qualities[i].append(q_score)
                
                # Count Q20 and Q30 bases
                if q_score >= 20:
                    stats.q20_bases += 1
                if q_score >= 30:
                    stats.q30_bases += 1
            
            # Update length distribution
            length_bin = (read.length // 10) * 10
            stats.length_distribution[length_bin] = \
                stats.length_distribution.get(length_bin, 0) + 1
        
        # Calculate summary statistics
        if stats.total_reads > 0:
            stats.mean_read_length = np.mean(all_lengths)
            stats.mean_quality = np.mean(all_qualities)
            stats.gc_content = np.mean(all_gc)
            stats.min_length = min(all_lengths)
            stats.max_length = max(all_lengths)
            
            # Quality distribution
            for q in range(0, 45):
                stats.quality_distribution[q] = all_qualities.count(q)
            
            # Per-position quality
            max_pos = max(position_qualities.keys()) + 1
            stats.per_position_quality = [
                np.mean(position_qualities[i]) if i in position_qualities else 0
                for i in range(max_pos)
            ]
        
        self.stats = stats
        self.logger.info(f"Processed {stats.total_reads} reads, {stats.total_bases} bases")
        
        return stats
    
    def filter_reads(
        self,
        output_path: str,
        min_quality: float = 20.0,
        min_length: int = 50,
        max_n_fraction: float = 0.1
    ) -> Dict[str, int]:
        """
        Filter reads based on quality thresholds.
        
        Parameters
        ----------
        output_path : str
            Path for filtered output FASTQ
        min_quality : float
            Minimum mean quality score (default: 20)
        min_length : int
            Minimum read length (default: 50)
        max_n_fraction : float
            Maximum fraction of N bases allowed (default: 0.1)
        
        Returns
        -------
        dict
            Filtering statistics:
            - total_reads: int
            - passed_reads: int
            - failed_quality: int
            - failed_length: int
            - failed_n_content: int
        
        Example
        -------
        >>> results = processor.filter_reads(
        ...     'filtered.fastq',
        ...     min_quality=25,
        ...     min_length=100
        ... )
        >>> print(f"Passed: {results['passed_reads']} / {results['total_reads']}")
        """
        self.logger.info(f"Filtering reads with Q>={min_quality}, length>={min_length}")
        
        # Ensure output directory exists
        ensure_directory(os.path.dirname(output_path) or '.')
        
        results = {
            'total_reads': 0,
            'passed_reads': 0,
            'failed_quality': 0,
            'failed_length': 0,
            'failed_n_content': 0
        }
        
        with open(output_path, 'w') as out:
            for read in self.read_fastq():
                results['total_reads'] += 1
                
                # Check length
                if read.length < min_length:
                    results['failed_length'] += 1
                    continue
                
                # Check quality
                if read.mean_quality < min_quality:
                    results['failed_quality'] += 1
                    continue
                
                # Check N content
                n_fraction = read.sequence.upper().count('N') / read.length
                if n_fraction > max_n_fraction:
                    results['failed_n_content'] += 1
                    continue
                
                # Write passing read
                out.write(read.to_fastq_string())
                results['passed_reads'] += 1
        
        pass_rate = results['passed_reads'] / results['total_reads'] * 100 \
            if results['total_reads'] > 0 else 0
        
        self.logger.info(
            f"Filtering complete: {results['passed_reads']}/{results['total_reads']} "
            f"reads passed ({pass_rate:.1f}%)"
        )
        
        return results
    
    def trim_reads(
        self,
        output_path: str,
        trim_5prime: int = 0,
        trim_3prime: int = 0,
        quality_trim: bool = True,
        quality_threshold: int = 20
    ) -> Dict[str, int]:
        """
        Trim reads from ends and by quality.
        
        Parameters
        ----------
        output_path : str
            Path for trimmed output FASTQ
        trim_5prime : int
            Bases to trim from 5' end (default: 0)
        trim_3prime : int
            Bases to trim from 3' end (default: 0)
        quality_trim : bool
            Perform quality-based 3' trimming (default: True)
        quality_threshold : int
            Quality threshold for trimming (default: 20)
        
        Returns
        -------
        dict
            Trimming statistics
        """
        self.logger.info(f"Trimming reads: 5'={trim_5prime}, 3'={trim_3prime}")
        
        results = {
            'total_reads': 0,
            'trimmed_reads': 0,
            'bases_removed': 0,
            'reads_discarded': 0
        }
        
        min_length = 30  # Minimum length after trimming
        
        with open(output_path, 'w') as out:
            for read in self.read_fastq():
                results['total_reads'] += 1
                original_length = read.length
                
                # Fixed trimming
                start = trim_5prime
                end = read.length - trim_3prime if trim_3prime > 0 else read.length
                
                # Quality-based 3' trimming
                if quality_trim:
                    while end > start:
                        q_score = ord(read.quality[end-1]) - 33
                        if q_score >= quality_threshold:
                            break
                        end -= 1
                
                # Check if read is long enough
                if end - start < min_length:
                    results['reads_discarded'] += 1
                    continue
                
                # Create trimmed read
                trimmed_seq = read.sequence[start:end]
                trimmed_qual = read.quality[start:end]
                
                trimmed_read = FASTQRead(
                    header=read.header,
                    sequence=trimmed_seq,
                    quality=trimmed_qual
                )
                
                out.write(trimmed_read.to_fastq_string())
                results['trimmed_reads'] += 1
                results['bases_removed'] += original_length - len(trimmed_seq)
        
        self.logger.info(
            f"Trimming complete: {results['trimmed_reads']} reads, "
            f"{results['bases_removed']} bases removed"
        )
        
        return results
    
    def generate_report(self, output_path: str) -> str:
        """
        Generate a text quality report.
        
        Parameters
        ----------
        output_path : str
            Path for output report file
        
        Returns
        -------
        str
            Path to generated report
        """
        if self.stats is None:
            self.get_quality_stats()
        
        report_lines = [
            "=" * 60,
            "FASTQ QUALITY CONTROL REPORT",
            "=" * 60,
            f"\nInput file: {self.filepath}",
            f"\n{'Metric':<30} {'Value':>20}",
            "-" * 52,
            f"{'Total reads':<30} {self.stats.total_reads:>20,}",
            f"{'Total bases':<30} {self.stats.total_bases:>20,}",
            f"{'Mean read length':<30} {self.stats.mean_read_length:>20.1f}",
            f"{'Min read length':<30} {self.stats.min_length:>20}",
            f"{'Max read length':<30} {self.stats.max_length:>20}",
            f"{'Mean quality score':<30} {self.stats.mean_quality:>20.1f}",
            f"{'GC content (%)':<30} {self.stats.gc_content:>20.1f}",
            f"{'Q20 bases (%)':<30} {self.stats.q20_percentage:>20.1f}",
            f"{'Q30 bases (%)':<30} {self.stats.q30_percentage:>20.1f}",
            "\n" + "=" * 60,
            "END OF REPORT",
            "=" * 60
        ]
        
        report_text = '\n'.join(report_lines)
        
        with open(output_path, 'w') as f:
            f.write(report_text)
        
        self.logger.info(f"Quality report saved to {output_path}")
        
        return output_path


# =============================================================================
# CONVENIENCE FUNCTIONS
# =============================================================================

def process_fastq(
    input_path: str,
    output_dir: str,
    min_quality: float = 20.0,
    min_length: int = 50
) -> Tuple[QualityStats, Dict[str, int]]:
    """
    Convenience function to process a FASTQ file in one step.
    
    Parameters
    ----------
    input_path : str
        Path to input FASTQ file
    output_dir : str
        Directory for output files
    min_quality : float
        Minimum quality threshold
    min_length : int
        Minimum read length
    
    Returns
    -------
    tuple
        (QualityStats, filtering_results dict)
    
    Example
    -------
    >>> stats, results = process_fastq(
    ...     'raw/sample.fastq',
    ...     'processed/',
    ...     min_quality=25
    ... )
    """
    # Setup output paths
    ensure_directory(output_dir)
    
    base_name = os.path.splitext(os.path.basename(input_path))[0]
    filtered_path = os.path.join(output_dir, f"{base_name}_filtered.fastq")
    report_path = os.path.join(output_dir, f"{base_name}_qc_report.txt")
    
    # Process
    processor = FASTQProcessor(input_path)
    stats = processor.get_quality_stats()
    
    filter_results = processor.filter_reads(
        filtered_path,
        min_quality=min_quality,
        min_length=min_length
    )
    
    processor.generate_report(report_path)
    
    return stats, filter_results


if __name__ == '__main__':
    # Example usage and testing
    print("FASTQ Preprocessing Module")
    print("-" * 40)
    
    # This would be run with actual data
    # processor = FASTQProcessor('data/raw/sample.fastq')
    # stats = processor.get_quality_stats()
    # print(stats.to_dict())
