"""
alignment.py - Sequence Alignment Module

This module handles read alignment to the Plasmodium falciparum reference genome.
It provides simulated alignment functionality for demonstration purposes and
interfaces with real alignment tools when available.

Author: Molecular Surveillance Team
Version: 1.0.0

Note: This module simulates alignment for demonstration. In production,
      use BWA, Bowtie2, or minimap2 for real alignments.
"""

import os
import random
from typing import Dict, List, Optional, Tuple, Generator
from dataclasses import dataclass, field
from datetime import datetime

from src.utils import (
    setup_logging,
    ensure_directory,
    reverse_complement,
    PF_CHROMOSOME_SIZES
)


# =============================================================================
# DATA CLASSES
# =============================================================================

@dataclass
class AlignmentRecord:
    """
    Represents a single read alignment (SAM record).
    
    Attributes
    ----------
    query_name : str
        Read identifier
    flag : int
        SAM flag (bitwise)
    reference_name : str
        Reference chromosome name
    position : int
        1-based leftmost mapping position
    mapq : int
        Mapping quality score
    cigar : str
        CIGAR string describing alignment
    sequence : str
        Read sequence
    quality : str
        Quality scores
    """
    query_name: str
    flag: int
    reference_name: str
    position: int
    mapq: int
    cigar: str
    sequence: str
    quality: str
    tags: Dict[str, str] = field(default_factory=dict)
    
    @property
    def is_mapped(self) -> bool:
        """Check if read is mapped."""
        return not (self.flag & 4)
    
    @property
    def is_reverse(self) -> bool:
        """Check if read is on reverse strand."""
        return bool(self.flag & 16)
    
    @property
    def is_secondary(self) -> bool:
        """Check if this is a secondary alignment."""
        return bool(self.flag & 256)
    
    @property
    def mapping_length(self) -> int:
        """Calculate alignment length from CIGAR."""
        import re
        length = 0
        for match in re.finditer(r'(\d+)([MIDNSHP=X])', self.cigar):
            count = int(match.group(1))
            op = match.group(2)
            if op in 'MDN=X':
                length += count
        return length
    
    def to_sam_line(self) -> str:
        """Convert to SAM format line."""
        tags_str = '\t'.join(f"{k}:{v}" for k, v in self.tags.items())
        
        fields = [
            self.query_name,
            str(self.flag),
            self.reference_name,
            str(self.position),
            str(self.mapq),
            self.cigar,
            '*',  # RNEXT (mate reference)
            '0',  # PNEXT (mate position)
            '0',  # TLEN (template length)
            self.sequence,
            self.quality
        ]
        
        line = '\t'.join(fields)
        if tags_str:
            line += '\t' + tags_str
        
        return line


@dataclass
class AlignmentStats:
    """
    Statistics from an alignment run.
    
    Attributes
    ----------
    total_reads : int
        Total reads processed
    mapped_reads : int
        Reads that aligned
    unmapped_reads : int
        Reads that did not align
    mapping_rate : float
        Percentage of reads mapped
    mean_mapq : float
        Average mapping quality
    """
    total_reads: int = 0
    mapped_reads: int = 0
    unmapped_reads: int = 0
    unique_mapped: int = 0
    multi_mapped: int = 0
    mean_mapq: float = 0.0
    coverage_by_chrom: Dict[str, float] = field(default_factory=dict)
    
    @property
    def mapping_rate(self) -> float:
        """Calculate mapping rate percentage."""
        if self.total_reads == 0:
            return 0.0
        return (self.mapped_reads / self.total_reads) * 100
    
    def to_dict(self) -> Dict:
        """Convert to dictionary."""
        return {
            'total_reads': self.total_reads,
            'mapped_reads': self.mapped_reads,
            'unmapped_reads': self.unmapped_reads,
            'mapping_rate': round(self.mapping_rate, 2),
            'unique_mapped': self.unique_mapped,
            'multi_mapped': self.multi_mapped,
            'mean_mapq': round(self.mean_mapq, 2)
        }


# =============================================================================
# REFERENCE GENOME HANDLER
# =============================================================================

class ReferenceGenome:
    """
    Handler for P. falciparum reference genome.
    
    This class manages reference genome data and provides methods
    for coordinate lookup and sequence retrieval.
    
    Parameters
    ----------
    reference_path : str, optional
        Path to reference FASTA file
    
    Note
    ----
    For demonstration, this uses simulated reference data.
    In production, load the actual Pf3D7 reference genome.
    """
    
    # Gene coordinates for key drug resistance genes (Pf3D7 v3)
    GENE_COORDINATES = {
        'pfcrt': {
            'chromosome': 'Pf3D7_07_v3',
            'start': 403612,
            'end': 406318,
            'strand': '+',
            'description': 'Chloroquine resistance transporter'
        },
        'pfmdr1': {
            'chromosome': 'Pf3D7_05_v3',
            'start': 957885,
            'end': 962144,
            'strand': '+',
            'description': 'Multidrug resistance protein 1'
        },
        'kelch13': {
            'chromosome': 'Pf3D7_13_v3',
            'start': 1724817,
            'end': 1726997,
            'strand': '-',
            'description': 'Kelch protein K13 (artemisinin resistance)'
        },
        'dhfr': {
            'chromosome': 'Pf3D7_04_v3',
            'start': 748068,
            'end': 749914,
            'strand': '+',
            'description': 'Dihydrofolate reductase-thymidylate synthase'
        },
        'dhps': {
            'chromosome': 'Pf3D7_08_v3',
            'start': 548200,
            'end': 550616,
            'strand': '+',
            'description': 'Dihydropteroate synthase'
        }
    }
    
    def __init__(self, reference_path: Optional[str] = None):
        """Initialize reference genome handler."""
        self.reference_path = reference_path
        self.chromosomes = PF_CHROMOSOME_SIZES
        self.logger = setup_logging()
        
        # Load reference if path provided
        if reference_path and os.path.exists(reference_path):
            self._load_reference()
        else:
            self.logger.info("Using simulated reference genome coordinates")
    
    def _load_reference(self):
        """Load reference genome from FASTA file."""
        # In production, use pysam or BioPython to load FASTA
        self.logger.info(f"Loading reference from {self.reference_path}")
        pass
    
    def get_gene_info(self, gene_name: str) -> Optional[Dict]:
        """
        Get coordinate information for a gene.
        
        Parameters
        ----------
        gene_name : str
            Gene name (case-insensitive)
        
        Returns
        -------
        dict or None
            Gene coordinates and information
        """
        return self.GENE_COORDINATES.get(gene_name.lower())
    
    def position_in_gene(self, chrom: str, position: int) -> Optional[str]:
        """
        Check if a genomic position falls within a known gene.
        
        Parameters
        ----------
        chrom : str
            Chromosome name
        position : int
            Genomic position
        
        Returns
        -------
        str or None
            Gene name if position is within a gene
        """
        for gene_name, info in self.GENE_COORDINATES.items():
            if info['chromosome'] == chrom:
                if info['start'] <= position <= info['end']:
                    return gene_name
        return None


# =============================================================================
# ALIGNMENT ENGINE
# =============================================================================

class ReadAligner:
    """
    Align sequencing reads to the P. falciparum reference genome.
    
    This class provides methods for simulating and performing read alignments.
    For demonstration purposes, it generates realistic simulated alignments.
    In production, it interfaces with BWA or other aligners.
    
    Parameters
    ----------
    reference : ReferenceGenome
        Reference genome object
    
    Example
    -------
    >>> ref = ReferenceGenome()
    >>> aligner = ReadAligner(ref)
    >>> stats = aligner.align_reads('filtered.fastq', 'aligned.sam')
    """
    
    def __init__(self, reference: ReferenceGenome):
        """Initialize read aligner."""
        self.reference = reference
        self.logger = setup_logging()
    
    def _simulate_alignment(
        self,
        read_id: str,
        sequence: str,
        quality: str
    ) -> AlignmentRecord:
        """
        Simulate a realistic alignment for demonstration.
        
        This generates alignments that target known drug resistance genes
        with realistic mapping quality and CIGAR strings.
        """
        # Probability of mapping
        maps = random.random() < 0.95  # 95% mapping rate
        
        if not maps:
            return AlignmentRecord(
                query_name=read_id,
                flag=4,  # Unmapped
                reference_name='*',
                position=0,
                mapq=0,
                cigar='*',
                sequence=sequence,
                quality=quality
            )
        
        # Select a chromosome (weighted toward drug resistance gene chromosomes)
        gene_chroms = [
            ('Pf3D7_07_v3', 'pfcrt'),
            ('Pf3D7_05_v3', 'pfmdr1'),
            ('Pf3D7_13_v3', 'kelch13'),
            ('Pf3D7_04_v3', 'dhfr'),
            ('Pf3D7_08_v3', 'dhps')
        ]
        
        # 50% chance of mapping to a drug resistance gene region
        if random.random() < 0.5:
            chrom, gene = random.choice(gene_chroms)
            gene_info = self.reference.GENE_COORDINATES[gene]
            # Position within or near the gene
            position = random.randint(
                max(1, gene_info['start'] - 1000),
                gene_info['end'] + 1000
            )
        else:
            # Random chromosome
            chrom = random.choice(list(self.reference.chromosomes.keys()))
            max_pos = self.reference.chromosomes[chrom] - len(sequence)
            position = random.randint(1, max(1, max_pos))
        
        # Generate CIGAR string (mostly perfect matches with occasional variants)
        read_len = len(sequence)
        if random.random() < 0.1:  # 10% have small indels
            if random.random() < 0.5:
                # Insertion
                ins_pos = random.randint(10, read_len - 10)
                cigar = f"{ins_pos}M1I{read_len - ins_pos - 1}M"
            else:
                # Deletion
                del_pos = random.randint(10, read_len - 10)
                cigar = f"{del_pos}M1D{read_len - del_pos}M"
        else:
            cigar = f"{read_len}M"
        
        # Mapping quality (higher for unique mappings)
        mapq = random.randint(30, 60)
        
        # Flag (0 for primary forward, 16 for reverse)
        flag = 16 if random.random() < 0.5 else 0
        
        return AlignmentRecord(
            query_name=read_id,
            flag=flag,
            reference_name=chrom,
            position=position,
            mapq=mapq,
            cigar=cigar,
            sequence=sequence if flag == 0 else reverse_complement(sequence),
            quality=quality,
            tags={'RG': 'Z:sample1', 'NM': f'i:{random.randint(0, 3)}'}
        )
    
    def align_reads(
        self,
        fastq_path: str,
        output_path: str,
        threads: int = 4
    ) -> AlignmentStats:
        """
        Align reads from FASTQ file to reference genome.
        
        Parameters
        ----------
        fastq_path : str
            Path to input FASTQ file
        output_path : str
            Path for output SAM file
        threads : int
            Number of threads (for real aligners)
        
        Returns
        -------
        AlignmentStats
            Alignment statistics
        
        Example
        -------
        >>> stats = aligner.align_reads('filtered.fastq', 'aligned.sam')
        >>> print(f"Mapping rate: {stats.mapping_rate}%")
        """
        self.logger.info(f"Starting alignment: {fastq_path} -> {output_path}")
        
        ensure_directory(os.path.dirname(output_path) or '.')
        
        stats = AlignmentStats()
        mapq_sum = 0
        chrom_counts = {}
        
        # Generate SAM header
        header_lines = [
            '@HD\tVN:1.6\tSO:coordinate',
            f'@PG\tID:genomic_pipeline\tPN:ReadAligner\tVN:1.0'
        ]
        
        # Add sequence dictionary
        for chrom, length in self.reference.chromosomes.items():
            header_lines.append(f'@SQ\tSN:{chrom}\tLN:{length}')
        
        # Read FASTQ and generate alignments
        with open(output_path, 'w') as out:
            # Write header
            out.write('\n'.join(header_lines) + '\n')
            
            # Process reads
            with open(fastq_path, 'r') as f:
                while True:
                    header = f.readline().strip()
                    if not header:
                        break
                    
                    sequence = f.readline().strip()
                    _ = f.readline()  # Plus line
                    quality = f.readline().strip()
                    
                    read_id = header.split()[0].lstrip('@')
                    
                    # Generate alignment
                    alignment = self._simulate_alignment(read_id, sequence, quality)
                    
                    # Update statistics
                    stats.total_reads += 1
                    
                    if alignment.is_mapped:
                        stats.mapped_reads += 1
                        mapq_sum += alignment.mapq
                        
                        if alignment.mapq >= 30:
                            stats.unique_mapped += 1
                        else:
                            stats.multi_mapped += 1
                        
                        chrom_counts[alignment.reference_name] = \
                            chrom_counts.get(alignment.reference_name, 0) + 1
                    else:
                        stats.unmapped_reads += 1
                    
                    # Write SAM record
                    out.write(alignment.to_sam_line() + '\n')
        
        # Calculate final statistics
        if stats.mapped_reads > 0:
            stats.mean_mapq = mapq_sum / stats.mapped_reads
        
        # Calculate per-chromosome coverage
        for chrom, count in chrom_counts.items():
            stats.coverage_by_chrom[chrom] = count
        
        self.logger.info(
            f"Alignment complete: {stats.mapped_reads}/{stats.total_reads} "
            f"reads mapped ({stats.mapping_rate:.1f}%)"
        )
        
        return stats
    
    def sam_to_bam(
        self,
        sam_path: str,
        bam_path: str,
        sort: bool = True,
        index: bool = True
    ) -> str:
        """
        Convert SAM to BAM format, optionally sort and index.
        
        Parameters
        ----------
        sam_path : str
            Input SAM file path
        bam_path : str
            Output BAM file path
        sort : bool
            Sort the BAM file (default: True)
        index : bool
            Create BAM index (default: True)
        
        Returns
        -------
        str
            Path to final BAM file
        
        Note
        ----
        In production, this uses pysam or samtools.
        For demonstration, we create a placeholder.
        """
        self.logger.info(f"Converting SAM to BAM: {sam_path} -> {bam_path}")
        
        # In production, use pysam:
        # import pysam
        # pysam.view('-bS', sam_path, '-o', bam_path, save_stdout=bam_path)
        # if sort:
        #     sorted_bam = bam_path.replace('.bam', '.sorted.bam')
        #     pysam.sort('-o', sorted_bam, bam_path)
        #     bam_path = sorted_bam
        # if index:
        #     pysam.index(bam_path)
        
        # For demonstration, create a placeholder BAM description
        with open(bam_path, 'w') as f:
            f.write(f"# Placeholder BAM file\n")
            f.write(f"# Source: {sam_path}\n")
            f.write(f"# Sorted: {sort}\n")
            f.write(f"# Indexed: {index}\n")
            f.write(f"# Note: In production, use pysam or samtools\n")
        
        self.logger.info("BAM conversion complete (simulated)")
        
        return bam_path


# =============================================================================
# CONVENIENCE FUNCTIONS
# =============================================================================

def align_fastq_to_reference(
    fastq_path: str,
    output_dir: str,
    sample_name: Optional[str] = None
) -> Tuple[str, AlignmentStats]:
    """
    Convenience function to align a FASTQ file to reference.
    
    Parameters
    ----------
    fastq_path : str
        Path to input FASTQ file
    output_dir : str
        Output directory for alignment files
    sample_name : str, optional
        Sample name (derived from filename if not provided)
    
    Returns
    -------
    tuple
        (path to SAM file, AlignmentStats)
    
    Example
    -------
    >>> sam_path, stats = align_fastq_to_reference(
    ...     'filtered.fastq',
    ...     'aligned/'
    ... )
    """
    ensure_directory(output_dir)
    
    if sample_name is None:
        sample_name = os.path.splitext(os.path.basename(fastq_path))[0]
    
    sam_path = os.path.join(output_dir, f"{sample_name}.sam")
    
    # Initialize reference and aligner
    reference = ReferenceGenome()
    aligner = ReadAligner(reference)
    
    # Perform alignment
    stats = aligner.align_reads(fastq_path, sam_path)
    
    return sam_path, stats


if __name__ == '__main__':
    print("Read Alignment Module for P. falciparum")
    print("-" * 50)
    
    # Show gene coordinates
    ref = ReferenceGenome()
    print("\nDrug Resistance Gene Coordinates:")
    for gene, info in ref.GENE_COORDINATES.items():
        print(f"  {gene}: {info['chromosome']}:{info['start']}-{info['end']}")
