"""
variant_calling.py - Variant Calling Module

This module detects single nucleotide polymorphisms (SNPs) and insertions/deletions
(INDELs) from aligned sequencing data. It generates VCF files following standard
specifications.

Author: Molecular Surveillance Team
Version: 1.0.0

Note: This module provides simulated variant calling for demonstration.
      In production, use GATK, FreeBayes, or bcftools for real variant calling.
"""

import os
import random
from typing import Dict, List, Optional, Tuple, Set
from dataclasses import dataclass, field
from datetime import datetime

from src.utils import setup_logging, ensure_directory
from src.alignment import ReferenceGenome


# =============================================================================
# DATA CLASSES
# =============================================================================

@dataclass
class Variant:
    """
    Represents a genetic variant (SNP or INDEL).
    
    Attributes
    ----------
    chrom : str
        Chromosome name
    pos : int
        1-based position
    id : str
        Variant identifier (or '.')
    ref : str
        Reference allele
    alt : str
        Alternative allele(s)
    qual : float
        Variant quality score
    filter : str
        Filter status (PASS or filter name)
    info : dict
        INFO field key-value pairs
    format_keys : list
        FORMAT field keys
    sample_data : dict
        Per-sample genotype data
    """
    chrom: str
    pos: int
    ref: str
    alt: str
    qual: float = 30.0
    id: str = '.'
    filter: str = 'PASS'
    info: Dict[str, str] = field(default_factory=dict)
    format_keys: List[str] = field(default_factory=lambda: ['GT', 'DP', 'AD'])
    sample_data: Dict[str, Dict[str, str]] = field(default_factory=dict)
    
    @property
    def is_snp(self) -> bool:
        """Check if variant is a SNP."""
        return len(self.ref) == 1 and len(self.alt) == 1
    
    @property
    def is_indel(self) -> bool:
        """Check if variant is an INDEL."""
        return len(self.ref) != len(self.alt)
    
    @property
    def variant_type(self) -> str:
        """Get variant type string."""
        if self.is_snp:
            return 'SNP'
        elif len(self.ref) < len(self.alt):
            return 'INS'
        else:
            return 'DEL'
    
    def to_vcf_line(self, sample_names: List[str] = None) -> str:
        """Convert to VCF format line."""
        if sample_names is None:
            sample_names = list(self.sample_data.keys()) or ['SAMPLE']
        
        # Build INFO string
        info_str = ';'.join(f"{k}={v}" for k, v in self.info.items()) or '.'
        
        # Build FORMAT string
        format_str = ':'.join(self.format_keys)
        
        # Build sample fields
        sample_fields = []
        for sample in sample_names:
            if sample in self.sample_data:
                values = [self.sample_data[sample].get(k, '.') for k in self.format_keys]
                sample_fields.append(':'.join(values))
            else:
                sample_fields.append('./.')
        
        fields = [
            self.chrom,
            str(self.pos),
            self.id,
            self.ref,
            self.alt,
            f"{self.qual:.1f}",
            self.filter,
            info_str,
            format_str
        ] + sample_fields
        
        return '\t'.join(fields)


@dataclass
class VariantCallStats:
    """
    Statistics from variant calling.
    
    Attributes
    ----------
    total_variants : int
        Total variants called
    snps : int
        Number of SNPs
    insertions : int
        Number of insertions
    deletions : int
        Number of deletions
    pass_filter : int
        Variants passing filters
    """
    total_variants: int = 0
    snps: int = 0
    insertions: int = 0
    deletions: int = 0
    pass_filter: int = 0
    low_quality: int = 0
    low_depth: int = 0
    variants_per_chrom: Dict[str, int] = field(default_factory=dict)
    
    def to_dict(self) -> Dict:
        """Convert to dictionary."""
        return {
            'total_variants': self.total_variants,
            'snps': self.snps,
            'insertions': self.insertions,
            'deletions': self.deletions,
            'pass_filter': self.pass_filter,
            'low_quality': self.low_quality,
            'low_depth': self.low_depth,
            'ti_tv_ratio': self.ti_tv_ratio
        }
    
    @property
    def ti_tv_ratio(self) -> float:
        """Calculate transition/transversion ratio (typical: 2.0-2.1 for coding regions)."""
        # This would be calculated from actual variant data
        return 2.1


# =============================================================================
# KNOWN RESISTANCE MUTATIONS
# =============================================================================

# Known drug resistance mutations in P. falciparum
# Format: {gene: {position: {ref_aa: {alt_aa: {'mutation': str, 'drugs': list}}}}}
KNOWN_RESISTANCE_MUTATIONS = {
    'pfcrt': {
        72: {'C': {'S': {'mutation': 'C72S', 'drugs': ['Chloroquine']}}},
        74: {'M': {'I': {'mutation': 'M74I', 'drugs': ['Chloroquine']}}},
        75: {'N': {'E': {'mutation': 'N75E', 'drugs': ['Chloroquine']}}},
        76: {'K': {'T': {'mutation': 'K76T', 'drugs': ['Chloroquine', 'Amodiaquine']}}},
        220: {'A': {'S': {'mutation': 'A220S', 'drugs': ['Chloroquine']}}},
        271: {'Q': {'E': {'mutation': 'Q271E', 'drugs': ['Chloroquine']}}},
        326: {'N': {'S': {'mutation': 'N326S', 'drugs': ['Chloroquine']}}},
        356: {'I': {'T': {'mutation': 'I356T', 'drugs': ['Chloroquine']}}},
    },
    'pfmdr1': {
        86: {'N': {'Y': {'mutation': 'N86Y', 'drugs': ['Chloroquine', 'Amodiaquine', 'Lumefantrine']}}},
        184: {'Y': {'F': {'mutation': 'Y184F', 'drugs': ['Lumefantrine', 'Mefloquine']}}},
        1034: {'S': {'C': {'mutation': 'S1034C', 'drugs': ['Mefloquine']}}},
        1042: {'N': {'D': {'mutation': 'N1042D', 'drugs': ['Mefloquine']}}},
        1246: {'D': {'Y': {'mutation': 'D1246Y', 'drugs': ['Mefloquine', 'Lumefantrine']}}},
    },
    'kelch13': {
        # Validated artemisinin resistance mutations (WHO 2021)
        449: {'G': {'A': {'mutation': 'G449A', 'drugs': ['Artemisinin']}}},
        458: {'N': {'Y': {'mutation': 'N458Y', 'drugs': ['Artemisinin']}}},
        469: {'Y': {'C': {'mutation': 'Y469C', 'drugs': ['Artemisinin']}}},
        476: {'M': {'I': {'mutation': 'M476I', 'drugs': ['Artemisinin']}}},
        493: {'Y': {'H': {'mutation': 'Y493H', 'drugs': ['Artemisinin']}}},
        539: {'R': {'T': {'mutation': 'R539T', 'drugs': ['Artemisinin']}}},
        543: {'I': {'T': {'mutation': 'I543T', 'drugs': ['Artemisinin']}}},
        553: {'P': {'L': {'mutation': 'P553L', 'drugs': ['Artemisinin']}}},
        561: {'R': {'H': {'mutation': 'R561H', 'drugs': ['Artemisinin']}}},
        574: {'P': {'L': {'mutation': 'P574L', 'drugs': ['Artemisinin']}}},
        580: {'C': {'Y': {'mutation': 'C580Y', 'drugs': ['Artemisinin']}}},
        622: {'A': {'S': {'mutation': 'A622S', 'drugs': ['Artemisinin']}}},
        675: {'A': {'V': {'mutation': 'A675V', 'drugs': ['Artemisinin']}}},
    },
    'dhfr': {
        16: {'A': {'V': {'mutation': 'A16V', 'drugs': ['Pyrimethamine']}}},
        51: {'N': {'I': {'mutation': 'N51I', 'drugs': ['Pyrimethamine', 'Sulfadoxine-Pyrimethamine']}}},
        59: {'C': {'R': {'mutation': 'C59R', 'drugs': ['Pyrimethamine', 'Sulfadoxine-Pyrimethamine']}}},
        108: {'S': {'N': {'mutation': 'S108N', 'drugs': ['Pyrimethamine', 'Sulfadoxine-Pyrimethamine']}}},
        164: {'I': {'L': {'mutation': 'I164L', 'drugs': ['Pyrimethamine', 'Sulfadoxine-Pyrimethamine']}}},
    },
    'dhps': {
        436: {'S': {'A': {'mutation': 'S436A', 'drugs': ['Sulfadoxine']}}},
        437: {'A': {'G': {'mutation': 'A437G', 'drugs': ['Sulfadoxine', 'Sulfadoxine-Pyrimethamine']}}},
        540: {'K': {'E': {'mutation': 'K540E', 'drugs': ['Sulfadoxine', 'Sulfadoxine-Pyrimethamine']}}},
        581: {'A': {'G': {'mutation': 'A581G', 'drugs': ['Sulfadoxine', 'Sulfadoxine-Pyrimethamine']}}},
        613: {'A': {'S': {'mutation': 'A613S', 'drugs': ['Sulfadoxine']}}},
    }
}


# =============================================================================
# VARIANT CALLER CLASS
# =============================================================================

class VariantCaller:
    """
    Call variants from aligned sequencing data.
    
    This class identifies SNPs and INDELs from BAM/SAM alignments
    and generates VCF output files.
    
    Parameters
    ----------
    reference : ReferenceGenome
        Reference genome object
    min_depth : int
        Minimum read depth to call variant (default: 10)
    min_qual : float
        Minimum quality score (default: 20)
    
    Example
    -------
    >>> ref = ReferenceGenome()
    >>> caller = VariantCaller(ref, min_depth=10, min_qual=20)
    >>> stats = caller.call_variants('aligned.sam', 'variants.vcf')
    """
    
    def __init__(
        self,
        reference: ReferenceGenome,
        min_depth: int = 10,
        min_qual: float = 20.0
    ):
        """Initialize variant caller."""
        self.reference = reference
        self.min_depth = min_depth
        self.min_qual = min_qual
        self.logger = setup_logging()
    
    def _generate_simulated_variants(
        self,
        sample_name: str
    ) -> List[Variant]:
        """
        Generate realistic simulated variants for demonstration.
        
        This creates variants enriched in drug resistance genes
        to demonstrate the pipeline functionality.
        """
        variants = []
        
        # Generate variants for each drug resistance gene
        for gene, gene_info in self.reference.GENE_COORDINATES.items():
            chrom = gene_info['chromosome']
            start = gene_info['start']
            end = gene_info['end']
            
            # Generate 5-15 variants per gene
            n_variants = random.randint(5, 15)
            
            for _ in range(n_variants):
                # 30% chance of being a known resistance mutation
                if random.random() < 0.3 and gene in KNOWN_RESISTANCE_MUTATIONS:
                    # Pick a known mutation
                    positions = list(KNOWN_RESISTANCE_MUTATIONS[gene].keys())
                    codon_pos = random.choice(positions)
                    
                    # Calculate genomic position (simplified)
                    pos = start + (codon_pos * 3) + random.randint(0, 2)
                    
                    # Get mutation info
                    ref_aas = KNOWN_RESISTANCE_MUTATIONS[gene][codon_pos]
                    ref_aa = list(ref_aas.keys())[0]
                    alt_aas = ref_aas[ref_aa]
                    alt_aa = list(alt_aas.keys())[0]
                    
                    # Simplified - in reality would be codon-based
                    ref_base = random.choice(['A', 'T', 'G', 'C'])
                    alt_base = random.choice([b for b in ['A', 'T', 'G', 'C'] if b != ref_base])
                else:
                    # Random variant
                    pos = random.randint(start, end)
                    ref_base = random.choice(['A', 'T', 'G', 'C'])
                    alt_base = random.choice([b for b in ['A', 'T', 'G', 'C'] if b != ref_base])
                
                # Generate realistic depth and quality
                depth = random.randint(15, 200)
                ref_count = random.randint(0, depth // 3)
                alt_count = depth - ref_count
                qual = random.uniform(25, 60)
                
                # Determine filter status
                filter_status = 'PASS'
                if depth < self.min_depth:
                    filter_status = 'LowDepth'
                elif qual < self.min_qual:
                    filter_status = 'LowQual'
                
                # Determine genotype
                af = alt_count / depth
                if af > 0.85:
                    gt = '1/1'  # Homozygous alt
                elif af > 0.15:
                    gt = '0/1'  # Heterozygous
                else:
                    gt = '0/0'  # Homozygous ref
                
                variant = Variant(
                    chrom=chrom,
                    pos=pos,
                    ref=ref_base,
                    alt=alt_base,
                    qual=qual,
                    filter=filter_status,
                    info={
                        'DP': str(depth),
                        'AF': f"{af:.3f}",
                        'GENE': gene
                    },
                    format_keys=['GT', 'DP', 'AD', 'GQ'],
                    sample_data={
                        sample_name: {
                            'GT': gt,
                            'DP': str(depth),
                            'AD': f"{ref_count},{alt_count}",
                            'GQ': str(int(qual))
                        }
                    }
                )
                
                variants.append(variant)
        
        # Add some variants in other regions
        for chrom in list(self.reference.chromosomes.keys())[:5]:
            if any(g['chromosome'] == chrom for g in self.reference.GENE_COORDINATES.values()):
                continue
            
            n_variants = random.randint(3, 8)
            chrom_size = self.reference.chromosomes[chrom]
            
            for _ in range(n_variants):
                pos = random.randint(1000, chrom_size - 1000)
                ref_base = random.choice(['A', 'T', 'G', 'C'])
                alt_base = random.choice([b for b in ['A', 'T', 'G', 'C'] if b != ref_base])
                depth = random.randint(15, 150)
                qual = random.uniform(20, 55)
                
                variant = Variant(
                    chrom=chrom,
                    pos=pos,
                    ref=ref_base,
                    alt=alt_base,
                    qual=qual,
                    filter='PASS' if qual >= self.min_qual else 'LowQual',
                    info={'DP': str(depth)},
                    sample_data={
                        sample_name: {
                            'GT': random.choice(['0/1', '1/1']),
                            'DP': str(depth),
                            'AD': f"{random.randint(0, depth//2)},{random.randint(depth//2, depth)}",
                            'GQ': str(int(qual))
                        }
                    }
                )
                variants.append(variant)
        
        # Sort by chromosome and position
        variants.sort(key=lambda v: (v.chrom, v.pos))
        
        return variants
    
    def _generate_vcf_header(
        self,
        sample_names: List[str],
        source: str = "GenomicVariantPipeline"
    ) -> List[str]:
        """Generate VCF header lines."""
        today = datetime.now().strftime('%Y%m%d')
        
        header = [
            '##fileformat=VCFv4.2',
            f'##fileDate={today}',
            f'##source={source}',
            '##reference=Pf3D7_v3',
            
            # Contigs
        ]
        
        for chrom, length in self.reference.chromosomes.items():
            header.append(f'##contig=<ID={chrom},length={length}>')
        
        # INFO fields
        header.extend([
            '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">',
            '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">',
            '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name if in coding region">',
        ])
        
        # FILTER fields
        header.extend([
            '##FILTER=<ID=PASS,Description="All filters passed">',
            '##FILTER=<ID=LowQual,Description="Quality below threshold">',
            '##FILTER=<ID=LowDepth,Description="Read depth below threshold">',
        ])
        
        # FORMAT fields
        header.extend([
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
            '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
            '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">',
            '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">',
        ])
        
        # Column header
        columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
        columns.extend(sample_names)
        header.append('\t'.join(columns))
        
        return header
    
    def call_variants(
        self,
        alignment_path: str,
        output_path: str,
        sample_name: Optional[str] = None
    ) -> VariantCallStats:
        """
        Call variants from alignment file.
        
        Parameters
        ----------
        alignment_path : str
            Path to input SAM/BAM file
        output_path : str
            Path for output VCF file
        sample_name : str, optional
            Sample name for VCF (derived from filename if not provided)
        
        Returns
        -------
        VariantCallStats
            Variant calling statistics
        
        Example
        -------
        >>> stats = caller.call_variants('aligned.sam', 'variants.vcf')
        >>> print(f"Total variants: {stats.total_variants}")
        """
        self.logger.info(f"Starting variant calling: {alignment_path} -> {output_path}")
        
        ensure_directory(os.path.dirname(output_path) or '.')
        
        if sample_name is None:
            sample_name = os.path.splitext(os.path.basename(alignment_path))[0]
        
        # Generate simulated variants
        variants = self._generate_simulated_variants(sample_name)
        
        # Calculate statistics
        stats = VariantCallStats()
        
        for variant in variants:
            stats.total_variants += 1
            
            if variant.is_snp:
                stats.snps += 1
            elif len(variant.ref) < len(variant.alt):
                stats.insertions += 1
            else:
                stats.deletions += 1
            
            if variant.filter == 'PASS':
                stats.pass_filter += 1
            elif variant.filter == 'LowQual':
                stats.low_quality += 1
            elif variant.filter == 'LowDepth':
                stats.low_depth += 1
            
            stats.variants_per_chrom[variant.chrom] = \
                stats.variants_per_chrom.get(variant.chrom, 0) + 1
        
        # Write VCF file
        with open(output_path, 'w') as f:
            # Header
            header = self._generate_vcf_header([sample_name])
            f.write('\n'.join(header) + '\n')
            
            # Variants
            for variant in variants:
                f.write(variant.to_vcf_line([sample_name]) + '\n')
        
        self.logger.info(
            f"Variant calling complete: {stats.total_variants} variants "
            f"({stats.snps} SNPs, {stats.insertions} INS, {stats.deletions} DEL)"
        )
        
        return stats
    
    def filter_vcf(
        self,
        input_vcf: str,
        output_vcf: str,
        min_qual: float = None,
        min_depth: int = None,
        pass_only: bool = True
    ) -> Tuple[int, int]:
        """
        Filter VCF file based on quality criteria.
        
        Parameters
        ----------
        input_vcf : str
            Path to input VCF
        output_vcf : str
            Path for filtered VCF
        min_qual : float, optional
            Minimum quality (uses class default if not specified)
        min_depth : int, optional
            Minimum depth (uses class default if not specified)
        pass_only : bool
            Keep only PASS variants (default: True)
        
        Returns
        -------
        tuple
            (total_variants, passed_variants)
        """
        min_qual = min_qual or self.min_qual
        min_depth = min_depth or self.min_depth
        
        total = 0
        passed = 0
        
        with open(input_vcf, 'r') as infile, open(output_vcf, 'w') as outfile:
            for line in infile:
                if line.startswith('#'):
                    outfile.write(line)
                    continue
                
                total += 1
                fields = line.strip().split('\t')
                
                # Check filter
                if pass_only and fields[6] != 'PASS':
                    continue
                
                # Check quality
                qual = float(fields[5])
                if qual < min_qual:
                    continue
                
                # Check depth from INFO
                info_dict = dict(item.split('=') for item in fields[7].split(';') if '=' in item)
                depth = int(info_dict.get('DP', 0))
                if depth < min_depth:
                    continue
                
                outfile.write(line)
                passed += 1
        
        self.logger.info(f"Filtered VCF: {passed}/{total} variants retained")
        
        return total, passed


# =============================================================================
# CONVENIENCE FUNCTIONS
# =============================================================================

def call_variants_from_alignment(
    alignment_path: str,
    output_dir: str,
    sample_name: Optional[str] = None,
    min_depth: int = 10,
    min_qual: float = 20.0
) -> Tuple[str, VariantCallStats]:
    """
    Convenience function to call variants from an alignment file.
    
    Parameters
    ----------
    alignment_path : str
        Path to SAM/BAM file
    output_dir : str
        Output directory for VCF file
    sample_name : str, optional
        Sample name
    min_depth : int
        Minimum read depth
    min_qual : float
        Minimum quality score
    
    Returns
    -------
    tuple
        (path to VCF file, VariantCallStats)
    """
    ensure_directory(output_dir)
    
    if sample_name is None:
        sample_name = os.path.splitext(os.path.basename(alignment_path))[0]
        sample_name = sample_name.replace('.sorted', '').replace('.filtered', '')
    
    vcf_path = os.path.join(output_dir, f"{sample_name}.vcf")
    
    reference = ReferenceGenome()
    caller = VariantCaller(reference, min_depth=min_depth, min_qual=min_qual)
    
    stats = caller.call_variants(alignment_path, vcf_path, sample_name)
    
    return vcf_path, stats


if __name__ == '__main__':
    print("Variant Calling Module for P. falciparum")
    print("-" * 50)
    
    # Show known resistance mutations
    print("\nKnown Drug Resistance Mutations:")
    for gene, positions in KNOWN_RESISTANCE_MUTATIONS.items():
        print(f"\n  {gene}:")
        for pos, info in list(positions.items())[:3]:
            ref_aa = list(info.keys())[0]
            alt_info = info[ref_aa]
            alt_aa = list(alt_info.keys())[0]
            mut = alt_info[alt_aa]
            print(f"    {mut['mutation']}: {', '.join(mut['drugs'])}")
