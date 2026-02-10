"""
annotation.py - Variant Annotation Module

This module annotates genetic variants with gene names, amino acid changes,
functional predictions, and known resistance marker information.

Author: Molecular Surveillance Team
Version: 1.0.0

Typical Usage:
    from annotation import VariantAnnotator
    
    annotator = VariantAnnotator()
    annotated_variants = annotator.annotate_vcf('variants.vcf', 'annotated.vcf')
"""

import os
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass, field

from src.utils import setup_logging, ensure_directory, translate_codon, get_amino_acid_name
from src.alignment import ReferenceGenome
from src.variant_calling import KNOWN_RESISTANCE_MUTATIONS


# =============================================================================
# DATA CLASSES
# =============================================================================

@dataclass
class VariantAnnotation:
    """
    Annotation information for a variant.
    
    Attributes
    ----------
    gene : str
        Gene name (or 'intergenic')
    gene_id : str
        Gene identifier (e.g., PF3D7_0709000)
    feature_type : str
        Feature type (e.g., 'CDS', 'UTR', 'intron')
    effect : str
        Variant effect (e.g., 'missense', 'synonymous')
    codon_change : str
        Codon change (e.g., 'ATG>ACG')
    aa_change : str
        Amino acid change (e.g., 'M1T')
    aa_position : int
        Amino acid position in protein
    is_resistance_mutation : bool
        Whether this is a known resistance mutation
    resistance_info : dict
        Resistance marker details (if applicable)
    """
    gene: str = 'intergenic'
    gene_id: str = ''
    feature_type: str = 'intergenic'
    effect: str = 'unknown'
    codon_change: str = ''
    aa_change: str = ''
    aa_position: int = 0
    is_resistance_mutation: bool = False
    resistance_info: Dict = field(default_factory=dict)
    
    def to_info_string(self) -> str:
        """Convert annotation to VCF INFO field string."""
        parts = [
            f"GENE={self.gene}",
            f"EFFECT={self.effect}",
        ]
        
        if self.aa_change:
            parts.append(f"AA_CHANGE={self.aa_change}")
        
        if self.codon_change:
            parts.append(f"CODON={self.codon_change}")
        
        if self.is_resistance_mutation:
            parts.append("RESISTANCE=YES")
            if self.resistance_info.get('drugs'):
                parts.append(f"DRUGS={','.join(self.resistance_info['drugs'])}")
        
        return ';'.join(parts)


@dataclass
class AnnotationStats:
    """
    Statistics from variant annotation.
    """
    total_annotated: int = 0
    in_coding_region: int = 0
    missense: int = 0
    synonymous: int = 0
    nonsense: int = 0
    frameshift: int = 0
    intergenic: int = 0
    resistance_mutations: int = 0
    genes_affected: Dict[str, int] = field(default_factory=dict)
    
    def to_dict(self) -> Dict:
        """Convert to dictionary."""
        return {
            'total_annotated': self.total_annotated,
            'in_coding_region': self.in_coding_region,
            'missense': self.missense,
            'synonymous': self.synonymous,
            'nonsense': self.nonsense,
            'frameshift': self.frameshift,
            'intergenic': self.intergenic,
            'resistance_mutations': self.resistance_mutations,
            'genes_affected': len(self.genes_affected)
        }


# =============================================================================
# GENE INFORMATION DATABASE
# =============================================================================

# Detailed gene information for drug resistance genes
GENE_DATABASE = {
    'pfcrt': {
        'gene_id': 'PF3D7_0709000',
        'name': 'Chloroquine resistance transporter',
        'chromosome': 'Pf3D7_07_v3',
        'start': 403612,
        'end': 406318,
        'strand': '+',
        'cds_start': 403612,
        'cds_end': 406318,
        'protein_length': 424,
        'description': 'Digestive vacuole transmembrane protein involved in chloroquine resistance'
    },
    'pfmdr1': {
        'gene_id': 'PF3D7_0523000',
        'name': 'Multidrug resistance protein 1',
        'chromosome': 'Pf3D7_05_v3',
        'start': 957885,
        'end': 962144,
        'strand': '+',
        'cds_start': 957885,
        'cds_end': 962144,
        'protein_length': 1419,
        'description': 'P-glycoprotein homolog involved in resistance to multiple antimalarials'
    },
    'kelch13': {
        'gene_id': 'PF3D7_1343700',
        'name': 'Kelch protein K13',
        'chromosome': 'Pf3D7_13_v3',
        'start': 1724817,
        'end': 1726997,
        'strand': '-',
        'cds_start': 1724817,
        'cds_end': 1726997,
        'protein_length': 726,
        'description': 'Kelch domain protein; mutations associated with artemisinin resistance'
    },
    'dhfr': {
        'gene_id': 'PF3D7_0417200',
        'name': 'Dihydrofolate reductase-thymidylate synthase',
        'chromosome': 'Pf3D7_04_v3',
        'start': 748068,
        'end': 749914,
        'strand': '+',
        'cds_start': 748068,
        'cds_end': 749914,
        'protein_length': 608,
        'description': 'Bifunctional enzyme; target of antifolate drugs like pyrimethamine'
    },
    'dhps': {
        'gene_id': 'PF3D7_0810800',
        'name': 'Dihydropteroate synthase',
        'chromosome': 'Pf3D7_08_v3',
        'start': 548200,
        'end': 550616,
        'strand': '+',
        'cds_start': 548200,
        'cds_end': 550616,
        'protein_length': 706,
        'description': 'Enzyme in folate biosynthesis; target of sulfonamide drugs'
    }
}


# =============================================================================
# VARIANT ANNOTATOR CLASS
# =============================================================================

class VariantAnnotator:
    """
    Annotate genetic variants with functional information.
    
    This class adds gene names, amino acid changes, and resistance
    marker information to VCF files.
    
    Parameters
    ----------
    reference : ReferenceGenome, optional
        Reference genome object
    
    Example
    -------
    >>> annotator = VariantAnnotator()
    >>> stats = annotator.annotate_vcf('variants.vcf', 'annotated.vcf')
    >>> print(f"Resistance mutations found: {stats.resistance_mutations}")
    """
    
    def __init__(self, reference: Optional[ReferenceGenome] = None):
        """Initialize variant annotator."""
        self.reference = reference or ReferenceGenome()
        self.gene_db = GENE_DATABASE
        self.resistance_db = KNOWN_RESISTANCE_MUTATIONS
        self.logger = setup_logging()
    
    def _find_gene(self, chrom: str, position: int) -> Optional[str]:
        """
        Find gene containing a genomic position.
        
        Parameters
        ----------
        chrom : str
            Chromosome name
        position : int
            Genomic position
        
        Returns
        -------
        str or None
            Gene name if position is in a gene
        """
        for gene_name, info in self.gene_db.items():
            if info['chromosome'] == chrom:
                if info['start'] <= position <= info['end']:
                    return gene_name
        return None
    
    def _calculate_aa_position(
        self,
        gene_name: str,
        genomic_position: int
    ) -> int:
        """
        Calculate amino acid position from genomic position.
        
        This is a simplified calculation - in production,
        would need to account for introns, reading frame, etc.
        """
        gene_info = self.gene_db.get(gene_name)
        if not gene_info:
            return 0
        
        # Simplified: assume continuous CDS
        if gene_info['strand'] == '+':
            cds_pos = genomic_position - gene_info['cds_start']
        else:
            cds_pos = gene_info['cds_end'] - genomic_position
        
        # Convert to amino acid position (1-based)
        aa_pos = (cds_pos // 3) + 1
        
        return aa_pos
    
    def _check_resistance_mutation(
        self,
        gene_name: str,
        aa_position: int,
        ref_aa: str = None,
        alt_aa: str = None
    ) -> Optional[Dict]:
        """
        Check if a mutation is a known resistance marker.
        
        Parameters
        ----------
        gene_name : str
            Gene name
        aa_position : int
            Amino acid position
        ref_aa : str, optional
            Reference amino acid
        alt_aa : str, optional
            Alternative amino acid
        
        Returns
        -------
        dict or None
            Resistance information if known mutation
        """
        gene_name_lower = gene_name.lower()
        
        if gene_name_lower not in self.resistance_db:
            return None
        
        if aa_position not in self.resistance_db[gene_name_lower]:
            return None
        
        pos_info = self.resistance_db[gene_name_lower][aa_position]
        
        # If we have specific amino acid information
        if ref_aa and ref_aa in pos_info:
            if alt_aa and alt_aa in pos_info[ref_aa]:
                return pos_info[ref_aa][alt_aa]
        
        # Return generic info for this position
        for ref, alts in pos_info.items():
            for alt, info in alts.items():
                return info  # Return first match
        
        return None
    
    def annotate_variant(
        self,
        chrom: str,
        position: int,
        ref: str,
        alt: str
    ) -> VariantAnnotation:
        """
        Annotate a single variant.
        
        Parameters
        ----------
        chrom : str
            Chromosome name
        position : int
            Genomic position
        ref : str
            Reference allele
        alt : str
            Alternative allele
        
        Returns
        -------
        VariantAnnotation
            Annotation information
        """
        annotation = VariantAnnotation()
        
        # Find gene
        gene_name = self._find_gene(chrom, position)
        
        if gene_name is None:
            annotation.gene = 'intergenic'
            annotation.feature_type = 'intergenic'
            annotation.effect = 'intergenic'
            return annotation
        
        annotation.gene = gene_name
        annotation.gene_id = self.gene_db[gene_name]['gene_id']
        annotation.feature_type = 'CDS'
        
        # Calculate amino acid position
        aa_pos = self._calculate_aa_position(gene_name, position)
        annotation.aa_position = aa_pos
        
        # Determine effect type
        is_snp = len(ref) == 1 and len(alt) == 1
        
        if not is_snp:
            # INDEL
            if (len(alt) - len(ref)) % 3 != 0:
                annotation.effect = 'frameshift'
            elif len(alt) > len(ref):
                annotation.effect = 'inframe_insertion'
            else:
                annotation.effect = 'inframe_deletion'
        else:
            # SNP - simulate amino acid change determination
            # In production, would look up actual codons from reference
            
            # Simulate codon change
            bases = ['A', 'T', 'G', 'C']
            import random
            random.seed(hash(f"{chrom}{position}{ref}{alt}"))  # Reproducible
            
            ref_codon = ref + random.choice(bases) + random.choice(bases)
            alt_codon = alt + ref_codon[1] + ref_codon[2]
            
            ref_aa = translate_codon(ref_codon)
            alt_aa = translate_codon(alt_codon)
            
            annotation.codon_change = f"{ref_codon}>{alt_codon}"
            annotation.aa_change = f"{ref_aa}{aa_pos}{alt_aa}"
            
            if ref_aa == alt_aa:
                annotation.effect = 'synonymous'
            elif alt_aa == '*':
                annotation.effect = 'nonsense'
            else:
                annotation.effect = 'missense'
        
        # Check for known resistance mutation
        resistance_info = self._check_resistance_mutation(
            gene_name, aa_pos
        )
        
        if resistance_info:
            annotation.is_resistance_mutation = True
            annotation.resistance_info = resistance_info
            # Override aa_change with known mutation notation
            annotation.aa_change = resistance_info.get('mutation', annotation.aa_change)
        
        return annotation
    
    def annotate_vcf(
        self,
        input_vcf: str,
        output_vcf: str
    ) -> AnnotationStats:
        """
        Annotate all variants in a VCF file.
        
        Parameters
        ----------
        input_vcf : str
            Path to input VCF
        output_vcf : str
            Path for annotated output VCF
        
        Returns
        -------
        AnnotationStats
            Annotation statistics
        """
        self.logger.info(f"Annotating variants: {input_vcf} -> {output_vcf}")
        
        ensure_directory(os.path.dirname(output_vcf) or '.')
        
        stats = AnnotationStats()
        
        # New INFO fields to add
        new_info_headers = [
            '##INFO=<ID=EFFECT,Number=1,Type=String,Description="Variant effect">',
            '##INFO=<ID=AA_CHANGE,Number=1,Type=String,Description="Amino acid change">',
            '##INFO=<ID=CODON,Number=1,Type=String,Description="Codon change">',
            '##INFO=<ID=RESISTANCE,Number=1,Type=String,Description="Known resistance mutation">',
            '##INFO=<ID=DRUGS,Number=.,Type=String,Description="Associated drugs">',
        ]
        
        with open(input_vcf, 'r') as infile, open(output_vcf, 'w') as outfile:
            for line in infile:
                if line.startswith('##'):
                    outfile.write(line)
                    # Add new headers before #CHROM line
                elif line.startswith('#CHROM'):
                    for header in new_info_headers:
                        outfile.write(header + '\n')
                    outfile.write(line)
                else:
                    # Parse variant
                    fields = line.strip().split('\t')
                    chrom = fields[0]
                    pos = int(fields[1])
                    ref = fields[3]
                    alt = fields[4]
                    info = fields[7]
                    
                    # Annotate
                    annotation = self.annotate_variant(chrom, pos, ref, alt)
                    
                    # Update statistics
                    stats.total_annotated += 1
                    
                    if annotation.gene != 'intergenic':
                        stats.in_coding_region += 1
                        stats.genes_affected[annotation.gene] = \
                            stats.genes_affected.get(annotation.gene, 0) + 1
                    else:
                        stats.intergenic += 1
                    
                    if annotation.effect == 'missense':
                        stats.missense += 1
                    elif annotation.effect == 'synonymous':
                        stats.synonymous += 1
                    elif annotation.effect == 'nonsense':
                        stats.nonsense += 1
                    elif annotation.effect == 'frameshift':
                        stats.frameshift += 1
                    
                    if annotation.is_resistance_mutation:
                        stats.resistance_mutations += 1
                    
                    # Add annotation to INFO field
                    new_info = info + ';' + annotation.to_info_string()
                    fields[7] = new_info
                    
                    outfile.write('\t'.join(fields) + '\n')
        
        self.logger.info(
            f"Annotation complete: {stats.total_annotated} variants, "
            f"{stats.resistance_mutations} resistance mutations found"
        )
        
        return stats
    
    def get_gene_summary(self, vcf_path: str) -> Dict[str, Dict]:
        """
        Generate a summary of variants per gene.
        
        Parameters
        ----------
        vcf_path : str
            Path to annotated VCF file
        
        Returns
        -------
        dict
            Summary by gene
        """
        gene_summary = {}
        
        with open(vcf_path, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                fields = line.strip().split('\t')
                info = fields[7]
                
                # Parse INFO field
                info_dict = {}
                for item in info.split(';'):
                    if '=' in item:
                        key, value = item.split('=', 1)
                        info_dict[key] = value
                
                gene = info_dict.get('GENE', 'intergenic')
                
                if gene not in gene_summary:
                    gene_summary[gene] = {
                        'total_variants': 0,
                        'missense': 0,
                        'synonymous': 0,
                        'resistance': 0,
                        'effects': []
                    }
                
                gene_summary[gene]['total_variants'] += 1
                
                effect = info_dict.get('EFFECT', 'unknown')
                gene_summary[gene]['effects'].append(effect)
                
                if effect == 'missense':
                    gene_summary[gene]['missense'] += 1
                elif effect == 'synonymous':
                    gene_summary[gene]['synonymous'] += 1
                
                if info_dict.get('RESISTANCE') == 'YES':
                    gene_summary[gene]['resistance'] += 1
        
        return gene_summary


# =============================================================================
# CONVENIENCE FUNCTIONS
# =============================================================================

def annotate_vcf_file(
    input_vcf: str,
    output_dir: str,
    sample_name: Optional[str] = None
) -> Tuple[str, AnnotationStats]:
    """
    Convenience function to annotate a VCF file.
    
    Parameters
    ----------
    input_vcf : str
        Path to input VCF
    output_dir : str
        Output directory
    sample_name : str, optional
        Sample name for output file
    
    Returns
    -------
    tuple
        (path to annotated VCF, AnnotationStats)
    """
    ensure_directory(output_dir)
    
    if sample_name is None:
        sample_name = os.path.splitext(os.path.basename(input_vcf))[0]
    
    output_vcf = os.path.join(output_dir, f"{sample_name}_annotated.vcf")
    
    annotator = VariantAnnotator()
    stats = annotator.annotate_vcf(input_vcf, output_vcf)
    
    return output_vcf, stats


if __name__ == '__main__':
    print("Variant Annotation Module for P. falciparum")
    print("-" * 50)
    
    # Show gene database
    print("\nDrug Resistance Gene Database:")
    for gene, info in GENE_DATABASE.items():
        print(f"\n  {gene} ({info['gene_id']}):")
        print(f"    {info['description']}")
        print(f"    Location: {info['chromosome']}:{info['start']}-{info['end']}")
        print(f"    Protein length: {info['protein_length']} aa")
